template <typename Scalar>
LocalSolution<Scalar> make_local_solution(
    const RecurrenceProblem<Scalar>& problem,
    AssembledResult<Scalar>&& assembled, LocalMetadata&& metadata) {
  LocalSolution<Scalar> solution;
  solution.chart = std::move(metadata.chart);
  solution.epsilon = {assembled.min_power, assembled.complete_max};
  solution.taylor_complete_max = problem.nmax;
  solution.dimension = problem.dimension;
  solution.prescriptions = std::move(metadata.prescriptions);
  solution.checkpoint_identity = std::move(metadata.checkpoint_identity);

  const auto sector_size = solution.sector_size();
  const auto sector_count = static_cast<std::size_t>(problem.log_max) + 1;
  if (sector_size > std::numeric_limits<std::size_t>::max() / sector_count ||
      assembled.coefficients.size() != sector_size * sector_count)
    throw std::invalid_argument(
        "assembled recurrence tensor cannot form the declared local sectors");
  auto cursor = assembled.coefficients.begin();
  for (std::uint32_t log = 0; log <= problem.log_max; ++log) {
    LocalSector<Scalar> sector;
    sector.a = metadata.a;
    sector.b = metadata.b;
    sector.log_power = log;
    sector.coefficients.reserve(sector_size);
    auto end = cursor + static_cast<std::ptrdiff_t>(sector_size);
    sector.coefficients.insert(
        sector.coefficients.end(), std::make_move_iterator(cursor),
        std::make_move_iterator(end));
    cursor = end;
    solution.sectors.push_back(std::move(sector));
  }
  return solution;
}

template <typename Scalar>
class RegularPhysicalEquationOwner final
    : public RegularPhysicalEquationOwnerBase {
 public:
  RegularPhysicalEquationOwner(
      std::string handle, std::string key, std::string exact_identity,
      std::string signature, std::string geometry_record,
      std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>> equation,
      slong precision_bits,
      std::string regular_value_relative_accuracy_max_exact)
      : handle_(std::move(handle)), key_(std::move(key)),
        exact_identity_(std::move(exact_identity)),
        signature_(std::move(signature)),
        geometry_record_(std::move(geometry_record)),
        equation_(std::move(equation)), precision_bits_(precision_bits),
        regular_value_relative_accuracy_max_exact_(
            std::move(regular_value_relative_accuracy_max_exact)) {
    if (handle_.empty() || key_.empty() || exact_identity_.empty() ||
        signature_.empty() || geometry_record_.empty() || !equation_ ||
        precision_bits_ < 64 ||
        regular_value_relative_accuracy_max_exact_.empty() ||
        equation_->owner_signature_identity != exact_identity_)
      throw std::invalid_argument(
          "regular physical equation owner lost an exact identity, geometry, or q/C payload");
  }

  const std::string& handle() const override { return handle_; }
  const std::string& key() const override { return key_; }
  const std::string& exact_identity() const override {
    return exact_identity_;
  }
  const std::string& signature() const override { return signature_; }
  const std::string& geometry_record() const override {
    return geometry_record_;
  }
  std::uint32_t dimension() const override { return equation_->dimension; }
  slong precision_bits() const override { return precision_bits_; }
  const std::string& regular_value_relative_accuracy_max_exact()
      const override {
    return regular_value_relative_accuracy_max_exact_;
  }
  const char* equation_scalar_domain() const override {
    if constexpr (std::is_same_v<Scalar, Rational>) return "rational";
    if constexpr (std::is_same_v<Scalar, ComplexBall>) return "acb";
    return "symbolic";
  }
  std::shared_ptr<const void> physical_ode_erased() const override {
    return std::static_pointer_cast<const void>(equation_);
  }
  const std::string& physical_payload_identity() const override {
    return equation_->payload_identity;
  }
  const std::string& physical_payload_record() const override {
    return equation_->exact_payload_record;
  }
  const std::string& owner_signature_identity() const override {
    return equation_->owner_signature_identity;
  }
  const std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>&
  physical_equation() const {
    return equation_;
  }

 private:
  std::string handle_;
  std::string key_;
  std::string exact_identity_;
  std::string signature_;
  std::string geometry_record_;
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>> equation_;
  slong precision_bits_ = 0;
  std::string regular_value_relative_accuracy_max_exact_;
};

template <typename Scalar>
class PreparedChart final : public PreparedChartBase {
 public:
  PreparedChart(std::string handle, std::string key,
                std::string exact_identity, std::string signature,
                std::optional<std::string> geometry_record,
                std::optional<std::string> principal_matrix_record,
                std::optional<std::string> native_scc_capabilities,
                std::optional<std::string>
                    regular_value_relative_accuracy_max_exact,
                SCCCertificate scc,
                PreparedRecurrenceOperator<Scalar>&& prepared,
                std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
                    physical_equation,
                slong precision_bits, std::vector<std::string> symbols,
                double prepare_parse_ms)
      : PreparedChartBase(std::move(handle), std::move(key),
                          std::move(exact_identity), std::move(signature),
                          std::move(geometry_record),
                          std::move(principal_matrix_record),
                          std::move(native_scc_capabilities),
                          std::move(regular_value_relative_accuracy_max_exact),
                          std::move(scc),
                          prepare_parse_ms),
        prepared_(std::move(prepared)),
        physical_equation_(std::move(physical_equation)),
        precision_bits_(precision_bits),
        symbols_(std::move(symbols)) {
    if (physical_equation_ &&
        physical_equation_->owner_signature_identity != exact_identity_)
      throw std::invalid_argument(
          "prepared physical q/C payload names a different chart owner identity");
    if constexpr (std::is_same_v<Scalar, Rational>) {
      try {
        exact_jordan_indicial_ =
            certify_exact_affine_jordan_operator(prepared_);
      } catch (const RecurrenceError& error) {
        exact_jordan_indicial_error_ = error.what();
      }
    }
  }

  json::object solve(const json::object& run, int output_digits) override {
    const auto parse_started = std::chrono::steady_clock::now();
    std::unique_ptr<AcbPrecisionLease> acb_lease;
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      acb_lease = std::make_unique<AcbPrecisionLease>(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
    }
    std::unique_lock<std::recursive_mutex> symbolic_lock;
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      symbolic_lock =
          std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
      SymbolicRational::configure(symbols_);
    }
    RecurrenceProblem<Scalar> problem;
    parse_run_state(run, prepared_, problem);
    const auto parse_ended = std::chrono::steady_clock::now();
    const auto run_parse_ms = std::chrono::duration<double, std::milli>(
        parse_ended - parse_started).count();
    auto result = run_prepared_problem(prepared_, problem, output_digits);
    const auto kernel_ms = result.at("elapsed_ms").as_double();
    const auto run_index = runs_.fetch_add(1) + 1;
    {
      std::lock_guard<std::mutex> lock(stats_mutex_);
      run_parse_ms_ += run_parse_ms;
      kernel_ms_ += kernel_ms;
    }
    result["persistent"] = json::object{
        {"run", run_index}, {"prepare_parse_ms", prepare_parse_ms_},
        {"run_parse_ms", run_parse_ms}, {"static_tensor_copies", 0},
        {"scc_components", scc_.component_count},
        {"scc_coupling_depth", scc_.coupling_depth}};
    return result;
  }

  NativeLocalRun<Scalar> solve_native(
      const json::object& run, const json::object& metadata_object) {
    return solve_native_impl(
        run, metadata_object, std::nullopt, false, std::nullopt, true);
  }

  NativeLocalRun<Scalar> solve_native_with_source(
      const json::object& run, const json::object& metadata_object,
      SourceData<Scalar>&& source) {
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC source injection rejects caller-supplied source data");
    return solve_native_impl(
        run, metadata_object, std::move(source), false, std::nullopt, true);
  }

  void record_native_local_success(
      const NativeLocalDiagnostics& diagnostics) {
    local_runs_.fetch_add(1);
    std::lock_guard<std::mutex> lock(stats_mutex_);
    local_run_parse_ms_ += diagnostics.parse_ms;
    local_kernel_ms_ += diagnostics.kernel_ms;
  }

  std::uint32_t dimension() const override { return prepared_.dimension; }
  std::int32_t frame_base() const override { return prepared_.frame_base; }
  std::uint32_t frame_width() const override { return prepared_.frame_width; }
  std::int32_t chop_digits() const { return prepared_.chop_digits; }
  const char* d0_inverse_mode() const override {
    return prepared_.d0_inverse_scalar.has_value()
        ? "retained-scalar" : "retained-frame";
  }
  std::string regular_value_tail_proxy_max_exact() const override {
    const auto structural_digits = std::min(
        std::max(prepared_.chop_digits / 2, 0), 24);
    Rational result(1);
    for (int digit = 0; digit < structural_digits + 2; ++digit)
      result = result / Rational(10);
    return result.str();
  }
  const char* equation_scalar_domain() const override {
    if constexpr (std::is_same_v<Scalar, Rational>) return "rational";
    if constexpr (std::is_same_v<Scalar, ComplexBall>) return "acb";
    return "symbolic";
  }
  std::shared_ptr<const void> physical_ode_erased() const override {
    return std::static_pointer_cast<const void>(physical_equation_);
  }
  const std::string& physical_payload_identity() const override {
    static const std::string empty;
    return physical_equation_ ? physical_equation_->payload_identity : empty;
  }
  const std::string& physical_payload_record() const override {
    static const std::string empty;
    return physical_equation_ ? physical_equation_->exact_payload_record : empty;
  }
  const std::string& owner_signature_identity() const override {
    static const std::string empty;
    return physical_equation_
        ? physical_equation_->owner_signature_identity : empty;
  }
  slong precision_bits() const { return precision_bits_; }
  bool has_identity_assembly() const {
    return prepared_.assembly_matrix.has_value() &&
           prepared_.assembly_matrix->identity;
  }
  bool uses_epsilon_regular_principal() const {
    return prepared_.epsilon_regular_principal;
  }
  bool has_assembly() const {
    return prepared_.assembly_matrix.has_value();
  }
  const std::string& assembly_exact_identity() const {
    static const std::string empty;
    return prepared_.assembly_matrix.has_value()
        ? prepared_.assembly_matrix->exact_identity : empty;
  }
  FiniteLaurentVector<ComplexBall> apply_assembly_to_acb_matching_vector(
      const FiniteLaurentVector<ComplexBall>& input) const {
    if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
      (void)input;
      throw std::logic_error(
          "only an Acb chart can specialize its matching assembly frame");
    } else {
      if (!prepared_.assembly_matrix.has_value() ||
          input.size() != prepared_.dimension)
        throw std::invalid_argument(
            "matching assembly requires one spectral frame per chart component");
      const auto work_top = matching_detail::checked_power(
          static_cast<std::int64_t>(prepared_.frame_base) +
              prepared_.frame_width - 1,
          "matching assembly work maximum");
      auto framed = detail::zero_block<ComplexBall>(
          prepared_.dimension, prepared_.frame_width);
      for (std::size_t component = 0; component < input.size(); ++component) {
        if (input[component].min_power() < prepared_.frame_base ||
            input[component].complete_max() > work_top)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              "matching assembly input component " +
                  std::to_string(component) + " window [" +
                  std::to_string(input[component].min_power()) + "," +
                  std::to_string(input[component].complete_max()) +
                  "] lies outside prepared [" +
                  std::to_string(prepared_.frame_base) + "," +
                  std::to_string(work_top) + "]");
        for (std::int64_t power = input[component].min_power();
             power <= input[component].complete_max(); ++power)
          framed[component][static_cast<std::size_t>(
              power - prepared_.frame_base)] =
              input[component].coefficient(
                  static_cast<std::int32_t>(power));
      }
      const auto op = recurrence_operator_view(prepared_);
      const auto output = detail::apply_prepared_matrix(
          *prepared_.assembly_matrix, framed, op);
      FiniteLaurentVector<ComplexBall> result;
      result.reserve(prepared_.dimension);
      for (std::uint32_t row = 0; row < prepared_.dimension; ++row) {
        std::optional<std::int32_t> minimum;
        auto complete_max = work_top;
        for (std::uint32_t column = 0; column < prepared_.dimension;
             ++column) {
          const auto valuation = prepared_.assembly_matrix->valuations[
              static_cast<std::size_t>(row) * prepared_.dimension + column];
          if (valuation == kCompleteInfinity) continue;
          const auto candidate_min = matching_detail::checked_power(
              static_cast<std::int64_t>(input[column].min_power()) +
                  valuation,
              "matching assembly output minimum");
          minimum = !minimum.has_value()
              ? candidate_min : std::min(*minimum, candidate_min);
          complete_max = std::min(
              complete_max, matching_detail::checked_power(
                  static_cast<std::int64_t>(
                      input[column].complete_max()) + valuation,
                  "matching assembly output maximum"));
        }
        if (!minimum.has_value() || *minimum < prepared_.frame_base ||
            complete_max < *minimum)
          throw MatchingArithmeticError(
              MatchingArithmeticErrorCode::InsufficientCompleteWindow,
              "matching assembly output does not fit its prepared epsilon frame");
        auto first = *minimum;
        while (first < complete_max &&
               output[row][static_cast<std::size_t>(
                   first - prepared_.frame_base)].is_zero())
          ++first;
        std::vector<ComplexBall> coefficients;
        coefficients.reserve(EpsilonWindow{first, complete_max}.width());
        for (std::int64_t power = first; power <= complete_max; ++power)
          coefficients.push_back(output[row][static_cast<std::size_t>(
              power - prepared_.frame_base)]);
        result.emplace_back(
            EpsilonWindow{first, complete_max}, std::move(coefficients));
      }
      return result;
    }
  }
  bool has_regular_singleton_partition() const {
    if (prepared_.blocks.size() != prepared_.dimension) return false;
    std::vector<std::uint8_t> seen(prepared_.dimension, 0);
    for (const auto& block : prepared_.blocks) {
      if (block.columns.size() != 1 ||
          block.columns.front() >= prepared_.dimension ||
          seen[block.columns.front()])
        return false;
      seen[block.columns.front()] = 1;
    }
    return true;
  }
  std::size_t jordan_block_count() const {
    return prepared_.blocks.size();
  }
  bool jordan_partition_matches(
      const ExactJordanIndicialCertificate& certificate) const {
    if (certificate.dimension != prepared_.dimension ||
        certificate.blocks.size() != prepared_.blocks.size())
      return false;
    for (std::size_t index = 0; index < prepared_.blocks.size(); ++index)
      if (certificate.blocks[index].block_index != index ||
          certificate.blocks[index].columns !=
              prepared_.blocks[index].columns)
        return false;
    return true;
  }
  const std::optional<ExactJordanIndicialCertificate>&
  exact_jordan_indicial() const {
    return exact_jordan_indicial_;
  }
  const std::optional<std::string>& exact_jordan_indicial_error() const {
    return exact_jordan_indicial_error_;
  }

  // Re-expand one retained regular solution about the next regular chart's
  // center.  This is the native equivalent of evaluating the incoming value
  // and running one d-vector recurrence; it avoids constructing and matching
  // a disposable d-column fundamental matrix.
  std::shared_ptr<StoredLocalBase> solve_regular_value_handoff(
      const std::string& local_handle, const json::object& run_prototype,
      json::object metadata_object,
      const std::shared_ptr<StoredLocalBase>& incoming_owner,
      const EpsilonVector& certified_handoff,
      const ExactAffineChart& receiving_geometry,
      const std::vector<Prescription>& receiving_prescriptions,
      EpsilonWindow requested_epsilon,
      std::int32_t required_complete_max,
      const std::string& result_checkpoint_identity,
      json::object derivation,
      std::shared_ptr<void> derivation_owner,
      std::shared_ptr<PhysicalEquationOwnerBase> equation_owner) {
    if (local_handle.empty() || result_checkpoint_identity.empty() ||
        !incoming_owner || !derivation_owner || !equation_owner ||
        equation_owner.get() != this)
      throw std::invalid_argument(
          "regular value handoff lost an identity or strong owner");
    (void)requested_epsilon.width();
    if (required_complete_max < requested_epsilon.min_power ||
        required_complete_max > requested_epsilon.complete_max)
      throw std::invalid_argument(
          "regular value handoff epsilon contract is inconsistent");
    if (!physical_equation_ || !has_identity_assembly() ||
        !has_regular_singleton_partition())
      throw std::invalid_argument(
          "regular value handoff receiver is not a physical identity-frame value solver");

    auto incoming =
        std::dynamic_pointer_cast<StoredLocal<Scalar>>(incoming_owner);
    if (!incoming || incoming->solution().dimension != prepared_.dimension)
      throw std::invalid_argument(
          "regular value handoff coefficient domain or dimension changed");
    const auto& source = incoming->solution();
    validate_local_solution(source, false);
    if (!source.error.empty() || source.sectors.size() != 1 ||
        source.sectors.front().a.domain != ExactDomain::Rational ||
        source.sectors.front().b.domain != ExactDomain::Rational ||
        !(Rational(source.sectors.front().a.canonical) == Rational(0)) ||
        !(Rational(source.sectors.front().b.canonical) == Rational(0)) ||
        source.sectors.front().log_power != 0)
      throw std::invalid_argument(
          "regular value handoff source is not one certified (0,0,0) sector");
    if (source.epsilon.min_power < prepared_.frame_base)
      throw std::domain_error(
          "regular value handoff receiver frame would discard certified lower epsilon coefficients");
    const auto source_complete_max = std::min(
        source.epsilon.complete_max, incoming->top_valid());
    if (source_complete_max < requested_epsilon.complete_max)
      throw std::domain_error(
          "regular value handoff source does not cover its requested complete epsilon window");

    AcbPrecisionLease acb_lease(precision_bits_);
    ComplexBall::set_precision(precision_bits_);
    auto run = run_prototype;
    const auto parse_started = std::chrono::steady_clock::now();
    RecurrenceProblem<Scalar> problem;
    parse_run_state(run, prepared_, problem);
    if (problem.log_max != 0 || !problem.has_initial ||
        problem.adaptive_lower_frame_probe || problem.source.has_value() ||
        problem.return_u || !ScalarTraits<Scalar>::is_zero(problem.a_target) ||
        !ScalarTraits<Scalar>::is_zero(problem.b_target) ||
        problem.a_shift_min != 0 ||
        problem.a_shifts.size() !=
            static_cast<std::size_t>(problem.nmax) + 1 ||
        problem.schedule.size() !=
            static_cast<std::size_t>(problem.nmax) + 1)
      throw std::invalid_argument(
          "regular value handoff prototype is not one homogeneous (0,0,0) value run");
    const auto scalar_identical = [](const Scalar& left,
                                     const Scalar& right) {
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        return acb_equal(left.raw(), right.raw());
      else
        return left == right;
    };
    for (std::size_t n = 0; n < problem.schedule.size(); ++n) {
      const auto expected_shift = ScalarTraits<Scalar>::integer(
          static_cast<long>(n));
      if (!scalar_identical(problem.a_shifts[n], expected_shift))
        throw std::invalid_argument(
            "regular value handoff prototype a-shifts are not the exact Taylor indices");
      if (problem.schedule[n].size() != prepared_.blocks.size())
        throw std::invalid_argument(
            "regular value handoff prototype schedule has the wrong block partition");
      for (const auto& step : problem.schedule[n]) {
        const auto expected = n == 0 ? StepCase::Resonant : StepCase::Taylor;
        if (step.kind != expected ||
            !scalar_identical(step.d_a, expected_shift) ||
            !ScalarTraits<Scalar>::is_zero(step.d_b))
          throw std::invalid_argument(
              "regular value handoff prototype schedule is not resonant at zero and Taylor by exact index");
      }
    }

    const auto frame_width = static_cast<std::size_t>(prepared_.frame_width);
    problem.initial.assign(
        static_cast<std::size_t>(prepared_.dimension) * frame_width,
        ScalarTraits<Scalar>::zero());
    problem.initial_validity.assign(prepared_.dimension,
                                    source_complete_max);
    if constexpr (std::is_same_v<Scalar, Rational>) {
      throw std::invalid_argument(
          "Rational regular value handoff requires an exact polynomial-tail-zero proof");
    } else {
      if (certified_handoff.dimension != prepared_.dimension ||
          certified_handoff.epsilon.min_power !=
              source.epsilon.min_power ||
          certified_handoff.epsilon.complete_max !=
              source.epsilon.complete_max ||
          !certified_handoff.error.empty())
        throw std::invalid_argument(
            "regular value handoff certified vector changed its source frame or retained an unabsorbed error envelope");
      for (std::uint32_t component = 0; component < prepared_.dimension;
           ++component)
        for (std::int64_t raw_power = source.epsilon.min_power;
             raw_power <= source_complete_max; ++raw_power) {
          const auto power = static_cast<std::int32_t>(raw_power);
          const auto frame = static_cast<std::size_t>(
              static_cast<std::int64_t>(power) - prepared_.frame_base);
          if (frame >= frame_width)
            throw std::domain_error(
                "regular value handoff source exceeds the receiver value frame");
          if (power >= certified_handoff.epsilon.min_power)
            problem.initial[
                static_cast<std::size_t>(component) * frame_width + frame] =
                certified_handoff.at(power, component);
        }
    }

    metadata_object["checkpoint_identity"] = result_checkpoint_identity;
    auto metadata = parse_local_metadata(metadata_object);
    const auto same_prescription = [](const Prescription& left,
                                      const Prescription& right) {
      return left.factor_exact == right.factor_exact &&
             left.sign == right.sign &&
             left.multiplicity == right.multiplicity &&
             left.leading_coefficient_sign ==
                 right.leading_coefficient_sign;
    };
    if (!(Rational(metadata.chart.center_exact) ==
              receiving_geometry.center) ||
        !(Rational(metadata.chart.scale_exact) == receiving_geometry.scale) ||
        metadata.chart.infinite_radius ||
        !acb_equal(metadata.chart.radius.raw(),
                   ComplexBall::from_strings(
                       receiving_geometry.radius.str()).raw()) ||
        metadata.prescriptions.size() != receiving_prescriptions.size() ||
        !std::equal(metadata.prescriptions.begin(),
                    metadata.prescriptions.end(),
                    receiving_prescriptions.begin(), same_prescription))
      throw std::invalid_argument(
          "regular value handoff prototype metadata differs from its retained plan chart");
    verify_tag_binding(metadata.a, problem.a_target,
                       "regular value handoff a");
    verify_tag_binding(metadata.b, problem.b_target,
                       "regular value handoff b");
    const auto parse_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - parse_started).count();

    const auto kernel_started = std::chrono::steady_clock::now();
    auto recurrence = RecurrenceSolver<Scalar>(problem, prepared_).run();
    auto assembled = assemble_recurrence(prepared_, problem, recurrence);
    auto solution = make_local_solution(
        problem, std::move(assembled), std::move(metadata));
    validate_local_solution(solution, false);
    if (solution.epsilon.complete_max < requested_epsilon.complete_max ||
        recurrence.top_valid < requested_epsilon.complete_max)
      throw std::domain_error(
          "regular value handoff solve did not preserve the requested complete epsilon window");
    auto tail_model = prepare_regular_homogeneous_tail_model(
        prepared_, problem, solution, exact_identity_);
    auto pseudo_hits = std::move(recurrence.hits);
    if (!pseudo_hits.empty())
      throw std::domain_error(
          "regular value handoff unexpectedly encountered a pseudo resonance");
    const auto kernel_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - kernel_started).count();
    const NativeLocalDiagnostics diagnostics{
        recurrence.top_valid, parse_ms, kernel_ms};

    derivation["evaluated_epsilon"] = json::object{
        {"min", source.epsilon.min_power},
        {"max", source_complete_max},
        {"required_complete_max", required_complete_max}};
    derivation["output"] = json::object{
        {"checkpoint_identity", solution.checkpoint_identity},
        {"chart", handle_},
        {"source_operator_identity", exact_identity_},
        {"epsilon", json::object{{"min", solution.epsilon.min_power},
                                  {"max", solution.epsilon.complete_max}}},
        {"taylor_complete_max", solution.taylor_complete_max},
        {"top_valid", encode_validity(recurrence.top_valid)},
        {"dimension", solution.dimension}};
    derivation["equation_owner_signature_identity"] =
        equation_owner->owner_signature_identity();
    derivation["equation_payload_identity"] =
        equation_owner->physical_payload_identity();
    derivation["provenance_identity"] = json::serialize(
        canonical_json_value(derivation));

    auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, handle_, exact_identity_, std::move(solution),
        precision_bits_, std::move(pseudo_hits), diagnostics, std::nullopt,
        std::move(derivation), std::move(derivation_owner),
        std::move(tail_model), std::nullopt, true, true,
        std::move(equation_owner), physical_equation_);
    record_native_local_success(diagnostics);
    return local;
  }
  ChartStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {runs_.load(), local_runs_.load(), prepare_parse_ms_,
            run_parse_ms_, kernel_ms_, local_run_parse_ms_,
            local_kernel_ms_};
  }

 private:
  NativeLocalRun<Scalar> solve_native_impl(
      const json::object& run, const json::object& metadata_object,
      std::optional<SourceData<Scalar>> native_source,
      bool attach_tail_model,
      std::optional<std::string> tail_operator_identity,
      bool allow_scc_deferred_completeness) {
    if (precision_bits_ < 64)
      throw std::invalid_argument(
          "native local solutions require at least 64 bits of Acb precision");
    // LocalSolution always carries numeric chart geometry and exact-tag
    // specializations even when its coefficient field is exact.  Lease the
    // output precision before parsing any such ball.
    AcbPrecisionLease acb_lease(precision_bits_);
    ComplexBall::set_precision(precision_bits_);
    std::unique_lock<std::recursive_mutex> symbolic_lock;
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      symbolic_lock =
          std::unique_lock<std::recursive_mutex>(symbolic_run_mutex());
      SymbolicRational::configure(symbols_);
    }

    const auto parse_started = std::chrono::steady_clock::now();
    RecurrenceProblem<Scalar> problem;
    parse_run_state(run, prepared_, problem);
    if (native_source.has_value()) {
      if (problem.source.has_value())
        throw std::invalid_argument(
            "native source injection found an unexpected parsed source");
      problem.source = std::move(*native_source);
    }
    if (!prepared_.assembly_matrix.has_value())
      throw std::invalid_argument(
          "local.solve requires a retained chart with native assembly");
    auto metadata = parse_local_metadata(metadata_object);
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      for (const auto* descriptor : {&metadata.a, &metadata.b})
        if (descriptor->domain == ExactDomain::SymbolicRational &&
            descriptor->symbols != symbols_)
          throw std::invalid_argument(
              "local exact-tag regulator field differs from its chart session");
    }
    verify_tag_binding(metadata.a, problem.a_target, "local a");
    verify_tag_binding(metadata.b, problem.b_target, "local b");
    const auto parse_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - parse_started).count();

    const auto kernel_started = std::chrono::steady_clock::now();
    auto recurrence = RecurrenceSolver<Scalar>(problem, prepared_).run();
    auto assembled = allow_scc_deferred_completeness
        ? assemble_scc_recurrence_candidate(prepared_, problem, recurrence)
        : SCCAssemblyCandidate<Scalar>{
              assemble_recurrence(prepared_, problem, recurrence), false,
              recurrence.top_valid};
    require_exact_domain_for_deferred_scc_candidate(assembled);
    auto solution = make_local_solution(
        problem, std::move(assembled.coefficients), std::move(metadata));
    validate_local_solution(solution, false);
    auto tail_model = unavailable_tail_model(
        attach_tail_model
            ? "symbolic coefficient locals do not support numeric certified tail bounds"
            : "tail model was not requested for this internal native block run");
    if (attach_tail_model) {
      if constexpr (std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>)
        tail_model = prepare_regular_homogeneous_tail_model(
            prepared_, problem, solution,
            tail_operator_identity.value_or(exact_identity_));
    }
    auto pseudo_hits = std::move(recurrence.hits);
    const auto kernel_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - kernel_started).count();
    NativeLocalDiagnostics diagnostics{
        recurrence.top_valid, parse_ms, kernel_ms};
    diagnostics.requires_parent_completeness_certificate =
        assembled.requires_parent_certificate;
    return {std::move(solution), std::move(pseudo_hits),
            diagnostics,
            std::move(tail_model)};
  }

  std::shared_ptr<StoredLocalBase> solve_local(
      const std::string& local_handle, const json::object& run,
      const json::object& metadata_object,
      std::shared_ptr<PhysicalEquationOwnerBase> equation_owner) override {
    if (!equation_owner)
      throw std::invalid_argument(
          "primitive local solve lost its physical equation owner");
    std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
        retained_physical_equation = physical_equation_;
    std::string source_handle = handle_;
    std::string source_identity = exact_identity_;
    if (equation_owner.get() != this) {
      auto regular = std::dynamic_pointer_cast<
          RegularPhysicalEquationOwner<Scalar>>(equation_owner);
      if (!regular || !physical_equation_ ||
          regular->dimension() != prepared_.dimension ||
          !geometry_record_.has_value() ||
          *geometry_record_ != regular->geometry_record() ||
          physical_equation_->payload_identity !=
              regular->physical_equation()->payload_identity)
        throw std::invalid_argument(
            "framed fallback chart differs from its regular physical equation owner");
      const auto payload_without_owner = [](const std::string& record) {
        auto object = as_object(json::parse(record),
                                "physical q/C payload identity");
        object.erase("owner_signature_identity");
        return json::serialize(canonical_json_value(object));
      };
      if (payload_without_owner(physical_equation_->exact_payload_record) !=
          payload_without_owner(
              regular->physical_equation()->exact_payload_record))
        throw std::invalid_argument(
            "framed fallback chart q/C payload differs from its regular physical equation owner");
      const auto metadata = parse_local_metadata(metadata_object);
      const auto retained_geometry = parse_retained_composite_geometry(
          json::parse(regular->geometry_record()));
      const auto same_prescription = [](const Prescription& left,
                                        const Prescription& right) {
        return left.factor_exact == right.factor_exact &&
               left.sign == right.sign &&
               left.multiplicity == right.multiplicity &&
               left.leading_coefficient_sign ==
                   right.leading_coefficient_sign;
      };
      if (metadata.chart.center_exact !=
              retained_geometry.chart.center_exact ||
          metadata.chart.scale_exact !=
              retained_geometry.chart.scale_exact ||
          metadata.chart.infinite_radius !=
              retained_geometry.chart.infinite_radius ||
          !acb_equal(metadata.chart.radius.raw(),
                     retained_geometry.chart.radius.raw()) ||
          metadata.prescriptions.size() !=
              retained_geometry.prescriptions.size() ||
          !std::equal(metadata.prescriptions.begin(),
                      metadata.prescriptions.end(),
                      retained_geometry.prescriptions.begin(),
                      same_prescription))
        throw std::invalid_argument(
            "framed fallback local metadata differs from its regular equation owner geometry");
      retained_physical_equation = regular->physical_equation();
      source_handle = regular->handle();
      source_identity = regular->exact_identity();
    }
    auto native = solve_native_impl(
        run, metadata_object, std::nullopt, true, source_identity, false);
    std::shared_ptr<PhysicalEquationOwnerBase> retained_equation_owner;
    std::string residual_unavailable_reason;
    const bool homogeneous = run.at("source").is_null();
    if constexpr (std::is_same_v<Scalar, Rational> ||
                  std::is_same_v<Scalar, ComplexBall>) {
      if (!homogeneous) {
        retained_physical_equation.reset();
        residual_unavailable_reason =
            "owner-bound residual is unsupported for sourced primitive locals until the exact physical source is retained";
      } else if (!retained_physical_equation) {
        residual_unavailable_reason =
            "owner-bound residual is unsupported: the prepared chart predates the physical q/C payload";
      } else {
        retained_equation_owner = std::move(equation_owner);
      }
    }
    auto local = make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, source_handle, source_identity,
        std::move(native.solution),
        precision_bits_,
        std::move(native.pseudo_hits), native.diagnostics, std::nullopt,
        std::nullopt, nullptr,
        std::move(native.tail_model), std::nullopt, true, true,
        std::move(retained_equation_owner),
        std::move(retained_physical_equation),
        std::move(residual_unavailable_reason));
    record_native_local_success(native.diagnostics);
    return local;
  }
  PreparedRecurrenceOperator<Scalar> prepared_;
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>> physical_equation_;
  std::optional<ExactJordanIndicialCertificate> exact_jordan_indicial_;
  std::optional<std::string> exact_jordan_indicial_error_;
  slong precision_bits_ = 256;
  std::vector<std::string> symbols_;
  std::atomic<std::uint64_t> runs_{0};
  std::atomic<std::uint64_t> local_runs_{0};
  mutable std::mutex stats_mutex_;
  double run_parse_ms_ = 0.0;
  double kernel_ms_ = 0.0;
  double local_run_parse_ms_ = 0.0;
  double local_kernel_ms_ = 0.0;
};

struct CompositeColumnSolveResult {
  std::shared_ptr<StoredLocalBase> local;
  json::array block_diagnostics;
  double elapsed_ms = 0.0;
};

struct CompositeWorkContract;

class CompositeSCCChartBase : public PhysicalEquationOwnerBase {
 public:
  CompositeSCCChartBase(std::string handle, std::string key,
                        std::string exact_identity, std::string signature,
                        std::string rational_shadow_identity)
      : handle_(std::move(handle)), key_(std::move(key)),
        exact_identity_(std::move(exact_identity)),
        signature_(std::move(signature)),
        rational_shadow_identity_(std::move(rational_shadow_identity)) {
    if (rational_shadow_identity_.empty())
      throw std::invalid_argument(
          "composite SCC rational-shadow identity cannot be empty");
  }
  virtual ~CompositeSCCChartBase() = default;

  virtual json::object stats_json() const = 0;
  virtual CompositeColumnSolveResult solve_column(
      const std::string& local_handle, const json::object& request,
      std::shared_ptr<CompositeSCCChartBase> equation_owner) = 0;
  virtual const char* column_execution_capability() const = 0;
  virtual const std::string& geometry_record() const = 0;
  virtual std::uint32_t dimension() const = 0;
  virtual const CompositeWorkContract& work_contract() const = 0;
  // For an epsilon-regular singular composite, express a physical parent
  // value in the exact block-spectral normal frame used at the residue.  The
  // same invertible left transformation is applied to the receiving basis
  // and incoming value, so matching weights are unchanged while the Laurent
  // solve no longer sees the confluent physical columns.
  virtual std::optional<std::pair<
      FiniteLaurentVector<ComplexBall>, std::string>>
  normalize_acb_matching_vector(
      const FiniteLaurentVector<ComplexBall>& physical) const override = 0;
  virtual std::vector<std::shared_ptr<PreparedChartBase>>
  dependency_charts() const = 0;
  const std::string& equation_owner_handle() const override {
    return handle_;
  }
  const std::string& equation_operator_identity() const override {
    return exact_identity_;
  }
  const char* equation_owner_kind() const override {
    return "composite-scc";
  }
  const std::string& handle() const { return handle_; }
  const std::string& key() const { return key_; }
  const std::string& exact_identity() const { return exact_identity_; }
  const std::string& signature() const { return signature_; }
  const std::string& rational_shadow_identity() const {
    return rational_shadow_identity_;
  }
  std::optional<std::uint32_t> matching_scc_dimension() const override {
    return dimension();
  }
  const char* matching_scc_column_execution_capability() const override {
    return column_execution_capability();
  }
  const std::string* matching_scc_rational_shadow_identity() const override {
    return &rational_shadow_identity_;
  }

 protected:
  std::string handle_;
  std::string key_;
  std::string exact_identity_;
  std::string signature_;
  std::string rational_shadow_identity_;
};

struct CompositeWorkContract {
  std::int32_t work_min = 0;
  std::optional<std::int32_t> cancellation_audit_base;
  std::int32_t requested_min = 0;
  std::int32_t requested_max = 0;
  std::int32_t work_complete_max = 0;
  std::uint32_t public_t_order = 0;
  std::uint32_t work_t_order = 0;
  std::uint32_t wolfram_coupling_depth = 0;
};

template <typename Scalar>
struct CompositeSpectralSourceTransform {
  bool identity = true;
  bool epsilon_unimodular = true;
  std::int32_t det_epsilon_valuation = 0;
  std::string producer_identity;
  std::string v_exact_identity;
  std::string vinv_exact_identity;
  std::string det_exact_identity;
  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
};

template <typename Scalar>
struct CompositeGaugeTransform {
  bool identity = true;
  std::string role;
  std::string producer_identity;
  std::string gauge_exact_identity;
  std::string gauge_inverse_exact_identity;
  std::string gauge_det_exact_identity;
  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
};

template <typename Scalar>
struct CompositeSCCBlock {
  std::uint32_t block = 0;
  std::vector<std::uint32_t> vertices;
  std::string source_handle;
  std::string principal_identity;
  bool regular = true;
  bool no_pseudo = false;
  CompositeSpectralSourceTransform<Scalar> source_transform;
  CompositeGaugeTransform<Scalar> to_physical;
  CompositeGaugeTransform<Scalar> to_reduced;
  std::optional<ExactJordanIndicialCertificate> exact_jordan_indicial;
  std::shared_ptr<PreparedChart<Scalar>> chart;
};

template <typename Scalar>
CompositeGaugeTransform<Scalar> parse_composite_gauge_transform(
    const json::value& raw, const std::string& expected_role,
    std::uint32_t dimension, std::uint32_t frame_width,
    const CompositeWorkContract& work, const std::string& domain,
    const std::vector<std::string>& symbols) {
  const auto& object = as_object(raw, "SCC gauge transform");
  require_exact_keys(object,
      {"schema", "role", "rows", "columns", "identity",
       "gauge_exact_identity", "gauge_inverse_exact_identity",
       "gauge_det_exact_identity", "exact_identity", "domain",
       "symbols", "entries"}, "SCC gauge transform");
  if (required_string(object, "schema") !=
          "diffexp2-scc-gauge-transform-v1" ||
      required_string(object, "role") != expected_role)
    throw std::invalid_argument("unsupported SCC gauge transform role/schema");
  const auto rows = as_u32(object.at("rows"), "gauge transform rows");
  const auto columns = as_u32(object.at("columns"), "gauge transform columns");
  if (dimension == 0 || rows != dimension || columns != dimension ||
      required_string(object, "domain") != domain ||
      parse_symbols(object) != symbols)
    throw std::invalid_argument("SCC gauge transform shape/field differs from its block");
  CompositeGaugeTransform<Scalar> transform;
  transform.identity = object.at("identity").as_bool();
  transform.role = expected_role;
  transform.producer_identity = required_string(object, "exact_identity");
  transform.gauge_exact_identity = required_string(object, "gauge_exact_identity");
  transform.gauge_inverse_exact_identity = required_string(
      object, "gauge_inverse_exact_identity");
  transform.gauge_det_exact_identity = required_string(
      object, "gauge_det_exact_identity");
  if (transform.producer_identity.empty() ||
      transform.gauge_exact_identity.empty() ||
      transform.gauge_inverse_exact_identity.empty() ||
      transform.gauge_det_exact_identity.empty())
    throw std::invalid_argument("SCC gauge transform exact identities must be nonempty");
  json::object proof;
  try {
    proof = as_object(json::parse(transform.producer_identity),
                      "SCC gauge transform identity");
  } catch (const std::exception&) {
    throw std::invalid_argument("SCC gauge transform identity is not structural JSON");
  }
  require_exact_keys(proof,
      {"schema", "role", "dimension", "identity",
       "gauge_exact_identity", "gauge_inverse_exact_identity",
       "gauge_det_exact_identity", "source_window", "entries"},
      "SCC gauge transform identity");
  if (required_string(proof, "schema") !=
          "diffexp2-scc-gauge-transform-identity-v1" ||
      required_string(proof, "role") != expected_role ||
      as_u32(proof.at("dimension"), "gauge identity dimension") != dimension ||
      proof.at("identity").as_bool() != transform.identity ||
      required_string(proof, "gauge_exact_identity") != transform.gauge_exact_identity ||
      required_string(proof, "gauge_inverse_exact_identity") != transform.gauge_inverse_exact_identity ||
      required_string(proof, "gauge_det_exact_identity") != transform.gauge_det_exact_identity)
    throw std::invalid_argument("SCC gauge transform differs from its exact identity");
  const auto& window = as_object(proof.at("source_window"),
                                 "SCC gauge transform window");
  require_exact_keys(window,
      {"epsilon_min", "epsilon_complete_max", "taylor_complete_max"},
      "SCC gauge transform window");
  if (as_i32(window.at("epsilon_min"), "gauge epsilon minimum") != work.work_min ||
      as_i32(window.at("epsilon_complete_max"), "gauge epsilon maximum") != work.work_complete_max ||
      as_u32(window.at("taylor_complete_max"), "gauge Taylor maximum") != work.work_t_order)
    throw std::invalid_argument("SCC gauge transform names a different work rectangle");
  transform.matrix.rows = rows;
  transform.matrix.columns = columns;
  transform.matrix.exact_identity = transform.producer_identity;
  std::optional<std::pair<std::uint32_t, std::uint32_t>> previous;
  std::vector<std::uint8_t> active_rows(rows, 0), active_columns(columns, 0);
  for (const auto& value : as_array(object.at("entries"), "SCC gauge entries")) {
    const auto& entry = as_object(value, "SCC gauge entry");
    require_exact_keys(entry, {"row", "column", "exact_entry", "multiplier"},
                       "SCC gauge entry");
    const auto row = as_u32(entry.at("row"), "gauge row");
    const auto column = as_u32(entry.at("column"), "gauge column");
    const auto location = std::make_pair(row, column);
    if (row >= rows || column >= columns ||
        (previous.has_value() && *previous >= location))
      throw std::invalid_argument("SCC gauge entries must be unique and row-major");
    previous = location;
    const auto exact_entry = required_string(entry, "exact_entry");
    const auto& raw_multiplier = as_object(entry.at("multiplier"),
                                           "SCC gauge multiplier");
    auto multiplier = parse_prepared_rational_taylor_multiplier<Scalar>(
        raw_multiplier, frame_width,
        static_cast<std::size_t>(work.work_t_order) + 1, true,
        "SCC gauge multiplier");
    if (multiplier.proven_zero || multiplier.exact_identity != exact_entry)
      throw std::invalid_argument("SCC gauge entry is zero or changed exact identity");
    active_rows[row] = 1; active_columns[column] = 1;
    transform.matrix.entries.push_back(
        typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
            row, column, std::move(multiplier)});
  }
  if (std::any_of(active_rows.begin(), active_rows.end(), [](auto x){return x==0;}) ||
      std::any_of(active_columns.begin(), active_columns.end(), [](auto x){return x==0;}))
    throw std::invalid_argument("SCC gauge transform has an empty row or column");
  const auto& identity_entries = as_array(proof.at("entries"),
                                          "SCC gauge identity entries");
  if (identity_entries.size() != transform.matrix.entries.size())
    throw std::invalid_argument("SCC gauge entries differ from exact identity");
  for (std::size_t i = 0; i < identity_entries.size(); ++i) {
    const auto& item = as_object(identity_entries[i], "SCC gauge identity entry");
    require_exact_keys(item, {"row", "column", "exact_entry",
      "epsilon_shift", "center_pole_order"}, "SCC gauge identity entry");
    const auto& prepared = transform.matrix.entries[i];
    if (as_u32(item.at("row"), "gauge identity row") != prepared.row ||
        as_u32(item.at("column"), "gauge identity column") != prepared.column ||
        required_string(item, "exact_entry") != prepared.multiplier.exact_identity ||
        as_i32(item.at("epsilon_shift"), "gauge identity shift") != prepared.multiplier.epsilon_shift ||
        as_u32(item.at("center_pole_order"), "gauge identity pole") != prepared.multiplier.center_pole_order)
      throw std::invalid_argument("prepared SCC gauge entry differs from exact identity");
  }
  return transform;
}

template <typename Scalar>
CompositeSpectralSourceTransform<Scalar>
parse_composite_spectral_source_transform(
    const json::value& raw, std::uint32_t dimension,
    std::uint32_t frame_width, const CompositeWorkContract& work,
    const std::string& domain, const std::vector<std::string>& symbols,
    bool identity_v, const std::string& assembly_exact_identity) {
  const auto& object = as_object(
      raw, "SCC spectral source transform");
  require_exact_keys(object,
      {"schema", "rows", "columns", "identity",
       "epsilon_unimodular", "det_epsilon_valuation",
       "v_exact_identity", "vinv_exact_identity", "det_exact_identity",
       "exact_identity", "domain", "symbols", "entries"},
      "SCC spectral source transform");
  if (required_string(object, "schema") !=
      "diffexp2-scc-spectral-source-transform-v1")
    throw std::invalid_argument(
        "unsupported SCC spectral source-transform schema");
  const auto rows = as_u32(object.at("rows"),
                           "spectral source-transform rows");
  const auto columns = as_u32(object.at("columns"),
                              "spectral source-transform columns");
  if (rows != dimension || columns != dimension || dimension == 0)
    throw std::invalid_argument(
        "SCC spectral source-transform dimensions differ from its block");
  if (required_string(object, "domain") != domain ||
      parse_symbols(object) != symbols)
    throw std::invalid_argument(
        "SCC spectral source-transform field differs from its session");

  CompositeSpectralSourceTransform<Scalar> transform;
  transform.identity = object.at("identity").as_bool();
  transform.epsilon_unimodular =
      object.at("epsilon_unimodular").as_bool();
  transform.det_epsilon_valuation = as_i32(
      object.at("det_epsilon_valuation"),
      "spectral determinant epsilon valuation");
  if (transform.identity != identity_v)
    throw std::invalid_argument(
        "SCC spectral source-transform identity claim differs from identity_v");
  if (!transform.epsilon_unimodular)
    throw std::invalid_argument(
        "SCC spectral source transform requires an exact Laurent-unimodular frame");
  transform.producer_identity = required_string(object, "exact_identity");
  transform.v_exact_identity = required_string(
      object, "v_exact_identity");
  transform.vinv_exact_identity = required_string(
      object, "vinv_exact_identity");
  transform.det_exact_identity = required_string(
      object, "det_exact_identity");
  if (transform.producer_identity.empty() ||
      transform.v_exact_identity.empty() ||
      transform.vinv_exact_identity.empty() ||
      transform.det_exact_identity.empty())
    throw std::invalid_argument(
        "SCC spectral source-transform exact identities must be nonempty");
  if (assembly_exact_identity.empty() ||
      transform.v_exact_identity != assembly_exact_identity)
    throw std::invalid_argument(
        "SCC spectral V identity differs from the retained assembly operator");

  json::object producer_record;
  try {
    producer_record = as_object(
        json::parse(transform.producer_identity),
        "SCC spectral source-transform exact identity");
  } catch (const std::exception&) {
    throw std::invalid_argument(
        "SCC spectral source-transform exact identity is not a valid structural JSON certificate");
  }
  require_exact_keys(producer_record,
      {"schema", "state_basis", "target_recurrence_basis", "dimension",
       "identity", "epsilon_unimodular", "det_epsilon_valuation",
       "v_exact_identity", "vinv_exact_identity", "det_exact_identity",
       "source_window", "serialization", "entries"},
      "SCC spectral source-transform exact identity");
  if (required_string(producer_record, "schema") !=
          "diffexp2-scc-spectral-source-transform-identity-v1" ||
      required_string(producer_record, "state_basis") !=
          "reduced-g-after-spectral-assembly" ||
      required_string(producer_record, "target_recurrence_basis") !=
          "spectral-u" ||
      as_u32(producer_record.at("dimension"),
             "spectral identity dimension") != dimension ||
      producer_record.at("identity").as_bool() != transform.identity ||
      !producer_record.at("epsilon_unimodular").as_bool() ||
      as_i32(producer_record.at("det_epsilon_valuation"),
             "spectral identity determinant valuation") !=
          transform.det_epsilon_valuation ||
      required_string(producer_record, "v_exact_identity") !=
          transform.v_exact_identity ||
      required_string(producer_record, "vinv_exact_identity") !=
          transform.vinv_exact_identity ||
      required_string(producer_record, "det_exact_identity") !=
          transform.det_exact_identity)
    throw std::invalid_argument(
        "SCC spectral source-transform fields differ from their exact structural identity");
  const auto& source_window = as_object(
      producer_record.at("source_window"),
      "spectral source-transform identity window");
  require_exact_keys(source_window,
      {"epsilon_min", "epsilon_complete_max", "taylor_complete_max"},
      "spectral source-transform identity window");
  if (as_i32(source_window.at("epsilon_min"),
             "spectral source epsilon minimum") != work.work_min ||
      as_i32(source_window.at("epsilon_complete_max"),
             "spectral source epsilon maximum") !=
          work.work_complete_max ||
      as_u32(source_window.at("taylor_complete_max"),
             "spectral source Taylor maximum") != work.work_t_order)
    throw std::invalid_argument(
        "SCC spectral source-transform exact identity names a different work rectangle");
  const auto& serialization = as_object(
      producer_record.at("serialization"),
      "spectral source-transform identity field");
  require_exact_keys(serialization, {"domain", "symbols"},
                     "spectral source-transform identity field");
  if (required_string(serialization, "domain") != domain)
    throw std::invalid_argument(
        "SCC spectral source-transform exact identity names a different scalar domain");
  const auto& identity_symbols = as_array(
      serialization.at("symbols"),
      "spectral source-transform identity symbols");
  if (identity_symbols.size() != symbols.size())
    throw std::invalid_argument(
        "SCC spectral source-transform exact identity has a different regulator field");
  for (std::size_t index = 0; index < symbols.size(); ++index) {
    const auto& symbol = as_object(
        identity_symbols[index],
        "spectral source-transform identity symbol");
    require_exact_keys(symbol, {"context", "name"},
                       "spectral source-transform identity symbol");
    if (required_string(symbol, "context") != "Global`" ||
        required_string(symbol, "name") != symbols[index])
      throw std::invalid_argument(
          "SCC spectral source-transform exact identity regulator differs from its session");
  }

  transform.matrix.rows = rows;
  transform.matrix.columns = columns;
  transform.matrix.exact_identity = transform.producer_identity;
  std::optional<std::pair<std::uint32_t, std::uint32_t>> previous;
  std::vector<std::uint8_t> active_rows(rows, 0);
  std::vector<std::uint8_t> active_columns(columns, 0);
  for (const auto& raw_entry_value : as_array(
           object.at("entries"), "spectral source-transform entries")) {
    const auto& raw_entry = as_object(
        raw_entry_value, "spectral source-transform entry");
    require_exact_keys(raw_entry,
        {"row", "column", "exact_entry", "multiplier"},
        "spectral source-transform entry");
    const auto row = as_u32(raw_entry.at("row"),
                            "spectral source-transform row");
    const auto column = as_u32(raw_entry.at("column"),
                               "spectral source-transform column");
    const auto location = std::make_pair(row, column);
    if (row >= rows || column >= columns ||
        (previous.has_value() && *previous >= location))
      throw std::invalid_argument(
          "SCC spectral source-transform entries must be unique, in range, and row-major");
    previous = location;
    const auto exact_entry = required_string(raw_entry, "exact_entry");
    if (exact_entry.empty())
      throw std::invalid_argument(
          "SCC spectral VInv entry identity must be nonempty");
    const auto& raw_multiplier = as_object(
        raw_entry.at("multiplier"),
        "prepared spectral source-transform multiplier");
    auto multiplier = parse_prepared_rational_taylor_multiplier<Scalar>(
        raw_multiplier, frame_width,
        static_cast<std::size_t>(work.work_t_order) + 1, true,
        "prepared spectral source-transform multiplier");
    const auto shifted_work_min = static_cast<std::int64_t>(work.work_min) +
        multiplier.epsilon_shift;
    const auto shifted_work_max =
        static_cast<std::int64_t>(work.work_complete_max) +
        multiplier.epsilon_shift;
    if (shifted_work_min < std::numeric_limits<std::int32_t>::min() ||
        shifted_work_min > std::numeric_limits<std::int32_t>::max() ||
        shifted_work_max < std::numeric_limits<std::int32_t>::min() ||
        shifted_work_max > std::numeric_limits<std::int32_t>::max())
      throw std::invalid_argument(
          "SCC spectral source-transform shift overflows the retained epsilon frame");
    if (multiplier.center_pole_order != 0)
      throw std::invalid_argument(
          "t-independent SCC spectral VInv cannot have a center pole");
    if (multiplier.exact_identity != exact_entry)
      throw std::invalid_argument(
          "prepared spectral multiplier identity differs from its exact VInv entry");
    if (multiplier.proven_zero)
      throw std::invalid_argument(
          "structurally zero SCC spectral source-transform entries must be omitted");
    active_rows[row] = 1;
    active_columns[column] = 1;
    transform.matrix.entries.push_back(
        typename PreparedSparseLocalMultiplierMatrix<Scalar>::Entry{
            row, column, std::move(multiplier)});
  }
  if (std::any_of(active_rows.begin(), active_rows.end(),
                  [](const auto active) { return active == 0; }) ||
      std::any_of(active_columns.begin(), active_columns.end(),
                  [](const auto active) { return active == 0; }))
    throw std::invalid_argument(
        "certified SCC spectral VInv has an empty structural row or column");
  const auto& identity_entries = as_array(
      producer_record.at("entries"),
      "spectral source-transform identity entries");
  if (identity_entries.size() != transform.matrix.entries.size())
    throw std::invalid_argument(
        "SCC spectral source-transform entries differ from their exact structural identity");
  for (std::size_t index = 0; index < identity_entries.size(); ++index) {
    const auto& identity_entry = as_object(
        identity_entries[index],
        "spectral source-transform identity entry");
    require_exact_keys(identity_entry,
        {"row", "column", "exact_entry", "epsilon_shift",
         "center_pole_order"},
        "spectral source-transform identity entry");
    const auto& prepared_entry = transform.matrix.entries[index];
    if (as_u32(identity_entry.at("row"),
               "spectral identity entry row") != prepared_entry.row ||
        as_u32(identity_entry.at("column"),
               "spectral identity entry column") != prepared_entry.column ||
        required_string(identity_entry, "exact_entry") !=
            prepared_entry.multiplier.exact_identity ||
        as_i32(identity_entry.at("epsilon_shift"),
               "spectral identity epsilon shift") !=
            prepared_entry.multiplier.epsilon_shift ||
        as_u32(identity_entry.at("center_pole_order"),
               "spectral identity center-pole order") !=
            prepared_entry.multiplier.center_pole_order)
      throw std::invalid_argument(
          "prepared SCC spectral VInv entry differs from its exact structural identity");
  }
  if (transform.identity) {
    if (transform.matrix.entries.size() != dimension)
      throw std::invalid_argument(
          "identity SCC spectral VInv must contain exactly one diagonal entry per component");
    const auto one = ScalarTraits<Scalar>::one();
    const auto scalar_identical = [&](const Scalar& left,
                                      const Scalar& right) {
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        return acb_equal(left.raw(), right.raw());
      else
        return left == right;
    };
    for (std::uint32_t component = 0; component < dimension; ++component) {
      const auto& entry = transform.matrix.entries[component];
      if (entry.row != component || entry.column != component ||
          entry.multiplier.epsilon_shift != 0 ||
          entry.multiplier.center_pole_order != 0)
        throw std::invalid_argument(
            "identity SCC spectral VInv is not a shift-zero structural diagonal");
      for (std::size_t epsilon = 0;
           epsilon < entry.multiplier.kernels.size(); ++epsilon)
        for (std::size_t taylor = 0;
             taylor < entry.multiplier.kernels[epsilon].size(); ++taylor) {
          const auto& coefficient =
              entry.multiplier.kernels[epsilon][taylor];
          const bool expected_one = epsilon == 0 && taylor == 0;
          if ((expected_one && !scalar_identical(coefficient, one)) ||
              (!expected_one &&
               !ScalarTraits<Scalar>::is_zero(coefficient)))
            throw std::invalid_argument(
                "identity SCC spectral VInv contains a nonidentity retained coefficient");
        }
    }
  }
  return transform;
}

ExactJordanIndicialCertificate parse_exact_jordan_indicial_record(
    const json::value& raw, std::uint32_t expected_dimension) {
  const auto& object = as_object(
      raw, "exact affine-Jordan indicial certificate");
  require_exact_keys(object, {"schema", "dimension", "blocks"},
                     "exact affine-Jordan indicial certificate");
  if (required_string(object, "schema") !=
      "diffexp2-exact-affine-jordan-indicial-v1")
    throw std::invalid_argument(
        "unsupported exact affine-Jordan indicial certificate schema");
  ExactJordanIndicialCertificate certificate;
  certificate.dimension = as_u32(
      object.at("dimension"), "exact affine-Jordan dimension");
  if (certificate.dimension == 0 ||
      certificate.dimension != expected_dimension)
    throw std::invalid_argument(
        "exact affine-Jordan certificate dimension differs from its block");
  certificate.block_of_column.assign(
      certificate.dimension, std::numeric_limits<std::uint32_t>::max());
  certificate.position_in_block.assign(
      certificate.dimension, std::numeric_limits<std::uint32_t>::max());
  const auto& blocks = as_array(
      object.at("blocks"), "exact affine-Jordan blocks");
  if (blocks.empty())
    throw std::invalid_argument(
        "exact affine-Jordan certificate has no Jordan blocks");
  certificate.blocks.reserve(blocks.size());
  for (std::size_t index = 0; index < blocks.size(); ++index) {
    const auto& block = as_object(
        blocks[index], "exact affine-Jordan block");
    require_exact_keys(block, {"block", "columns", "a", "b"},
                       "exact affine-Jordan block");
    const auto block_index = as_u32(
        block.at("block"), "exact affine-Jordan block index");
    if (block_index != index)
      throw std::invalid_argument(
          "exact affine-Jordan blocks are not in deterministic order");
    std::vector<std::uint32_t> columns;
    for (const auto& raw_column : as_array(
             block.at("columns"), "exact affine-Jordan columns")) {
      const auto column = as_u32(
          raw_column, "exact affine-Jordan column");
      if (column >= certificate.dimension)
        throw std::invalid_argument(
            "exact affine-Jordan column is outside its dimension");
      columns.push_back(column);
    }
    if (columns.empty())
      throw std::invalid_argument(
          "exact affine-Jordan certificate contains an empty block");
    for (std::size_t position = 0; position < columns.size(); ++position) {
      const auto column = columns[position];
      if (certificate.block_of_column[column] !=
          std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(
            "exact affine-Jordan blocks contain a duplicate column");
      certificate.block_of_column[column] = block_index;
      certificate.position_in_block[column] =
          static_cast<std::uint32_t>(position);
    }
    certificate.blocks.push_back(ExactJordanBlockCertificate{
        block_index, std::move(columns),
        ExactAffineIndicialRoot{
            Rational(required_string(block, "a")),
            Rational(required_string(block, "b"))}});
  }
  if (std::any_of(
          certificate.block_of_column.begin(),
          certificate.block_of_column.end(), [](std::uint32_t value) {
            return value == std::numeric_limits<std::uint32_t>::max();
          }))
    throw std::invalid_argument(
        "exact affine-Jordan blocks do not partition their dimension");
  return certificate;
}

bool same_exact_jordan_indicial(
    const ExactJordanIndicialCertificate& left,
    const ExactJordanIndicialCertificate& right) {
  if (left.dimension != right.dimension ||
      left.block_of_column != right.block_of_column ||
      left.position_in_block != right.position_in_block ||
      left.blocks.size() != right.blocks.size())
    return false;
  for (std::size_t index = 0; index < left.blocks.size(); ++index) {
    const auto& a = left.blocks[index];
    const auto& b = right.blocks[index];
    if (a.block_index != b.block_index || a.columns != b.columns ||
        !(a.root == b.root))
      return false;
  }
  return true;
}

struct CompositeCouplingIdentity {
  std::uint32_t source_vertex = 0;
  std::uint32_t target_vertex = 0;
  std::string exact_original_entry;
  std::string exact_theta_entry;
  bool proven_zero = false;
};

template <typename Scalar>
struct CompositeSCCCoupling {
  std::uint32_t source_block = 0;
  std::uint32_t target_block = 0;
  std::string producer_identity;
  PreparedSparseLocalMultiplierMatrix<Scalar> matrix;
  std::vector<CompositeCouplingIdentity> identities;
};

template <typename Scalar>
LocalSolution<Scalar> cap_composite_public_local(
    const LocalSolution<Scalar>& input, std::int32_t complete_max,
    std::uint32_t taylor_complete_max, const ChartGeometry& parent_chart,
    const std::vector<Prescription>& parent_prescriptions,
    std::string checkpoint_identity) {
  validate_local_solution(input, false);
  if (!input.error.empty())
    throw std::invalid_argument(
        "native SCC delivery cannot discard an error envelope");
  if (input.epsilon.complete_max < complete_max ||
      input.epsilon.min_power > complete_max)
    throw std::invalid_argument(
        "native SCC work state cannot deliver requested epsilon maximum " +
        std::to_string(complete_max) + " from work window [" +
        std::to_string(input.epsilon.min_power) + "," +
        std::to_string(input.epsilon.complete_max) + "]");
  if (input.taylor_complete_max < taylor_complete_max)
    throw std::invalid_argument(
        "native SCC work state cannot deliver the requested Taylor order");

  LocalSolution<Scalar> output;
  output.chart = parent_chart;
  output.epsilon = {input.epsilon.min_power, complete_max};
  output.taylor_complete_max = taylor_complete_max;
  output.dimension = input.dimension;
  output.prescriptions = parent_prescriptions;
  output.checkpoint_identity = std::move(checkpoint_identity);
  output.sectors.reserve(input.sectors.size());
  for (const auto& sector : input.sectors) {
    LocalSector<Scalar> capped;
    capped.a = sector.a;
    capped.b = sector.b;
    capped.log_power = sector.log_power;
    capped.coefficients.assign(output.sector_size(),
                               ScalarTraits<Scalar>::zero());
    for (std::int64_t power = output.epsilon.min_power;
         power <= output.epsilon.complete_max; ++power) {
      const auto input_epsilon = static_cast<std::size_t>(
          power - input.epsilon.min_power);
      const auto output_epsilon = static_cast<std::size_t>(
          power - output.epsilon.min_power);
      for (std::size_t n = 0; n < output.taylor_width(); ++n)
        for (std::uint32_t component = 0; component < output.dimension;
             ++component)
          capped.coefficients[local_algebra_detail::flat_index(
              output_epsilon, n, component, output.taylor_width(),
              output.dimension)] = sector.coefficients[
                  local_algebra_detail::flat_index(
                      input_epsilon, n, component, input.taylor_width(),
                      input.dimension)];
    }
    output.sectors.push_back(std::move(capped));
  }
  validate_local_solution(output, false);
  return output;
}

template <typename Scalar>
SourceData<Scalar> local_solution_source_data(
    const LocalSolution<Scalar>& source, std::uint32_t nmax,
    std::uint32_t log_max, std::int32_t frame_base,
    std::uint32_t frame_width) {
  validate_local_solution(source, false);
  if (source.dimension == 0 || !source.error.empty())
    throw std::invalid_argument(
        "native SCC source injection requires an uncertified nonempty local");
  const auto frame_top_i64 = static_cast<std::int64_t>(frame_base) +
      frame_width - 1;
  if (frame_top_i64 > std::numeric_limits<std::int32_t>::max() ||
      source.epsilon.min_power < frame_base ||
      source.epsilon.complete_max > frame_top_i64)
    throw std::invalid_argument(
        "native SCC source window lies outside the target retained frame");
  if (source.taylor_complete_max < nmax)
    throw std::invalid_argument(
        "native SCC source has insufficient Taylor order for its target");
  const auto taylor_points = static_cast<std::size_t>(nmax) + 1;
  const auto log_points = static_cast<std::size_t>(log_max) + 1;
  if (frame_width == 0 ||
      taylor_points > std::numeric_limits<std::size_t>::max() / log_points)
    throw std::overflow_error("native SCC source point tensor size overflow");
  const auto points = taylor_points * log_points;
  if (points > std::numeric_limits<std::size_t>::max() / source.dimension)
    throw std::overflow_error("native SCC source validity size overflow");
  const auto component_points = points * source.dimension;
  if (component_points >
      std::numeric_limits<std::size_t>::max() / frame_width)
    throw std::overflow_error("native SCC source tensor size overflow");
  SourceData<Scalar> data;
  data.frames.assign(component_points * frame_width,
                     ScalarTraits<Scalar>::zero());
  data.validity.assign(component_points, kCompleteInfinity);
  data.present.assign(points, 0);
  for (const auto& sector : source.sectors) {
    if (sector.log_power > log_max)
      throw std::invalid_argument(
          "native SCC source log sector exceeds the target run depth");
    for (std::uint32_t n = 0; n <= nmax; ++n) {
      const auto point = static_cast<std::size_t>(n) * log_points +
          sector.log_power;
      if (data.present[point])
        throw std::invalid_argument(
            "native SCC source contains duplicate exact log sectors");
      data.present[point] = 1;
      // Coefficients above CompleteMax are unknown, not certified zeros.  The
      // dense work-frame storage remains zero there while this finite validity
      // bound propagates independently for every target-block component.
      for (std::uint32_t component = 0;
           component < source.dimension; ++component) {
        const auto component_point = point * source.dimension + component;
        data.validity[component_point] = source.epsilon.complete_max;
        for (std::int64_t power = source.epsilon.min_power;
             power <= source.epsilon.complete_max; ++power) {
          const auto input_epsilon = static_cast<std::size_t>(
              power - source.epsilon.min_power);
          const auto output_epsilon = static_cast<std::size_t>(
              power - frame_base);
          data.frames[component_point * frame_width + output_epsilon] =
              sector.coefficients[local_algebra_detail::flat_index(
                  input_epsilon, n, component, source.taylor_width(),
                  source.dimension)];
        }
      }
    }
  }
  return data;
}

bool exact_nonnegative_integer(const Rational& value, bool include_zero) {
  if (value.sign() < 0 || (!include_zero && value.is_zero())) return false;
  return value.str().find('/') == std::string::npos;
}

std::uint32_t exact_log_ceiling(
    const ExactJordanIndicialCertificate& indicial,
    const Rational& a, const Rational& b, std::uint32_t base,
    bool include_zero_offset) {
  std::uint64_t result = base;
  for (const auto& block : indicial.blocks) {
    if (!(block.root.b == b)) continue;
    const auto offset = block.root.a - a;
    if (exact_nonnegative_integer(offset, include_zero_offset))
      result += block.size();
  }
  if (result > std::numeric_limits<std::uint32_t>::max())
    throw RecurrenceError(
        "E5", "derived exact Jordan log ceiling exceeds uint32 range");
  return static_cast<std::uint32_t>(result);
}

template <typename Scalar>
json::object exact_derived_run(
    const json::object& prototype, const PreparedChart<Scalar>& chart,
    const ExactJordanIndicialCertificate& indicial,
    const Rational& a, const Rational& b, std::uint32_t base_log,
    bool homogeneous, std::optional<std::uint32_t> seed_component) {
  const auto nmax = as_u32(prototype.at("nmax"),
                           "derived pseudo-compensation Taylor order");
  const auto dimension = chart.dimension();
  const auto frame_base = chart.frame_base();
  const auto frame_width = chart.frame_width();
  const auto frame_top_i64 = static_cast<std::int64_t>(frame_base) +
                             frame_width - 1;
  if (frame_top_i64 < std::numeric_limits<std::int32_t>::min() ||
      frame_top_i64 > std::numeric_limits<std::int32_t>::max())
    throw RecurrenceError(
        "E5", "derived pseudo-compensation frame exceeds int32 range");
  const auto frame_top = static_cast<std::int32_t>(frame_top_i64);

  std::uint32_t position = 0;
  if (homogeneous) {
    if (!seed_component.has_value() || *seed_component >= dimension)
      throw RecurrenceError(
          "E5", "derived pseudo-compensation seed component is out of range");
    position = indicial.position_in_block[*seed_component];
    base_log = std::max(base_log, position);
  } else if (seed_component.has_value()) {
    throw RecurrenceError(
        "E5", "derived particular run unexpectedly carries a seed component");
  }
  const auto log_max = exact_log_ceiling(
      indicial, a, b, base_log, !homogeneous);

  json::array a_shifts;
  json::array schedule;
  a_shifts.reserve(static_cast<std::size_t>(nmax) + 1);
  schedule.reserve(static_cast<std::size_t>(nmax) + 1);
  for (std::uint32_t n = 0; n <= nmax; ++n) {
    const auto a_n = a + Rational(std::to_string(n));
    a_shifts.push_back(json::string(a_n.str()));
    json::array row;
    row.reserve(indicial.blocks.size());
    for (const auto& block : indicial.blocks) {
      const auto d_a = a_n - block.root.a;
      const auto d_b = b - block.root.b;
      const auto step = singular_indicial_detail::classify_step(d_a, d_b);
      row.push_back(json::object{
          {"case", singular_indicial_detail::step_name(step)},
          {"da", d_a.str()}, {"db", d_b.str()}});
    }
    schedule.push_back(std::move(row));
  }

  json::array initial;
  json::array validity;
  const auto points = static_cast<std::size_t>(log_max) + 1;
  if (points > std::numeric_limits<std::size_t>::max() / dimension ||
      points * dimension >
          std::numeric_limits<std::size_t>::max() / frame_width)
    throw RecurrenceError(
        "E5", "derived pseudo-compensation seed tensor size overflows");
  initial.reserve(points * dimension * frame_width);
  validity.reserve(points * dimension);
  for (std::uint32_t log = 0; log <= log_max; ++log) {
    std::optional<std::uint32_t> expected_component;
    std::optional<std::size_t> expected_epsilon;
    if (homogeneous && log <= position) {
      const auto block_index = indicial.block_of_column[*seed_component];
      const auto& block = indicial.blocks[block_index];
      expected_component = block.columns[position - log];
      const auto epsilon_i64 = -static_cast<std::int64_t>(log) - frame_base;
      if (epsilon_i64 < 0 || epsilon_i64 >= frame_width)
        throw RecurrenceError(
            "E4", "derived canonical Jordan seed exceeds the retained lower epsilon frame",
            frame_base, -static_cast<std::int32_t>(log));
      expected_epsilon = static_cast<std::size_t>(epsilon_i64);
    }
    for (std::uint32_t component = 0; component < dimension; ++component) {
      for (std::uint32_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const bool unit = expected_component.has_value() &&
                          component == *expected_component &&
                          epsilon == *expected_epsilon;
        initial.push_back(unit ? "1" : "0");
      }
      validity.push_back(homogeneous ? json::value(frame_top)
                                     : json::value(nullptr));
    }
  }

  auto run = prototype;
  run["p"] = log_max;
  run["has_initial"] = homogeneous;
  // A homogeneous CASE-P target is a composite-owned mathematical object,
  // not a replay of whichever public column happened to discover it first.
  // Submitted seed templates disable the lower-frame probe while derived
  // particular templates enable it, so inheriting this diagnostic policy
  // made the same (block, component) target acquire two cache contracts.
  // Probe mode only changes when lower-edge cancellation is checked (matrix
  // assembly before a negative epsilon shift); the enabled form is the
  // stronger, cancellation-aware evaluation and is deterministic for every
  // cached homogeneous target. Particular runs retain their caller policy.
  run["adaptive_probe"] = homogeneous
      ? json::value(true)
      : prototype.at("adaptive_probe");
  run["a_target"] = a.str();
  run["b_target"] = b.str();
  run["a_shift_min"] = 0;
  run["a_shifts"] = std::move(a_shifts);
  run["schedule"] = std::move(schedule);
  run["initial"] = std::move(initial);
  run["initial_validity"] = std::move(validity);
  run["source"] = nullptr;
  run["return_u"] = false;
  return run;
}

json::object exact_derived_metadata(
    const json::object& prototype, const Rational& a, const Rational& b,
    std::uint32_t log_max, const std::string& suffix) {
  auto metadata = prototype;
  auto& tag = metadata.at("tag").as_object();
  tag["a"] = json::object{{"domain", "rational"},
                           {"canonical", a.str()}};
  tag["b"] = json::object{{"domain", "rational"},
                           {"canonical", b.str()}};
  tag["p"] = json::object{{"domain", "integer"},
                           {"canonical", std::to_string(log_max)}};
  metadata["checkpoint_identity"] =
      required_string(prototype, "checkpoint_identity") + suffix;
  return metadata;
}

struct ExactFormalKey {
  std::string t_power;
  std::int32_t epsilon_power = 0;
  std::uint32_t log_power = 0;
  std::uint32_t component = 0;

  friend bool operator<(const ExactFormalKey& left,
                        const ExactFormalKey& right) {
    return std::tie(left.t_power, left.epsilon_power, left.log_power,
                    left.component) <
           std::tie(right.t_power, right.epsilon_power, right.log_power,
                    right.component);
  }
};

// Expand an exact local slab in the formal basis
//
//   t^(a+n) eps^K Log(t)^p
//
// using t^(b eps) = Sum_j (b eps Log(t))^j/j!.  This is the
// Rational-domain CASE-P certificate: neither an Acb midpoint nor a
// tolerance can decide whether a polar coefficient cancels.  The optional
// t ceiling removes only the unmatched high-Taylor tail introduced when a
// target root starts n>0 orders above the source; every overlapping stored
// coefficient remains an exact proof obligation.
void add_exact_formal_below(
    const LocalSolution<Rational>& solution, std::int32_t exclusive_top,
    const std::optional<Rational>& maximum_t_power,
    std::map<ExactFormalKey, Rational>& coefficients) {
  validate_local_solution(solution, false);
  for (const auto& sector : solution.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw RecurrenceError(
          "E5", "exact pseudo-compensation certificate requires rational sector tags");
    const Rational a(sector.a.canonical);
    const Rational b(sector.b.canonical);
    Rational log_normalization(1);
    for (std::uint32_t divisor = 2; divisor <= sector.log_power; ++divisor)
      log_normalization =
          log_normalization / Rational(std::to_string(divisor));
    for (std::size_t n = 0; n < solution.taylor_width(); ++n) {
      const auto total_t_power = a + Rational(std::to_string(n));
      if (maximum_t_power.has_value() &&
          total_t_power > *maximum_t_power)
        continue;
      for (std::int64_t epsilon = solution.epsilon.min_power;
           epsilon <= solution.epsilon.complete_max; ++epsilon) {
        const auto epsilon_index = static_cast<std::size_t>(
            epsilon - solution.epsilon.min_power);
        const auto base_power_i64 = epsilon + sector.log_power;
        if (base_power_i64 >= exclusive_top) continue;
        if (base_power_i64 < std::numeric_limits<std::int32_t>::min())
          throw RecurrenceError(
              "E5", "exact pseudo-compensation epsilon power underflows int32");
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component) {
          const auto& value = sector.coefficients[
              local_algebra_detail::flat_index(
                  epsilon_index, n, component, solution.taylor_width(),
                  solution.dimension)];
          if (value.is_zero()) continue;
          Rational exponential_factor(1);
          for (std::uint64_t j = 0;
               base_power_i64 + static_cast<std::int64_t>(j) <
                   exclusive_top;
               ++j) {
            if (j > std::numeric_limits<std::uint32_t>::max() -
                        sector.log_power)
              throw RecurrenceError(
                  "E5", "exact pseudo-compensation log degree overflows uint32");
            ExactFormalKey key{
                total_t_power.str(),
                static_cast<std::int32_t>(base_power_i64 +
                                          static_cast<std::int64_t>(j)),
                sector.log_power + static_cast<std::uint32_t>(j),
                component};
            auto found = coefficients.try_emplace(key, Rational(0)).first;
            found->second += value * log_normalization * exponential_factor;
            if (found->second.is_zero()) coefficients.erase(found);
            if (b.is_zero()) break;
            if (j == std::numeric_limits<std::uint32_t>::max())
              throw RecurrenceError(
                  "E5", "exact pseudo-compensation exponential order overflows uint32");
            exponential_factor = exponential_factor * b /
                Rational(std::to_string(j + 1));
          }
        }
      }
    }
  }
}

// Numeric counterpart of the exact formal expansion above.  Exact Rational
// tags still determine the formal monomials; only their coefficients live in
// Acb.  This lets an enclosure prove CASE-P cancellation at the configured
// chop floor without promoting the whole recurrence to arbitrary-precision
// rational arithmetic.
void add_acb_formal_below(
    const LocalSolution<ComplexBall>& solution,
    std::int32_t exclusive_top,
    const std::optional<Rational>& maximum_t_power,
    std::map<ExactFormalKey, ComplexBall>& coefficients) {
  validate_local_solution(solution, false);
  for (const auto& sector : solution.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw RecurrenceError(
          "E5", "Acb pseudo-compensation requires exact rational sector tags");
    const Rational a(sector.a.canonical);
    const Rational b(sector.b.canonical);
    Rational exact_log_normalization(1);
    for (std::uint32_t divisor = 2; divisor <= sector.log_power; ++divisor)
      exact_log_normalization =
          exact_log_normalization / Rational(std::to_string(divisor));
    const auto log_normalization =
        ComplexBall::from_strings(exact_log_normalization.str());
    const auto b_ball = ComplexBall::from_strings(b.str());
    for (std::size_t n = 0; n < solution.taylor_width(); ++n) {
      const auto total_t_power = a + Rational(std::to_string(n));
      if (maximum_t_power.has_value() &&
          total_t_power > *maximum_t_power)
        continue;
      for (std::int64_t epsilon = solution.epsilon.min_power;
           epsilon <= solution.epsilon.complete_max; ++epsilon) {
        const auto epsilon_index = static_cast<std::size_t>(
            epsilon - solution.epsilon.min_power);
        const auto base_power_i64 = epsilon + sector.log_power;
        if (base_power_i64 >= exclusive_top) continue;
        if (base_power_i64 < std::numeric_limits<std::int32_t>::min())
          throw RecurrenceError(
              "E5", "Acb pseudo-compensation epsilon power underflows int32");
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component) {
          const auto& value = sector.coefficients[
              local_algebra_detail::flat_index(
                  epsilon_index, n, component, solution.taylor_width(),
                  solution.dimension)];
          if (value.is_zero()) continue;
          auto exponential_factor = ComplexBall::from_strings("1");
          for (std::uint64_t j = 0;
               base_power_i64 + static_cast<std::int64_t>(j) <
                   exclusive_top;
               ++j) {
            if (j > std::numeric_limits<std::uint32_t>::max() -
                        sector.log_power)
              throw RecurrenceError(
                  "E5", "Acb pseudo-compensation log degree overflows uint32");
            ExactFormalKey key{
                total_t_power.str(),
                static_cast<std::int32_t>(base_power_i64 +
                                          static_cast<std::int64_t>(j)),
                sector.log_power + static_cast<std::uint32_t>(j),
                component};
            auto found = coefficients.try_emplace(
                key, ScalarTraits<ComplexBall>::zero()).first;
            found->second += value * log_normalization * exponential_factor;
            if (b.is_zero()) break;
            if (j == std::numeric_limits<std::uint32_t>::max())
              throw RecurrenceError(
                  "E5", "Acb pseudo-compensation exponential order overflows uint32");
            exponential_factor *= b_ball /
                ComplexBall(static_cast<long>(j + 1));
          }
        }
      }
    }
  }
}

void canonicalize_acb_formal_coefficients(
    std::map<ExactFormalKey, ComplexBall>& coefficients,
    std::int32_t chop_digits, const char* context) {
  for (auto iterator = coefficients.begin();
       iterator != coefficients.end();) {
    auto canonical = ScalarTraits<ComplexBall>::canonicalized(
        iterator->second, chop_digits);
    if (canonical.is_zero()) {
      iterator = coefficients.erase(iterator);
      continue;
    }
    if (canonical.contains_zero())
      throw RecurrenceError(
          "E5", std::string(context) +
              " is numerically ambiguous; requires the exact Rational shadow");
    iterator->second = std::move(canonical);
    ++iterator;
  }
}

void certify_acb_pseudo_value_floor(
    const LocalSolution<ComplexBall>& solution,
    std::int32_t allowed_floor, const Rational& source_a,
    std::uint32_t source_taylor_max, std::int32_t chop_digits) {
  const auto checked_top = std::min<std::int32_t>(0, allowed_floor);
  std::map<ExactFormalKey, ComplexBall> coefficients;
  add_acb_formal_below(
      solution, checked_top,
      source_a + Rational(std::to_string(source_taylor_max)), coefficients);
  canonicalize_acb_formal_coefficients(
      coefficients, chop_digits, "native Acb CASE-P value-floor cancellation");
  if (coefficients.empty()) return;
  const auto& witness = coefficients.begin()->first;
  // The Acb owner encloses the finite-precision coefficients serialized by
  // Wolfram, not the original exact rational system.  A certified remnant
  // therefore proves that this numeric owner cannot authorize the required
  // cancellation; it does not prove that the exact CASE-P identity fails.
  // Name the existing one-SCC Rational-shadow fallback explicitly so the
  // bridge can retry without weakening or chopping this certificate.
  throw RecurrenceError(
      "E5", "Acb CASE-P compensation leaves a certified value pole below the work frame at eps^" +
                std::to_string(witness.epsilon_power) + ", t_power=" +
                witness.t_power + ", log_power=" +
                std::to_string(witness.log_power) + ", component=" +
                std::to_string(witness.component) +
                "; requires the exact Rational shadow",
      solution.epsilon.min_power, witness.epsilon_power);
}

EpsilonLatticeSaturationResult<Rational>
rational_shadow_formal_saturation(
    const std::vector<std::shared_ptr<const RationalShadowColumnWitness>>&
        witnesses,
    EpsilonWindow window, const std::string& context) {
  const auto dimension = witnesses.size();
  if (dimension == 0 || window.min_power > window.complete_max)
    throw std::invalid_argument(
        context + ": empty Rational-shadow formal saturation request");
  using FormalRow = std::tuple<std::string, std::uint32_t, std::uint32_t>;
  std::vector<std::map<ExactFormalKey, Rational>> columns(dimension);
  std::set<FormalRow> row_set;
  const auto exclusive_top_i64 =
      static_cast<std::int64_t>(window.complete_max) + 1;
  if (exclusive_top_i64 > std::numeric_limits<std::int32_t>::max())
    throw std::overflow_error(
        "Rational-shadow formal saturation top overflows int32");
  for (std::size_t column = 0; column < dimension; ++column) {
    if (!witnesses[column] || !witnesses[column]->solution)
      throw std::invalid_argument(
          context + ": Rational-shadow formal column is absent");
    add_exact_formal_below(
        *witnesses[column]->solution,
        static_cast<std::int32_t>(exclusive_top_i64), std::nullopt,
        columns[column]);
    for (const auto& [key, coefficient] : columns[column]) {
      if (coefficient.is_zero()) continue;
      if (key.epsilon_power < window.min_power)
        throw MatchingArithmeticError(
            MatchingArithmeticErrorCode::InsufficientCompleteWindow,
            context +
                ": exact formal shadow contains content below the matching frame",
            std::nullopt, column, key.epsilon_power);
      if (key.epsilon_power <= window.complete_max)
        row_set.emplace(key.t_power, key.log_power, key.component);
    }
  }
  if (row_set.size() < dimension)
    throw MatchingArithmeticError(
        MatchingArithmeticErrorCode::SingularOrIncompleteSystem,
        context + ": exact formal shadow has fewer rows than basis columns");
  const std::vector<FormalRow> rows(row_set.begin(), row_set.end());
  FiniteLaurentMatrix<Rational> rectangular;
  rectangular.reserve(rows.size());
  for (const auto& [t_power, log_power, component] : rows) {
    FiniteLaurentVector<Rational> row;
    row.reserve(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      std::vector<Rational> coefficients(window.width(), Rational(0));
      for (std::int64_t power = window.min_power;
           power <= window.complete_max; ++power) {
        const auto found = columns[column].find(ExactFormalKey{
            t_power, static_cast<std::int32_t>(power), log_power,
            component});
        if (found != columns[column].end())
          coefficients[static_cast<std::size_t>(
              power - window.min_power)] = found->second;
      }
      row.emplace_back(window, std::move(coefficients));
    }
    rectangular.push_back(std::move(row));
  }

  return saturate_rectangular_finite_laurent_basis(
      rectangular, context + ":full-formal-module");
}

void certify_exact_pseudo_value_floor(
    const LocalSolution<Rational>& solution, std::int32_t allowed_floor,
    const Rational& source_a, std::uint32_t source_taylor_max) {
  // CASE-P compensation may introduce genuine dimensional poles, including
  // poles absent from the inhomogeneous source.  The parent SCC residual is
  // the proof of correctness; this local guard only prevents compensation
  // from silently escaping the bounded work frame.
  const auto checked_top = std::min<std::int32_t>(0, allowed_floor);
  std::map<ExactFormalKey, Rational> coefficients;
  add_exact_formal_below(
      solution, checked_top,
      source_a + Rational(std::to_string(source_taylor_max)), coefficients);
  if (coefficients.empty()) return;
  const auto& witness = coefficients.begin()->first;
  std::string sector_witnesses;
  for (const auto& sector : solution.sectors) {
    auto single = solution;
    single.sectors = {sector};
    std::map<ExactFormalKey, Rational> part;
    add_exact_formal_below(
        single, checked_top,
        source_a + Rational(std::to_string(source_taylor_max)), part);
    const auto found = part.find(witness);
    if (found != part.end()) {
      if (!sector_witnesses.empty()) sector_witnesses += ";";
      sector_witnesses += "a=" + sector.a.canonical + ",b=" +
          sector.b.canonical + ",p=" +
          std::to_string(sector.log_power) + ":" + found->second.str();
    }
  }
  throw RecurrenceError(
      "E5", "exact CASE-P compensation leaves a value pole below the work frame at eps^" +
                std::to_string(witness.epsilon_power) + ", t_power=" +
                witness.t_power + ", log_power=" +
                std::to_string(witness.log_power) + ", component=" +
                std::to_string(witness.component) + ", coefficient=" +
                coefficients.begin()->second.str() + ", sectors=[" +
                sector_witnesses + "]",
      solution.epsilon.min_power, witness.epsilon_power);
}

std::string exact_local_sector_summary(
    const LocalSolution<Rational>& solution) {
  validate_local_solution(solution, false);
  std::string result = "{";
  constexpr std::size_t kSectorLimit = 32;
  for (std::size_t sector_index = 0;
       sector_index < solution.sectors.size() &&
       sector_index < kSectorLimit; ++sector_index) {
    if (sector_index != 0) result += ";";
    const auto& sector = solution.sectors[sector_index];
    result += "a=" + sector.a.canonical + ",b=" + sector.b.canonical +
        ",p=" + std::to_string(sector.log_power);
    std::size_t nonzero = 0;
    std::string first;
    for (std::int32_t power = solution.epsilon.min_power;
         power <= solution.epsilon.complete_max; ++power) {
      const auto epsilon_index = static_cast<std::size_t>(
          power - solution.epsilon.min_power);
      for (std::uint32_t n = 0;
           n <= solution.taylor_complete_max; ++n)
        for (std::uint32_t component = 0;
             component < solution.dimension; ++component) {
          const auto& coefficient = sector.coefficients[
              local_algebra_detail::flat_index(
                  epsilon_index, n, component, solution.taylor_width(),
                  solution.dimension)];
          if (coefficient.is_zero()) continue;
          ++nonzero;
          if (first.empty())
            first = "@" + std::to_string(n) + "," +
                std::to_string(power) + "," +
                std::to_string(component) + "=" +
                scc_completeness_detail::bounded_rational_witness(
                    coefficient);
        }
    }
    result += ",nz=" + std::to_string(nonzero);
    if (!first.empty()) result += first;
  }
  if (solution.sectors.size() > kSectorLimit)
    result += ";...[sectors=" +
        std::to_string(solution.sectors.size()) + "]";
  result += "}";
  return result;
}

template <typename Scalar>
std::vector<LocalSolution<Scalar>> split_exact_rational_tags(
    const LocalSolution<Scalar>& source, const std::string& identity) {
  validate_local_solution(source, false);
  std::map<std::pair<std::string, std::string>,
           std::vector<LocalSector<Scalar>>> grouped;
  for (const auto& sector : source.sectors) {
    if (sector.a.domain != ExactDomain::Rational ||
        sector.b.domain != ExactDomain::Rational)
      throw RecurrenceError(
          "E5", "native SCC pseudo propagation requires exact rational source tags");
    grouped[{Rational(sector.a.canonical).str(),
             Rational(sector.b.canonical).str()}].push_back(sector);
  }
  std::vector<LocalSolution<Scalar>> result;
  result.reserve(grouped.size());
  for (auto& [tag, sectors] : grouped) {
    auto group = source;
    group.sectors = std::move(sectors);
    group.checkpoint_identity = identity + ":tag:" + tag.first + ":" +
                                tag.second;
    result.push_back(canonicalize_identical_local_sectors(std::move(group)));
  }
  return result;
}

template <typename Scalar>
struct CachedPseudoTarget {
  LocalSolution<Scalar> solution;
  NativeLocalDiagnostics diagnostics;
};

template <typename Scalar>
class PseudoTargetCache {
 public:
  using Key = std::pair<std::uint32_t, std::uint32_t>;
  using Core = detail::ImmutableRecursiveCache<Key, CachedPseudoTarget<Scalar>>;
  using Lookup = Core::Lookup;
  using Stats = Core::Stats;

  explicit PseudoTargetCache(std::string owner_identity)
      : owner_identity_(std::move(owner_identity)) {
    if (owner_identity_.empty())
      throw std::invalid_argument(
          "exact CASE-P target cache requires an immutable owner identity");
  }

  std::string target_checkpoint_identity(std::uint32_t block,
                                         std::uint32_t component) const {
    return owner_identity_ + ":casep-homogeneous:block:" +
        std::to_string(block) + ":component:" +
        std::to_string(component);
  }

  template <typename Builder>
  Lookup get_or_build(std::uint32_t block, std::uint32_t component,
                      const std::string& contract, Builder&& builder) {
    try {
      return cache_.get_or_build(
          {block, component}, contract, std::forward<Builder>(builder));
    } catch (const detail::ImmutableCacheContractError&) {
      throw RecurrenceError(
          "E5", "exact CASE-P homogeneous target cache contract changed within one retained composite");
    } catch (const detail::ImmutableCacheCycleError&) {
      throw RecurrenceError(
          "E5", "exact CASE-P compensation dependency is cyclic; the retained family ordering is not well founded");
    }
  }

  Stats stats() const { return cache_.stats(); }

 private:
  std::string owner_identity_;
  Core cache_;
};

template <typename Scalar>
class PseudoCompensator {
 public:
  PseudoCompensator(PreparedChart<Scalar>& chart,
                    const ExactJordanIndicialCertificate& indicial,
                    std::uint32_t block_index,
                    std::string identity,
                    PseudoTargetCache<Scalar>& target_cache)
      : chart_(chart), block_index_(block_index),
        indicial_(indicial), identity_(std::move(identity)),
        target_cache_(target_cache) {}

  NativeLocalRun<Scalar> solve(
      const json::object& run, const json::object& metadata,
      std::optional<SourceData<Scalar>> source,
      std::int32_t allowed_value_floor) {
    NativeLocalDiagnostics diagnostics;
    auto result = solve_impl(run, metadata, std::move(source),
                             allowed_value_floor, diagnostics);
    result.diagnostics = diagnostics;
    return result;
  }

 private:
  static void accumulate(NativeLocalDiagnostics& total,
                         const NativeLocalDiagnostics& current) {
    total.top_valid = std::min(total.top_valid, current.top_valid);
    total.parse_ms += current.parse_ms;
    total.kernel_ms += current.kernel_ms;
    total.requires_parent_completeness_certificate =
        total.requires_parent_completeness_certificate ||
        current.requires_parent_completeness_certificate;
    total.pseudo_hits += current.pseudo_hits;
    total.pseudo_compensations += current.pseudo_compensations;
    total.max_pseudo_depth = std::max(
        total.max_pseudo_depth, current.max_pseudo_depth);
    total.pseudo_value_certified =
        total.pseudo_value_certified && current.pseudo_value_certified;
  }

  std::optional<PreparedRationalTaylorMultiplier<Scalar>> polar_weight(
      const PseudoHit<Scalar>& hit, std::size_t row,
      const LocalSolution<Scalar>& target) const {
    if (row >= hit.gamma_frames.size() || row >= hit.gamma_validity.size())
      throw RecurrenceError(
          "E5", "CASE-P hit has inconsistent gamma frame dimensions");
    if (hit.gamma_validity[row] != kCompleteInfinity &&
        hit.gamma_validity[row] < -1)
      throw RecurrenceError(
          "E4", "CASE-P polar frame is not complete through eps^-1",
          chart_.frame_base(), hit.gamma_validity[row]);
    const auto& gamma = hit.gamma_frames[row];
    if (gamma.size() != chart_.frame_width())
      throw RecurrenceError(
          "E5", "CASE-P gamma frame differs from its retained chart width");
    std::optional<std::int32_t> minimum;
    for (std::size_t index = 0; index < gamma.size(); ++index) {
      const auto power = chart_.frame_base() +
                         static_cast<std::int32_t>(index);
      if (power >= 0) break;
      auto coefficient = ScalarTraits<Scalar>::canonicalized(
          gamma[index], chart_.chop_digits());
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        if (!coefficient.is_zero() && coefficient.contains_zero())
          throw RecurrenceError(
              "E5", "native Acb CASE-P polar weight is numerically ambiguous; requires the exact Rational shadow");
      if (!ScalarTraits<Scalar>::is_zero(coefficient)) {
        minimum = power;
        break;
      }
    }
    if (!minimum.has_value()) return std::nullopt;
    PreparedRationalTaylorMultiplier<Scalar> multiplier;
    multiplier.epsilon_shift = *minimum;
    multiplier.center_pole_order = 0;
    multiplier.exact_identity = identity_ + ":casep:" +
        std::to_string(hit.n) + ":" + std::to_string(row);
    multiplier.kernels.assign(
        target.epsilon.width(),
        std::vector<Scalar>(target.taylor_width(),
                            ScalarTraits<Scalar>::zero()));
    for (std::size_t kernel = 0; kernel < multiplier.kernels.size();
         ++kernel) {
      const auto power = static_cast<std::int64_t>(*minimum) +
                         static_cast<std::int64_t>(kernel);
      if (power >= 0) break;
      const auto gamma_index = power - chart_.frame_base();
      if (gamma_index < 0 ||
          gamma_index >= static_cast<std::int64_t>(gamma.size()))
        throw RecurrenceError(
            "E4", "CASE-P polar weight lies outside its retained gamma frame",
            chart_.frame_base(), static_cast<std::int32_t>(power));
      auto coefficient = ScalarTraits<Scalar>::canonicalized(
          gamma[static_cast<std::size_t>(gamma_index)],
          chart_.chop_digits());
      if constexpr (std::is_same_v<Scalar, ComplexBall>)
        if (!coefficient.is_zero() && coefficient.contains_zero())
          throw RecurrenceError(
              "E5", "native Acb CASE-P polar coefficient is numerically ambiguous; requires the exact Rational shadow");
      multiplier.kernels[kernel][0] = -coefficient;
    }
    return multiplier;
  }

  const LocalSolution<Scalar>& homogeneous_target(
      std::uint32_t component, const json::object& prototype_run,
      const json::object& prototype_metadata,
      std::int32_t allowed_value_floor,
      NativeLocalDiagnostics& diagnostics) {
    if (const auto found = homogeneous_cache_.find(component);
        found != homogeneous_cache_.end())
      return found->second->solution;
    if (!active_components_.insert(component).second)
      throw RecurrenceError(
          "E5", "exact CASE-P compensation dependency is cyclic; the retained family ordering is not well founded");
    try {
      if (component >= indicial_.dimension)
        throw RecurrenceError(
            "E5", "CASE-P target component has no retained exact Jordan root");
      const auto& block =
          indicial_.blocks[indicial_.block_of_column[component]];
      auto run = exact_derived_run(
          prototype_run, chart_, indicial_, block.root.a, block.root.b,
          indicial_.position_in_block[component], true, component);
      auto metadata = exact_derived_metadata(
          prototype_metadata, block.root.a, block.root.b,
          as_u32(run.at("p"), "derived homogeneous log maximum"),
          ":casep-target:" + std::to_string(component));
      auto contract_metadata = metadata;
      contract_metadata.erase("checkpoint_identity");
      const auto contract = json::serialize(canonical_json_value(
          json::object{{"schema", "diffexp2-casep-target-cache-v1"},
                       {"block", block_index_},
                       {"component", component},
                       {"run", run},
                       {"metadata", std::move(contract_metadata)}}));
      const auto target_identity =
          target_cache_.target_checkpoint_identity(block_index_, component);
      metadata["checkpoint_identity"] = target_identity;
      auto lookup = target_cache_.get_or_build(
          block_index_, component, contract, [&] {
            NativeLocalDiagnostics target_diagnostics;
            PseudoCompensator<Scalar> target_builder(
                chart_, indicial_, block_index_, target_identity,
                target_cache_);
            auto solved = target_builder.solve_impl(
                run, metadata, std::nullopt, allowed_value_floor,
                target_diagnostics);
            solved.solution.checkpoint_identity = target_identity;
            validate_local_solution(solved.solution, false);
            return CachedPseudoTarget<Scalar>{
                std::move(solved.solution), target_diagnostics};
          });
      accumulate(diagnostics, lookup.value->diagnostics);
      auto [stored, inserted] = homogeneous_cache_.emplace(
          component, std::move(lookup.value));
      if (!inserted)
        throw std::logic_error(
            "CASE-P homogeneous target cache insertion failed");
      active_components_.erase(component);
      return stored->second->solution;
    } catch (...) {
      active_components_.erase(component);
      throw;
    }
  }

  NativeLocalRun<Scalar> solve_impl(
      const json::object& run, const json::object& metadata,
      std::optional<SourceData<Scalar>> source,
      std::int32_t allowed_value_floor,
      NativeLocalDiagnostics& diagnostics) {
    auto raw = source.has_value()
        ? chart_.solve_native_with_source(
              run, metadata, std::move(*source))
        : chart_.solve_native(run, metadata);
    accumulate(diagnostics, raw.diagnostics);
    if (raw.pseudo_hits.empty()) return raw;

    const auto source_a = parse_scalar<Rational>(run.at("a_target"));
    const auto nmax = as_u32(run.at("nmax"),
                             "CASE-P certificate Taylor order");
    diagnostics.pseudo_hits += raw.pseudo_hits.size();
    std::vector<LocalSolution<Scalar>> terms;
    terms.push_back(std::move(raw.solution));
    for (const auto& hit : raw.pseudo_hits) {
      diagnostics.max_pseudo_depth = std::max<std::uint32_t>(
          diagnostics.max_pseudo_depth,
          static_cast<std::uint32_t>(hit.columns.size()));
      if (hit.columns.size() != hit.gamma_frames.size() ||
          hit.columns.size() != hit.gamma_validity.size())
        throw RecurrenceError(
            "E5", "CASE-P hit target and gamma dimensions disagree");
      for (std::size_t row = 0; row < hit.columns.size(); ++row) {
        const auto& target = homogeneous_target(
            hit.columns[row], run, metadata, allowed_value_floor,
            diagnostics);
        auto multiplier = polar_weight(hit, row, target);
        if (!multiplier.has_value()) continue;
        auto product = multiply_prepared_rational(
            target, *multiplier,
            identity_ + ":casep-product:" + std::to_string(hit.n) + ":" +
                std::to_string(hit.columns[row]));
        terms.push_back(std::move(product));
        ++diagnostics.pseudo_compensations;
      }
    }
    auto compensated = terms.size() == 1
        ? std::move(terms.front())
        : combine_local_solutions(
              terms, identity_ + ":casep-compensated");
    if constexpr (std::is_same_v<Scalar, Rational>) {
      certify_exact_pseudo_value_floor(
          compensated, allowed_value_floor, source_a, nmax);
    } else {
      static_assert(std::is_same_v<Scalar, ComplexBall>);
      certify_acb_pseudo_value_floor(
          compensated, allowed_value_floor, source_a, nmax,
          chart_.chop_digits());
    }
    raw.solution = std::move(compensated);
    raw.pseudo_hits.clear();
    return raw;
  }

  PreparedChart<Scalar>& chart_;
  std::uint32_t block_index_ = 0;
  const ExactJordanIndicialCertificate& indicial_;
  std::string identity_;
  PseudoTargetCache<Scalar>& target_cache_;
  std::map<std::uint32_t,
           std::shared_ptr<const CachedPseudoTarget<Scalar>>>
      homogeneous_cache_;
  std::set<std::uint32_t> active_components_;
};

template <typename Scalar>
class CompositeSCCChart final : public CompositeSCCChartBase {
 public:
  CompositeSCCChart(std::string handle, std::string key,
                    std::string exact_identity, std::string signature,
                    std::string rational_shadow_identity,
                    std::uint32_t dimension, SCCCertificate graph,
                    std::string exact_system_record,
                    std::string exact_theta_record,
                    std::string geometry_record,
                    RetainedCompositeGeometry retained_geometry,
                    CompositeWorkContract work,
                    std::vector<CompositeSCCBlock<Scalar>> blocks,
                    std::vector<CompositeSCCCoupling<Scalar>> couplings,
                    std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
                        physical_equation)
      : CompositeSCCChartBase(std::move(handle), std::move(key),
                              std::move(exact_identity),
                              std::move(signature),
                              std::move(rational_shadow_identity)),
        dimension_(dimension), graph_(std::move(graph)),
        exact_system_record_(std::move(exact_system_record)),
        exact_theta_record_(std::move(exact_theta_record)),
        geometry_record_(std::move(geometry_record)),
        retained_geometry_(std::move(retained_geometry)), work_(work),
        blocks_(std::move(blocks)), couplings_(std::move(couplings)),
        physical_equation_(std::move(physical_equation)) {
    if (physical_equation_ &&
        physical_equation_->owner_signature_identity != exact_identity_)
      throw std::invalid_argument(
          "composite physical q/C payload names a different full parent owner identity");
    if constexpr (std::is_same_v<Scalar, Rational> ||
                  std::is_same_v<Scalar, ComplexBall>)
      pseudo_target_cache_ =
          std::make_unique<PseudoTargetCache<Scalar>>(exact_identity_);
  }

  CompositeColumnSolveResult solve_column(
      const std::string& local_handle,
      const json::object& request,
      std::shared_ptr<CompositeSCCChartBase> equation_owner) override {
    if (!equation_owner || equation_owner.get() != this)
      throw std::invalid_argument(
          "completed SCC solve received a different physical equation owner");
    if constexpr (std::is_same_v<Scalar, ComplexBall>) {
      return solve_regular_acb_column(
          local_handle, request, std::move(equation_owner));
    } else if constexpr (!std::is_same_v<Scalar, Rational>) {
      throw std::invalid_argument(
          "native SCC column execution supports exact rational columns and regular Acb columns only");
    } else {
      const auto started = std::chrono::steady_clock::now();
      const bool regular_execution = regular_block_column_ready();
      const bool regular_singular_execution =
          regular_singular_jordan_column_ready();
      if (!regular_execution && !regular_singular_execution)
        throw std::invalid_argument(
            "retained SCC chart does not satisfy a native exact-rational SCC column capability");
      const auto checkpoint_identity = required_string(
          request, "checkpoint_identity");
      const auto& seed_request = as_object(
          request.at("seed"), "native SCC seed request");
      const auto seed_block = as_u32(
          seed_request.at("block"), "native SCC seed block");
      if (seed_block >= blocks_.size())
        throw std::invalid_argument("native SCC seed block is out of range");

      std::vector<std::uint8_t> reachable(blocks_.size(), 0);
      reachable[seed_block] = 1;
      for (const auto block : graph_.topological_order) {
        if (!reachable[block]) continue;
        for (const auto [source, target] : graph_.condensation_edges)
          if (source == block) reachable[target] = 1;
      }
      std::vector<std::uint32_t> expected_targets;
      for (const auto block : graph_.topological_order)
        if (block != seed_block && reachable[block])
          expected_targets.push_back(block);
      const auto& target_requests = as_array(
          request.at("targets"), "native SCC target requests");
      if (target_requests.size() != expected_targets.size())
        throw std::invalid_argument(
            "native SCC targets must cover every reachable descendant exactly");
      for (std::size_t index = 0; index < target_requests.size(); ++index) {
        const auto& target = as_object(
            target_requests[index], "native SCC target request");
        if (as_u32(target.at("block"), "native SCC target block") !=
            expected_targets[index])
          throw std::invalid_argument(
              "native SCC targets are not in deterministic topological order");
      }

      std::vector<std::optional<LocalSolution<Scalar>>> state(blocks_.size());
      json::array diagnostics;
      NativeLocalDiagnostics aggregate;
      aggregate.top_valid = kCompleteInfinity;
      std::vector<ValidatedScheduleEvidence> validated_schedules;
      std::vector<std::unique_ptr<PseudoCompensator<Rational>>>
          pseudo_compensators(blocks_.size());
      const auto pseudo_compensator = [&](std::uint32_t block)
          -> PseudoCompensator<Rational>& {
        if (!pseudo_target_cache_)
          throw std::logic_error(
              "exact Rational SCC composite has no CASE-P target cache");
        if (!blocks_[block].exact_jordan_indicial.has_value())
          throw std::logic_error(
              "exact Rational SCC block has no CASE-P Jordan certificate");
        if (!pseudo_compensators[block])
          pseudo_compensators[block] =
              std::make_unique<PseudoCompensator<Rational>>(
                  *blocks_[block].chart,
                  *blocks_[block].exact_jordan_indicial, block,
                  checkpoint_identity + ":block:" +
                      std::to_string(block),
                  *pseudo_target_cache_);
        return *pseudo_compensators[block];
      };

      const auto& seed_run = checked_column_run(
          seed_request, seed_block, true, nullptr,
          regular_singular_execution, &validated_schedules);
      const auto seed_local_component = seed_component_from_run(
          seed_run, seed_block, regular_singular_execution);
      const auto basis_index =
          blocks_[seed_block].vertices[seed_local_component];
      auto seed_native = regular_singular_execution
          ? pseudo_compensator(seed_block).solve(
                seed_run,
                as_object(seed_request.at("metadata"),
                          "native SCC seed metadata"),
                std::nullopt, work_.work_min)
          : blocks_[seed_block].chart->solve_native(
                seed_run, as_object(seed_request.at("metadata"),
                                    "native SCC seed metadata"));
      validate_block_result(seed_native, seed_block, true,
                            regular_singular_execution);
      blocks_[seed_block].chart->record_native_local_success(
          seed_native.diagnostics);
      accumulate_diagnostics(aggregate, seed_native.diagnostics);
      diagnostics.push_back(block_diagnostic(
          seed_block, "seed", {}, nullptr, seed_native));
      state[seed_block] = std::move(seed_native.solution);

      for (std::size_t target_index = 0;
           target_index < target_requests.size(); ++target_index) {
        const auto target_block = expected_targets[target_index];
        const auto& target_request = as_object(
            target_requests[target_index], "native SCC target request");
        std::vector<LocalSolution<Scalar>> incoming;
        std::vector<std::uint32_t> predecessors;
        for (const auto& coupling : couplings_) {
          if (coupling.target_block != target_block ||
              !state[coupling.source_block].has_value())
            continue;
          validate_column_coupling(
              coupling, regular_singular_execution);
          auto source_physical = block_gauge_transform(
              *state[coupling.source_block], coupling.source_block, true,
              checkpoint_identity + ":source-gauge:" +
                  std::to_string(coupling.source_block));
          auto contribution = apply_prepared_sparse_local_matrix(
              coupling.matrix, source_physical,
              checkpoint_identity + ":source:" +
                  std::to_string(coupling.source_block) + ":" +
                  std::to_string(target_block));
          if (!contribution.has_value())
            throw std::logic_error(
                "an exact nonzero SCC edge produced no structural source");
          predecessors.push_back(coupling.source_block);
          incoming.push_back(std::move(*contribution));
        }
        if (incoming.empty())
          throw std::invalid_argument(
              "reachable native SCC target has no available predecessor source");
        auto source = incoming.size() == 1
            ? std::move(incoming.front())
            : combine_local_solutions(
                  incoming, checkpoint_identity + ":combined-source:" +
                                std::to_string(target_block));
        source = target_recurrence_source(
            std::move(source), target_block,
            checkpoint_identity + ":target-vinv:" +
                std::to_string(target_block));
        // A signed-shift halo is a property of the complete target source.
        // Restricting predecessors or individual VInv entries separately
        // would reject exact below-frame pieces which cancel only after the
        // complete physical source is summed and transformed.
        source = restrict_local_epsilon_frame_strict_lower(
            source, work_.work_min, work_.work_complete_max,
            checkpoint_identity + ":source-frame:" +
                std::to_string(target_block));
        require_work_local(source, "combined coupling source");
        if (!regular_singular_execution) {
          const auto& target_run = checked_column_run(
              target_request, target_block, false, &source, false,
              &validated_schedules);
          require_source_tag_matches_run(source, target_run);
          auto source_data = local_solution_source_data(
              source, as_u32(target_run.at("nmax"), "target nmax"),
              as_u32(target_run.at("p"), "target log maximum"),
              blocks_[target_block].chart->frame_base(),
              blocks_[target_block].chart->frame_width());
          auto target_native =
              blocks_[target_block].chart->solve_native_with_source(
                  target_run,
                  as_object(target_request.at("metadata"),
                            "native SCC target metadata"),
                  std::move(source_data));
          validate_block_result(target_native, target_block, false, false);
          blocks_[target_block].chart->record_native_local_success(
              target_native.diagnostics);
          accumulate_diagnostics(aggregate, target_native.diagnostics);
          diagnostics.push_back(block_diagnostic(
              target_block, "particular", predecessors, &source,
              target_native));
          state[target_block] = std::move(target_native.solution);
          continue;
        }

        // CASE-P compensation introduces homogeneous target-root sectors.
        // They are linearly independent exact tags and must remain separate
        // through every later SCC edge.  Solve one exact particular per tag,
        // then recombine; collapsing them into a single recurrence would
        // silently apply the wrong affine schedule to all but one sector.
        auto groups = split_exact_rational_tags(
            source, checkpoint_identity + ":source-groups:" +
                        std::to_string(target_block));
        const auto& submitted_run = as_object(
            target_request.at("run"), "native SCC target run");
        const Rational submitted_a =
            parse_scalar<Rational>(submitted_run.at("a_target"));
        const Rational submitted_b =
            parse_scalar<Rational>(submitted_run.at("b_target"));
        bool submitted_used = false;
        std::vector<LocalSolution<Rational>> target_parts;
        target_parts.reserve(groups.size());
        for (std::size_t group_index = 0; group_index < groups.size();
             ++group_index) {
          auto& group = groups[group_index];
          if (group.sectors.empty())
            throw std::logic_error("exact SCC source tag group is empty");
          const Rational group_a(group.sectors.front().a.canonical);
          const Rational group_b(group.sectors.front().b.canonical);
          const bool use_submitted = group_a == submitted_a &&
                                     group_b == submitted_b;
          json::object derived_entry;
          const json::object* entry = &target_request;
          if (!use_submitted) {
            std::uint32_t source_log = 0;
            for (const auto& sector : group.sectors)
              source_log = std::max(source_log, sector.log_power);
            auto run = exact_derived_run(
                submitted_run, *blocks_[target_block].chart,
                *blocks_[target_block].exact_jordan_indicial,
                group_a, group_b, source_log, false, std::nullopt);
            auto metadata = exact_derived_metadata(
                as_object(target_request.at("metadata"),
                          "native SCC target metadata"),
                group_a, group_b,
                as_u32(run.at("p"), "derived particular log maximum"),
                ":derived-tag:" + std::to_string(group_index));
            derived_entry = json::object{
                {"block", target_block}, {"run", std::move(run)},
                {"metadata", std::move(metadata)}};
            entry = &derived_entry;
          } else {
            if (submitted_used)
              throw std::logic_error(
                  "exact SCC source contains duplicate submitted tag groups");
            submitted_used = true;
          }

          const auto& target_run = checked_column_run(
              *entry, target_block, false, &group, true,
              &validated_schedules);
          require_source_tag_matches_run(group, target_run);
          auto source_data = local_solution_source_data(
              group, as_u32(target_run.at("nmax"), "target nmax"),
              as_u32(target_run.at("p"), "target log maximum"),
              blocks_[target_block].chart->frame_base(),
              blocks_[target_block].chart->frame_width());
          auto target_native = pseudo_compensator(target_block).solve(
              target_run,
              as_object(entry->at("metadata"),
                        "native SCC target metadata"),
              std::move(source_data), work_.work_min);
          validate_block_result(target_native, target_block, false, true);
          blocks_[target_block].chart->record_native_local_success(
              target_native.diagnostics);
          accumulate_diagnostics(aggregate, target_native.diagnostics);
          auto diagnostic = block_diagnostic(
              target_block, "particular-tag", predecessors, &group,
              target_native);
          diagnostic["source_a"] = group_a.str();
          diagnostic["source_b"] = group_b.str();
          diagnostics.push_back(std::move(diagnostic));
          target_parts.push_back(std::move(target_native.solution));
        }
        // A polar cross-SCC multiplier may shift every incoming source away
        // from the submitted prototype tag.  In that case every particular
        // above used an exact derived run and the prototype was correctly
        // unused.  Requiring its tag to occur would reject a closed exact
        // source decomposition merely because the wire request supplies one
        // deterministic schedule template.
        auto target_state = target_parts.size() == 1
            ? std::move(target_parts.front())
            : combine_local_solutions(
                  target_parts,
                  checkpoint_identity + ":target-tag-sum:" +
                      std::to_string(target_block));
        require_work_local(target_state, "combined target tag particulars");
        state[target_block] = std::move(target_state);
      }

      std::vector<LocalSolution<Scalar>> embedded;
      for (std::uint32_t block = 0; block < state.size(); ++block) {
        if (!state[block].has_value()) continue;
        auto physical = block_gauge_transform(
            *state[block], block, true,
            checkpoint_identity + ":final-gauge:" + std::to_string(block));
        embedded.push_back(local_algebra_detail::embedded_components(
            physical, blocks_[block].vertices, dimension_));
      }
      if (embedded.empty())
        throw std::logic_error("native SCC column produced no block state");
      auto parent = embedded.size() == 1
          ? std::move(embedded.front())
          : combine_local_solutions(
                embedded, checkpoint_identity + ":work-parent");
      require_work_local(parent, "combined parent work state");
      auto retained_match_max = static_cast<std::int32_t>(
          std::min<std::int64_t>(
              parent.epsilon.complete_max,
              work_.work_complete_max));
      if (!aggregate.requires_parent_completeness_certificate &&
          aggregate.top_valid != kCompleteInfinity)
        retained_match_max =
            std::min(retained_match_max, aggregate.top_valid);
      if (retained_match_max < work_.requested_max)
        throw std::domain_error(
            "exact SCC column does not retain the requested public epsilon maximum");
      const auto parent_completeness = certify_parent_completeness(
          parent, aggregate, seed_block, seed_local_component, basis_index,
          retained_match_max, regular_singular_execution, reachable,
          expected_targets, state, validated_schedules);
      if (parent_completeness.has_value()) {
        const auto& certified_epsilon = as_object(
            parent_completeness->at("epsilon"),
            "SCC parent completeness epsilon");
        retained_match_max = as_i32(
            certified_epsilon.at("complete_max"),
            "SCC parent completeness maximum");
        aggregate.top_valid = retained_match_max;
      }
      parent = cap_composite_public_local(
          parent, retained_match_max, work_.public_t_order,
          retained_geometry_.chart, retained_geometry_.prescriptions,
          checkpoint_identity);
      validate_local_solution(parent, false);
      const auto scalar_execution = scalar_block_shape();
      json::object column_identity_record{
          {"schema", regular_singular_execution
               ? (scalar_execution
                     ? "diffexp2-native-scc-regular-singular-scalar-column-v1"
                     : "diffexp2-native-scc-regular-singular-jordan-column-v2")
               : (scalar_execution
                     ? "diffexp2-native-scc-column-v1"
                     : "diffexp2-native-scc-column-v2")},
          {"scc_exact_identity", exact_identity_},
          {"basis_index", basis_index},
          {"seed", seed_request},
          {"targets", target_requests}};
      if (!scalar_execution)
        column_identity_record["seed_local_component"] =
            seed_local_component;
      column_identity_record["pseudo_compensation"] =
          aggregate.pseudo_hits == 0
              ? "none"
              : "exact-rational-derived-jordan-targets-v1";
      if (parent_completeness.has_value())
        column_identity_record["parent_completeness_certificate"] =
            *parent_completeness;
      SCCColumnProvenance column_provenance{
          handle_, exact_identity_, seed_block,
          basis_index,
          json::serialize(canonical_json_value(column_identity_record))};
      auto local = retain_completed_parent_local(
          local_handle, seed_block, std::move(parent), aggregate,
          std::move(column_provenance), std::move(equation_owner));
      const auto elapsed_ms = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      column_solves_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(column_stats_mutex_);
        column_solve_ms_ += elapsed_ms;
      }
      return {std::move(local), std::move(diagnostics), elapsed_ms};
    }
  }

  double column_solve_ms() const {
    std::lock_guard<std::mutex> lock(column_stats_mutex_);
    return column_solve_ms_;
  }

  const char* column_execution_capability() const override {
    if (regular_singular_jordan_column_ready())
      return regular_singular_column_capability();
    if (!regular_block_column_ready())
      return "unsupported-native-scc-column";
    return regular_column_capability();
  }

  const std::string& geometry_record() const override {
    return geometry_record_;
  }
  std::uint32_t dimension() const override { return dimension_; }
  const CompositeWorkContract& work_contract() const override {
    return work_;
  }

  std::optional<std::pair<FiniteLaurentVector<ComplexBall>, std::string>>
  normalize_acb_matching_vector(
      const FiniteLaurentVector<ComplexBall>& physical) const override {
    if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
      (void)physical;
      return std::nullopt;
    } else {
      if (physical.size() != dimension_)
        throw std::invalid_argument(
            "SCC matching normal-frame vector has the wrong parent dimension");
      const bool epsilon_regular = std::any_of(
          blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.chart->uses_epsilon_regular_principal();
          });
      if (!epsilon_regular) return std::nullopt;

      std::vector<std::optional<EpsilonFrame<ComplexBall>>> rows(
          dimension_);
      std::string frame_identity =
          exact_identity_ + ":epsilon-regular-block-spectral-match-frame";
      for (const auto& block : blocks_) {
        // A t-dependent gauge needs an exact point specialization rather
        // than a finite Taylor replay.  Keep this first normal-frame scope
        // deliberately to identity-gauge blocks; the banana confluent chart
        // and its surrounding SCC satisfy that stronger certificate.
        if (!block.to_reduced.identity)
          return std::nullopt;
        const auto& transform = block.source_transform;
        if (!spectral_frame_ready(block))
          throw std::logic_error(
              "epsilon-regular SCC matching lost its certified spectral frame");
        frame_identity += ":" + transform.producer_identity;
        if (transform.identity) {
          for (std::size_t component = 0;
               component < block.vertices.size(); ++component) {
            const auto vertex = block.vertices[component];
            rows[vertex] = rows[vertex].has_value()
                ? *rows[vertex] + physical[vertex]
                : physical[vertex];
          }
          continue;
        }
        for (const auto& entry : transform.matrix.entries) {
          if (entry.row >= block.vertices.size() ||
              entry.column >= block.vertices.size())
            throw std::logic_error(
                "spectral matching transform entry is outside its block");
          const auto& multiplier = entry.multiplier;
          std::vector<ComplexBall> coefficients;
          coefficients.reserve(multiplier.kernels.size());
          for (const auto& kernel : multiplier.kernels) {
            if (kernel.empty())
              throw std::logic_error(
                  "spectral matching transform has an empty Taylor kernel");
            for (std::size_t taylor = 1; taylor < kernel.size(); ++taylor)
              if (!kernel[taylor].is_zero())
                throw std::logic_error(
                    "spectral matching transform unexpectedly depends on the chart variable");
            coefficients.push_back(kernel.front());
          }
          auto term = EpsilonFrame<ComplexBall>(
              multiplier.epsilon_shift, std::move(coefficients)) *
              physical[block.vertices[entry.column]];
          const auto output = block.vertices[entry.row];
          rows[output] = rows[output].has_value()
              ? *rows[output] + term
              : std::move(term);
        }
      }
      FiniteLaurentVector<ComplexBall> normalized;
      normalized.reserve(dimension_);
      for (std::uint32_t row = 0; row < dimension_; ++row) {
        if (!rows[row].has_value())
          throw std::logic_error(
              "spectral matching normal frame has an empty parent row");
        normalized.push_back(std::move(*rows[row]));
      }
      return std::pair{std::move(normalized), std::move(frame_identity)};
    }
  }

  std::optional<FiniteLaurentMatrix<ComplexBall>>
  right_normalize_acb_matching_basis(
      const FiniteLaurentMatrix<ComplexBall>& left_normalized)
      const override {
    if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
      (void)left_normalized;
      return std::nullopt;
    } else {
      if (left_normalized.size() != dimension_ ||
          std::any_of(left_normalized.begin(), left_normalized.end(),
              [&](const auto& row) { return row.size() != dimension_; }))
        throw std::invalid_argument(
            "SCC right matching normalization requires a square parent basis");
      if (!std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.chart->uses_epsilon_regular_principal();
          }))
        return std::nullopt;

      FiniteLaurentMatrix<ComplexBall> normalized(
          dimension_, FiniteLaurentVector<ComplexBall>());
      for (auto& row : normalized) row.reserve(dimension_);
      for (std::uint32_t parent_row = 0; parent_row < dimension_;
           ++parent_row) {
        std::vector<std::optional<EpsilonFrame<ComplexBall>>> columns(
            dimension_);
        for (const auto& block : blocks_) {
          if (block.source_transform.identity) {
            for (const auto vertex : block.vertices)
              columns[vertex] = left_normalized[parent_row][vertex];
            continue;
          }
          const auto frame_base = block.chart->frame_base();
          const auto frame_top = matching_detail::checked_power(
              static_cast<std::int64_t>(frame_base) +
                  block.chart->frame_width() - 1,
              "SCC spectral assembly matching maximum");
          if (frame_base > 0 || frame_top < 0)
            throw std::logic_error(
                "SCC spectral assembly frame does not contain epsilon^0");
          FiniteLaurentMatrix<ComplexBall> assembly(
              block.vertices.size(), FiniteLaurentVector<ComplexBall>());
          for (auto& row : assembly) row.reserve(block.vertices.size());
          for (std::size_t spectral_column = 0;
               spectral_column < block.vertices.size(); ++spectral_column) {
            FiniteLaurentVector<ComplexBall> unit;
            unit.reserve(block.vertices.size());
            for (std::size_t component = 0;
                 component < block.vertices.size(); ++component) {
              std::vector<ComplexBall> coefficients(
                  static_cast<std::size_t>(frame_top) + 1,
                  ScalarTraits<ComplexBall>::zero());
              if (component == spectral_column)
                coefficients.front() = ScalarTraits<ComplexBall>::one();
              unit.emplace_back(0, std::move(coefficients));
            }
            auto assembled =
                block.chart->apply_assembly_to_acb_matching_vector(unit);
            for (std::size_t row = 0; row < block.vertices.size(); ++row)
              assembly[row].push_back(std::move(assembled[row]));
          }
          for (std::size_t spectral_column = 0;
               spectral_column < block.vertices.size(); ++spectral_column) {
            std::optional<EpsilonFrame<ComplexBall>> value;
            for (std::size_t physical_column = 0;
                 physical_column < block.vertices.size(); ++physical_column) {
              auto term = left_normalized[parent_row]
                              [block.vertices[physical_column]] *
                  assembly[physical_column][spectral_column];
              value = value.has_value()
                  ? *value + term : std::move(term);
            }
            columns[block.vertices[spectral_column]] = std::move(*value);
          }
        }
        for (auto& column : columns) {
          if (!column.has_value())
            throw std::logic_error(
                "right spectral matching normalization left an empty column");
          normalized[parent_row].push_back(std::move(*column));
        }
      }
      return normalized;
    }
  }

  std::optional<FiniteLaurentVector<ComplexBall>>
  denormalize_acb_matching_weights(
      const FiniteLaurentVector<ComplexBall>& normalized) const override {
    if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
      (void)normalized;
      return std::nullopt;
    } else {
      if (normalized.size() != dimension_)
        throw std::invalid_argument(
            "SCC matching weight denormalization has the wrong dimension");
      if (!std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.chart->uses_epsilon_regular_principal();
          }))
        return std::nullopt;
      std::vector<std::optional<EpsilonFrame<ComplexBall>>> physical(
          dimension_);
      for (const auto& block : blocks_) {
        FiniteLaurentVector<ComplexBall> local;
        local.reserve(block.vertices.size());
        for (const auto vertex : block.vertices)
          local.push_back(normalized[vertex]);
        if (!block.source_transform.identity) {
          const auto frame_base = block.chart->frame_base();
          const auto frame_top = matching_detail::checked_power(
              static_cast<std::int64_t>(frame_base) +
                  block.chart->frame_width() - 1,
              "SCC matching weight assembly maximum");
          if (frame_base > 0 || frame_top < 0)
            throw std::logic_error(
                "SCC matching weight assembly does not contain epsilon^0");
          FiniteLaurentMatrix<ComplexBall> assembly(
              block.vertices.size(), FiniteLaurentVector<ComplexBall>());
          for (auto& row : assembly) row.reserve(block.vertices.size());
          for (std::size_t column = 0; column < block.vertices.size();
               ++column) {
            FiniteLaurentVector<ComplexBall> unit;
            unit.reserve(block.vertices.size());
            for (std::size_t component = 0;
                 component < block.vertices.size(); ++component) {
              std::vector<ComplexBall> coefficients(
                  static_cast<std::size_t>(frame_top) + 1,
                  ScalarTraits<ComplexBall>::zero());
              if (component == column)
                coefficients.front() = ScalarTraits<ComplexBall>::one();
              unit.emplace_back(0, std::move(coefficients));
            }
            auto assembled =
                block.chart->apply_assembly_to_acb_matching_vector(unit);
            for (std::size_t row = 0; row < block.vertices.size(); ++row)
              assembly[row].push_back(std::move(assembled[row]));
          }
          FiniteLaurentVector<ComplexBall> assembled_weights;
          assembled_weights.reserve(block.vertices.size());
          for (std::size_t row = 0; row < block.vertices.size(); ++row) {
            std::optional<EpsilonFrame<ComplexBall>> value;
            for (std::size_t column = 0; column < block.vertices.size();
                 ++column) {
              auto term = assembly[row][column] * local[column];
              value = value.has_value()
                  ? *value + term : std::move(term);
            }
            assembled_weights.push_back(std::move(*value));
          }
          local = std::move(assembled_weights);
        }
        for (std::size_t component = 0; component < block.vertices.size();
             ++component)
          physical[block.vertices[component]] = std::move(local[component]);
      }
      FiniteLaurentVector<ComplexBall> result;
      result.reserve(dimension_);
      for (auto& frame : physical) {
        if (!frame.has_value())
          throw std::logic_error(
              "SCC matching weight denormalization left an empty component");
        result.push_back(std::move(*frame));
      }
      return result;
    }
  }

  const char* equation_scalar_domain() const override {
    if constexpr (std::is_same_v<Scalar, Rational>) return "rational";
    if constexpr (std::is_same_v<Scalar, ComplexBall>) return "acb";
    return "symbolic";
  }

  std::shared_ptr<const void> physical_ode_erased() const override {
    return std::static_pointer_cast<const void>(physical_equation_);
  }

  const std::string& physical_payload_identity() const override {
    static const std::string empty;
    return physical_equation_ ? physical_equation_->payload_identity : empty;
  }

  const std::string& physical_payload_record() const override {
    static const std::string empty;
    return physical_equation_
        ? physical_equation_->exact_payload_record : empty;
  }

  const std::string& owner_signature_identity() const override {
    static const std::string empty;
    return physical_equation_
        ? physical_equation_->owner_signature_identity : empty;
  }

  std::vector<std::shared_ptr<PreparedChartBase>>
  dependency_charts() const override {
    std::vector<std::shared_ptr<PreparedChartBase>> result;
    result.reserve(blocks_.size());
    for (const auto& block : blocks_) result.push_back(block.chart);
    return result;
  }

  json::object stats_json() const override {
    std::size_t active_entries = 0, proven_zero_entries = 0;
    std::optional<std::int32_t> min_coupling_shift;
    std::optional<std::int32_t> max_coupling_shift;
    json::array block_handles;
    block_handles.reserve(blocks_.size());
    for (const auto& block : blocks_) {
      json::object block_record{
          {"block", block.block}, {"chart", block.source_handle},
          {"dimension", block.chart->dimension()},
          {"regular", block.regular},
          {"no_pseudo", block.no_pseudo},
          {"principal_identity", block.principal_identity}};
      if (const auto& indicial = block.exact_jordan_indicial;
          indicial.has_value()) {
        json::array indicial_blocks;
        std::uint32_t max_jordan_size = 0;
        for (const auto& spectral_block : indicial->blocks) {
          max_jordan_size = std::max(max_jordan_size,
                                     spectral_block.size());
          indicial_blocks.push_back(json::object{
              {"block", spectral_block.block_index},
              {"columns", encode_indices(spectral_block.columns)},
              {"jordan_size", spectral_block.size()},
              {"a", spectral_block.root.a.str()},
              {"b", spectral_block.root.b.str()}});
        }
        block_record["exact_affine_jordan_indicial"] = json::object{
            {"dimension", indicial->dimension},
            {"blocks", std::move(indicial_blocks)},
            {"max_jordan_size", max_jordan_size}};
        // Preserve the scalar-v1 diagnostic field while deriving it from the
        // same complete exact proof used by multidimensional admission.
        if (indicial->dimension == 1 && indicial->blocks.size() == 1)
          block_record["affine_indicial_root"] = json::object{
              {"a", indicial->blocks.front().root.a.str()},
              {"b", indicial->blocks.front().root.b.str()}};
      } else if (block.chart->exact_jordan_indicial_error().has_value()) {
        block_record["exact_affine_jordan_indicial_error"] =
            *block.chart->exact_jordan_indicial_error();
      }
      block_handles.push_back(std::move(block_record));
    }
    for (const auto& coupling : couplings_)
      for (std::size_t index = 0; index < coupling.identities.size(); ++index) {
        const auto& identity = coupling.identities[index];
        identity.proven_zero ? ++proven_zero_entries : ++active_entries;
        if (!identity.proven_zero) {
          const auto shift = coupling.matrix.entries[index]
                                 .multiplier.epsilon_shift;
          min_coupling_shift = min_coupling_shift.has_value()
              ? std::min(*min_coupling_shift, shift) : shift;
          max_coupling_shift = max_coupling_shift.has_value()
              ? std::max(*max_coupling_shift, shift) : shift;
        }
      }
    const auto regular_ready = regular_block_column_ready();
    const auto regular_singular_ready =
        regular_singular_jordan_column_ready();
    const auto execution_ready = regular_ready || regular_singular_ready;
    const auto scalar_shape = scalar_block_shape();
    const auto scalar_ready = scalar_column_ready();
    const auto pseudo_cache_stats = pseudo_target_cache_
        ? pseudo_target_cache_->stats()
        : typename PseudoTargetCache<Scalar>::Stats{};
    json::object result{
        {"scc", handle_}, {"key", key_}, {"identity", exact_identity_},
        {"rational_shadow_identity", rational_shadow_identity_},
        {"dimension", dimension_}, {"blocks", blocks_.size()},
        {"physical_ode_owner",
         physical_equation_ ? "full-parent" : "unsupported"},
        {"physical_payload_identity",
         physical_equation_
             ? json::value(physical_equation_->payload_identity)
             : json::value(nullptr)},
        {"coupling_groups", couplings_.size()},
        {"coupling_entries", active_entries + proven_zero_entries},
        {"active_coupling_entries", active_entries},
        {"proven_zero_coupling_entries", proven_zero_entries},
        {"frame_base", work_.work_min},
        {"cancellation_audit_base", work_.cancellation_audit_base.has_value()
             ? json::value(*work_.cancellation_audit_base)
             : json::value(nullptr)},
        {"frame_width", static_cast<std::uint32_t>(
             static_cast<std::int64_t>(work_.work_complete_max) -
             work_.work_min + 1)},
        {"requested_min", work_.requested_min},
        {"requested_max", work_.requested_max},
        {"work_complete_max", work_.work_complete_max},
        {"public_t_order", work_.public_t_order},
        {"work_t_order", work_.work_t_order},
        {"wolfram_coupling_depth", work_.wolfram_coupling_depth},
        {"native_coupling_depth", graph_.coupling_depth},
        {"min_coupling_shift", min_coupling_shift.has_value()
             ? json::value(*min_coupling_shift) : json::value(nullptr)},
        {"max_coupling_shift", max_coupling_shift.has_value()
             ? json::value(*max_coupling_shift) : json::value(nullptr)},
        {"execution_mode", "BlockSequentialStrict"},
        {"execution_implemented", execution_ready},
        {"execution_scope", regular_singular_ready
             ? regular_singular_column_capability()
             : (regular_ready
                   ? regular_column_capability()
                   : "unsupported")},
        {"general_scc_execution", false},
        {"scalar_block_dag_column_execution", scalar_ready},
        {"regular_singular_scalar_block_dag_column_execution",
         regular_singular_ready && scalar_shape},
        {"regular_singular_jordan_block_dag_column_execution",
         regular_singular_ready},
        {"scc_column_solves", column_solves_.load()},
        {"scc_column_solve_ms", column_solve_ms()},
        {"casep_homogeneous_target_cache_scope",
         pseudo_target_cache_
             ? "immutable-composite" : "not-applicable"},
        {"casep_homogeneous_target_cache_entries",
         pseudo_cache_stats.entries},
        {"casep_homogeneous_target_cache_builds",
         pseudo_cache_stats.builds},
        {"casep_homogeneous_target_cache_hits",
         pseudo_cache_stats.hits},
        {"capability_evidence", json::object{
             {"identity_v",
              "native-retained-spectral-assembly-and-target-inverse"},
             {"regular", regular_ready
                  ? "collision-bound-producer-certificate"
                  : "not-required-by-selected-scope"},
             {"regular_or_regular_singular",
              "collision-bound-producer-certificate"},
             {"identity_gauge",
              "optional-fast-path-exact-directional-gauge-frame"},
             {"exact_gauge",
              "retained-gauge-and-gauge-inverse-multiplier-certificates"},
             {"no_pseudo", regular_singular_ready
                  ? (std::is_same_v<Scalar, ComplexBall>
                        ? "runtime-exact-schedule-case-p-gate"
                        : "producer-provenance-only-execution-revalidated-by-exact-schedule-certificate")
                  : "collision-bound-producer-certificate"},
             {"jordan_indicial",
              regular_singular_ready
                  ? "retained-exact-rational-full-matrix-certificate"
                  : "not-required-by-selected-scope"},
             {"pseudo_schedule_execution",
              regular_singular_ready
                  ? (std::is_same_v<Scalar, ComplexBall>
                        ? "exact-rational-certificate-case-p-rejected-for-acb"
                        : "exact-rational-joint-compensation-and-formal-overlap-certificate")
                  : "not-required-by-selected-scope"},
             {"resonance_schedule",
              regular_singular_ready
                  ? "retained-affine-jordan-verified-exact-captured-run"
                  : "retained-affine-root-verified-exact-captured-run"}}},
        {"execution_must_revalidate_producer_capabilities", true},
        {"block_charts", std::move(block_handles)}};
    if (!scalar_shape)
      result["regular_block_dag_column_execution"] = regular_ready;
    return result;
  }

 private:
  struct ValidatedScheduleEvidence {
    std::uint32_t block = 0;
    std::uint32_t taylor_rows = 0;
    bool exact_affine_jordan = false;
    std::size_t validation_runs = 0;
  };

  json::object certify_column_uniqueness_path(
      std::uint32_t seed_block, std::uint32_t seed_local_component,
      std::uint32_t basis_index, bool regular_singular_execution,
      const std::vector<std::uint8_t>& submitted_reachable,
      const std::vector<std::uint32_t>& submitted_targets,
      const std::vector<std::optional<LocalSolution<Scalar>>>& state,
      const std::vector<ValidatedScheduleEvidence>& schedules) const {
    if (seed_block >= blocks_.size() ||
        seed_local_component >= blocks_[seed_block].vertices.size() ||
        blocks_[seed_block].vertices[seed_local_component] != basis_index)
      throw std::domain_error(
          "deferred SCC completeness lost its validated canonical seed column");
    if (submitted_reachable.size() != blocks_.size() ||
        state.size() != blocks_.size() ||
        graph_.component_count != blocks_.size() ||
        graph_.topological_order.size() != blocks_.size())
      throw std::domain_error(
          "deferred SCC completeness lost its retained exact block-DAG shape");
    const std::set<std::pair<std::uint32_t, std::uint32_t>> graph_edges(
        graph_.condensation_edges.begin(), graph_.condensation_edges.end());
    std::set<std::pair<std::uint32_t, std::uint32_t>> coupling_edges;
    for (const auto& coupling : couplings_)
      if (!coupling_edges.emplace(coupling.source_block,
                                  coupling.target_block).second)
        throw std::domain_error(
            "deferred SCC completeness has duplicate coupling groups for one block-DAG edge");
    if (coupling_edges != graph_edges)
      throw std::domain_error(
          "deferred SCC completeness coupling set contradicts its exact block DAG");

    std::vector<std::uint8_t> reachable(blocks_.size(), 0);
    reachable[seed_block] = 1;
    for (const auto block : graph_.topological_order) {
      if (block >= blocks_.size())
        throw std::domain_error(
            "deferred SCC completeness has an invalid topological block");
      if (!reachable[block]) continue;
      for (const auto [source, target] : graph_.condensation_edges) {
        if (source >= blocks_.size() || target >= blocks_.size())
          throw std::domain_error(
              "deferred SCC completeness has an invalid condensation edge");
        if (source == block) reachable[target] = 1;
      }
    }
    std::vector<std::uint32_t> expected_targets;
    std::size_t reachable_count = 0;
    for (const auto block : graph_.topological_order) {
      if (!reachable[block]) continue;
      ++reachable_count;
      if (block != seed_block) expected_targets.push_back(block);
    }
    if (reachable != submitted_reachable ||
        expected_targets != submitted_targets)
      throw std::domain_error(
          "deferred SCC completeness execution is not the retained deterministic block-DAG path");
    for (std::size_t block = 0; block < state.size(); ++block)
      if (state[block].has_value() != static_cast<bool>(reachable[block]))
        throw std::domain_error(
            "deferred SCC completeness parent does not cover exactly its reachable block DAG");

    std::vector<std::uint8_t> schedule_covered(blocks_.size(), 0);
    std::size_t validated_schedule_runs = 0;
    for (const auto& evidence : schedules) {
      if (evidence.block >= blocks_.size() ||
          !reachable[evidence.block] ||
          evidence.taylor_rows !=
              static_cast<std::size_t>(work_.work_t_order) + 1 ||
          evidence.exact_affine_jordan != regular_singular_execution ||
          evidence.validation_runs == 0 || schedule_covered[evidence.block])
        throw std::domain_error(
            "deferred SCC completeness has inconsistent validated recurrence-schedule evidence");
      schedule_covered[evidence.block] = 1;
      validated_schedule_runs += evidence.validation_runs;
    }
    if (schedules.size() != reachable_count)
      throw std::domain_error(
          "deferred SCC completeness schedule evidence does not cover each reachable block exactly once");
    for (std::size_t block = 0; block < blocks_.size(); ++block) {
      if (!reachable[block]) continue;
      if (!schedule_covered[block])
        throw std::domain_error(
            "deferred SCC completeness has a reachable block without a validated recurrence schedule");
      if (regular_singular_execution &&
          !blocks_[block].exact_jordan_indicial.has_value())
        throw std::domain_error(
            "deferred SCC completeness lost an exact affine-Jordan certificate");
      if (!regular_singular_execution &&
          !blocks_[block].chart->has_regular_singleton_partition())
        throw std::domain_error(
            "deferred SCC completeness lost an exact regular singleton schedule");
    }
    return json::object{
        {"schema", "diffexp2-scc-column-uniqueness-path-v1"},
        {"scc_exact_identity", exact_identity_},
        {"canonical_seed",
         json::object{{"block", seed_block},
                      {"local_component", seed_local_component},
                      {"basis_index", basis_index},
                      {"normalization",
                       regular_singular_execution
                           ? "validated-exact-affine-jordan-log-unit"
                           : "validated-exact-regular-epsilon-unit"}}},
        {"reachable_blocks", reachable_count},
        {"validated_schedule_blocks", schedules.size()},
        {"validated_schedule_runs", validated_schedule_runs},
        {"schedule_proof",
         regular_singular_execution
             ? "certify_exact_affine_jordan_schedule"
             : "exact-regular-resonant-zero-and-Taylor-index"},
        {"block_dag_proof",
         "validated-SCC-condensation-and-deterministic-reachable-cover"}};
  }

  std::optional<json::object> certify_parent_completeness(
      const LocalSolution<Scalar>& parent,
      const NativeLocalDiagnostics& diagnostics,
      std::uint32_t seed_block, std::uint32_t seed_local_component,
      std::uint32_t basis_index, std::int32_t claimed_complete_max,
      bool regular_singular_execution,
      const std::vector<std::uint8_t>& reachable,
      const std::vector<std::uint32_t>& expected_targets,
      const std::vector<std::optional<LocalSolution<Scalar>>>& state,
      const std::vector<ValidatedScheduleEvidence>& schedules) const {
    if (!diagnostics.requires_parent_completeness_certificate)
      return std::nullopt;
    if constexpr (!std::is_same_v<Scalar, Rational>) {
      throw std::domain_error(
          "deferred SCC completeness reached a non-Rational parent after its exact-domain gate");
    } else {
      if (!physical_equation_)
        throw std::domain_error(
            "deferred SCC completeness requires the retained full-parent physical q/C owner");
      if (parent.epsilon.min_power > claimed_complete_max)
        throw std::domain_error(
            "deferred SCC completeness has no nonempty claimed epsilon window");
      const EpsilonWindow claimed{
          parent.epsilon.min_power, claimed_complete_max};
      SCCFormalResidualCertificate residual;
      try {
        residual = certify_scc_parent_exact_formal_residual(
            *physical_equation_, parent, claimed, work_.public_t_order,
            true);
      } catch (const std::domain_error& error) {
        std::string state_summary;
        for (std::uint32_t block = 0; block < state.size(); ++block) {
          if (!state[block].has_value()) continue;
          auto physical = block_gauge_transform(
              *state[block], block, true,
              "failed-parent-residual:block:" + std::to_string(block));
          if (!state_summary.empty()) state_summary += ";";
          state_summary += "block=" + std::to_string(block) +
              ",reduced=" + exact_local_sector_summary(*state[block]) +
              ",physical=" + exact_local_sector_summary(physical);
        }
        throw std::domain_error(
            "SCC parent completeness failed for basis_index=" +
            std::to_string(basis_index) + ", seed_block=" +
            std::to_string(seed_block) + ", seed_local_component=" +
            std::to_string(seed_local_component) + ": " + error.what() +
            "; state=[" + state_summary + "]");
      }
      const auto uniqueness = certify_column_uniqueness_path(
          seed_block, seed_local_component, basis_index,
          regular_singular_execution, reachable, expected_targets, state,
          schedules);
      if (residual.epsilon.complete_max < work_.requested_max) {
        // Recover the exact first failing coefficient in the public contract;
        // the prefix probe above is used only to retain honest private
        // reservoir beyond that mandatory edge.
        (void)certify_scc_parent_exact_formal_residual(
            *physical_equation_, parent,
            {parent.epsilon.min_power, work_.requested_max},
            work_.public_t_order);
        throw std::logic_error(
            "SCC parent prefix probe returned below the requested maximum without a residual witness");
      }
      return json::object{
          {"schema", "diffexp2-scc-parent-exact-formal-completeness-v1"},
          {"owner_signature_identity",
           physical_equation_->owner_signature_identity},
          {"physical_payload_identity",
           physical_equation_->payload_identity},
          {"epsilon", json::object{{"min", residual.epsilon.min_power},
                                    {"complete_max",
                                     residual.epsilon.complete_max}}},
          {"taylor_complete_max", residual.taylor_complete_max},
          {"exact_tag_count", residual.exact_tag_count},
          {"coefficient_rows", residual.coefficient_rows},
          {"residual",
           "coefficientwise-exact-rational-q-theta-minus-C-equality"},
          {"reservoir",
           "every-active-q/C-product-covered-through-claimed-complete-max"},
          {"uniqueness", uniqueness},
          {"seed_block", seed_block},
          {"seed_local_component", seed_local_component},
          {"basis_index", basis_index},
          {"recurrence_top_valid", encode_validity(diagnostics.top_valid)}};
    }
  }

  std::shared_ptr<StoredLocalBase> retain_completed_parent_local(
      const std::string& local_handle, std::uint32_t seed_block,
      LocalSolution<Scalar> parent, const NativeLocalDiagnostics& diagnostics,
      SCCColumnProvenance column_provenance,
      std::shared_ptr<CompositeSCCChartBase> equation_owner) {
    std::shared_ptr<PhysicalEquationOwnerBase> retained_owner;
    std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
        retained_equation;
    std::string unavailable_reason;
    if constexpr (std::is_same_v<Scalar, Rational> ||
                  std::is_same_v<Scalar, ComplexBall>) {
      if (physical_equation_) {
        retained_owner = std::move(equation_owner);
        retained_equation = physical_equation_;
      } else {
        unavailable_reason =
            "owner-bound residual is unsupported: retained composite SCC has no full parent physical q/C payload";
      }
    }
    return make_retained_typed_shared<Scalar, StoredLocal<Scalar>>(
        local_handle, handle_, exact_identity_, std::move(parent),
        blocks_.at(seed_block).chart->precision_bits(),
        std::vector<PseudoHit<Scalar>>{}, diagnostics,
        std::move(column_provenance), std::nullopt, nullptr,
        unavailable_tail_model(
            "tail model is unavailable for an assembled SCC parent local"),
        std::nullopt, true, true, std::move(retained_owner),
        std::move(retained_equation), std::move(unavailable_reason));
  }

  CompositeColumnSolveResult solve_regular_acb_column(
      const std::string& local_handle, const json::object& request,
      std::shared_ptr<CompositeSCCChartBase> equation_owner) {
    static_assert(std::is_same_v<Scalar, ComplexBall>);
    const bool regular_singular_execution =
        regular_singular_jordan_column_ready();
    if (!regular_block_column_ready() && !regular_singular_execution)
      throw std::invalid_argument(
          "retained Acb SCC chart does not satisfy a strict regular or exact affine-Jordan no-CASE-P block-DAG column capability");

    // Coupling products and block recombination happen between the individual
    // PreparedChart solves, so the complete composite operation must retain
    // the same Arb precision lease rather than relying on each nested solve's
    // shorter lease.
    AcbPrecisionLease acb_lease(blocks_.front().chart->precision_bits());
    ComplexBall::set_precision(blocks_.front().chart->precision_bits());
    const auto started = std::chrono::steady_clock::now();
    const auto checkpoint_identity = required_string(
        request, "checkpoint_identity");
    const auto& seed_request = as_object(
        request.at("seed"), "native Acb SCC seed request");
    const auto seed_block = as_u32(
        seed_request.at("block"), "native Acb SCC seed block");
    if (seed_block >= blocks_.size())
      throw std::invalid_argument("native Acb SCC seed block is out of range");

    std::vector<std::uint8_t> reachable(blocks_.size(), 0);
    reachable[seed_block] = 1;
    for (const auto block : graph_.topological_order) {
      if (!reachable[block]) continue;
      for (const auto [source, target] : graph_.condensation_edges)
        if (source == block) reachable[target] = 1;
    }
    std::vector<std::uint32_t> expected_targets;
    for (const auto block : graph_.topological_order)
      if (block != seed_block && reachable[block])
        expected_targets.push_back(block);
    const auto& target_requests = as_array(
        request.at("targets"), "native Acb SCC target requests");
    if (target_requests.size() != expected_targets.size())
      throw std::invalid_argument(
          "native Acb SCC targets must cover every reachable descendant exactly");
    for (std::size_t index = 0; index < target_requests.size(); ++index) {
      const auto& target = as_object(
          target_requests[index], "native Acb SCC target request");
      if (as_u32(target.at("block"), "native Acb SCC target block") !=
          expected_targets[index])
        throw std::invalid_argument(
            "native Acb SCC targets are not in deterministic topological order");
    }

    std::vector<std::optional<LocalSolution<Scalar>>> state(blocks_.size());
    json::array diagnostics;
    NativeLocalDiagnostics aggregate;
    aggregate.top_valid = kCompleteInfinity;
    std::vector<std::unique_ptr<PseudoCompensator<ComplexBall>>>
        pseudo_compensators(blocks_.size());
    const auto pseudo_compensator = [&](std::uint32_t block)
        -> PseudoCompensator<ComplexBall>& {
      if (!pseudo_target_cache_)
        throw std::logic_error(
            "Acb SCC composite has no CASE-P target cache");
      if (!blocks_[block].exact_jordan_indicial.has_value())
        throw std::logic_error(
            "Acb SCC block has no CASE-P Jordan certificate");
      if (!pseudo_compensators[block])
        pseudo_compensators[block] =
            std::make_unique<PseudoCompensator<ComplexBall>>(
                *blocks_[block].chart,
                *blocks_[block].exact_jordan_indicial, block,
                checkpoint_identity + ":block:" +
                    std::to_string(block),
                *pseudo_target_cache_);
      return *pseudo_compensators[block];
    };

    const auto& seed_run = checked_column_run(
        seed_request, seed_block, true, nullptr,
        regular_singular_execution);
    const auto seed_local_component = seed_component_from_run(
        seed_run, seed_block, regular_singular_execution);
    const auto basis_index =
        blocks_[seed_block].vertices[seed_local_component];
    auto seed_native = regular_singular_execution
        ? pseudo_compensator(seed_block).solve(
              seed_run,
              as_object(seed_request.at("metadata"),
                        "native Acb SCC seed metadata"),
              std::nullopt, work_.work_min)
        : blocks_[seed_block].chart->solve_native(
              seed_run, as_object(seed_request.at("metadata"),
                                  "native Acb SCC seed metadata"));
    validate_block_result(seed_native, seed_block, true,
                          regular_singular_execution);
    blocks_[seed_block].chart->record_native_local_success(
        seed_native.diagnostics);
    accumulate_diagnostics(aggregate, seed_native.diagnostics);
    diagnostics.push_back(block_diagnostic(
        seed_block, "seed", {}, nullptr, seed_native));
    state[seed_block] = std::move(seed_native.solution);

    for (std::size_t target_index = 0;
         target_index < target_requests.size(); ++target_index) {
      const auto target_block = expected_targets[target_index];
      const auto& target_request = as_object(
          target_requests[target_index], "native Acb SCC target request");
      std::vector<LocalSolution<Scalar>> incoming;
      std::vector<std::uint32_t> predecessors;
      for (const auto& coupling : couplings_) {
        if (coupling.target_block != target_block ||
            !state[coupling.source_block].has_value())
          continue;
        validate_column_coupling(coupling, regular_singular_execution);
        auto source_physical = block_gauge_transform(
            *state[coupling.source_block], coupling.source_block, true,
            checkpoint_identity + ":source-gauge:" +
                std::to_string(coupling.source_block));
        auto contribution = apply_prepared_sparse_local_matrix(
            coupling.matrix, source_physical,
            checkpoint_identity + ":source:" +
                std::to_string(coupling.source_block) + ":" +
                std::to_string(target_block));
        if (!contribution.has_value())
          throw std::logic_error(
              "a certified nonzero Acb SCC edge produced no structural source");
        predecessors.push_back(coupling.source_block);
        incoming.push_back(std::move(*contribution));
      }
      if (incoming.empty())
        throw std::invalid_argument(
            "reachable native Acb SCC target has no available predecessor source");
      auto source = incoming.size() == 1
          ? std::move(incoming.front())
          : combine_local_solutions(
                incoming, checkpoint_identity + ":combined-source:" +
                              std::to_string(target_block));
      source = target_recurrence_source(
          std::move(source), target_block,
          checkpoint_identity + ":target-vinv:" +
              std::to_string(target_block));
      source = restrict_local_epsilon_frame_strict_lower(
          source, work_.work_min, work_.work_complete_max,
          checkpoint_identity + ":source-frame:" +
              std::to_string(target_block));
      require_work_local(source, "combined Acb coupling source");

      if (!regular_singular_execution) {
        const auto& target_run = checked_column_run(
            target_request, target_block, false, &source, false);
        require_source_tag_matches_run(source, target_run);
        auto source_data = local_solution_source_data(
            source, as_u32(target_run.at("nmax"), "target nmax"),
            as_u32(target_run.at("p"), "target log maximum"),
            blocks_[target_block].chart->frame_base(),
            blocks_[target_block].chart->frame_width());
        auto target_native =
            blocks_[target_block].chart->solve_native_with_source(
                target_run,
                as_object(target_request.at("metadata"),
                          "native Acb SCC target metadata"),
                std::move(source_data));
        validate_block_result(target_native, target_block, false, false);
        blocks_[target_block].chart->record_native_local_success(
            target_native.diagnostics);
        accumulate_diagnostics(aggregate, target_native.diagnostics);
        diagnostics.push_back(block_diagnostic(
            target_block, "particular", predecessors, &source,
            target_native));
        state[target_block] = std::move(target_native.solution);
        continue;
      }

      // Exact tags determine separate affine schedules even though the
      // coefficient field is Acb.  This is the same decomposition used by
      // the Rational fallback, with cancellation accepted only when the
      // resulting enclosure is below the configured chop floor.
      auto groups = split_exact_rational_tags(
          source, checkpoint_identity + ":source-groups:" +
                      std::to_string(target_block));
      const auto& submitted_run = as_object(
          target_request.at("run"), "native Acb SCC target run");
      const Rational submitted_a =
          parse_scalar<Rational>(submitted_run.at("a_target"));
      const Rational submitted_b =
          parse_scalar<Rational>(submitted_run.at("b_target"));
      bool submitted_used = false;
      std::vector<LocalSolution<ComplexBall>> target_parts;
      target_parts.reserve(groups.size());
      for (std::size_t group_index = 0; group_index < groups.size();
           ++group_index) {
        auto& group = groups[group_index];
        if (group.sectors.empty())
          throw std::logic_error("Acb SCC source tag group is empty");
        const Rational group_a(group.sectors.front().a.canonical);
        const Rational group_b(group.sectors.front().b.canonical);
        const bool use_submitted = group_a == submitted_a &&
                                   group_b == submitted_b;
        json::object derived_entry;
        const json::object* entry = &target_request;
        if (!use_submitted) {
          std::uint32_t source_log = 0;
          for (const auto& sector : group.sectors)
            source_log = std::max(source_log, sector.log_power);
          auto run = exact_derived_run(
              submitted_run, *blocks_[target_block].chart,
              *blocks_[target_block].exact_jordan_indicial,
              group_a, group_b, source_log, false, std::nullopt);
          auto metadata = exact_derived_metadata(
              as_object(target_request.at("metadata"),
                        "native Acb SCC target metadata"),
              group_a, group_b,
              as_u32(run.at("p"), "derived particular log maximum"),
              ":derived-tag:" + std::to_string(group_index));
          derived_entry = json::object{
              {"block", target_block}, {"run", std::move(run)},
              {"metadata", std::move(metadata)}};
          entry = &derived_entry;
        } else {
          if (submitted_used)
            throw std::logic_error(
                "Acb SCC source contains duplicate submitted tag groups");
          submitted_used = true;
        }

        const auto& target_run = checked_column_run(
            *entry, target_block, false, &group, true);
        require_source_tag_matches_run(group, target_run);
        auto source_data = local_solution_source_data(
            group, as_u32(target_run.at("nmax"), "target nmax"),
            as_u32(target_run.at("p"), "target log maximum"),
            blocks_[target_block].chart->frame_base(),
            blocks_[target_block].chart->frame_width());
        auto target_native = pseudo_compensator(target_block).solve(
            target_run,
            as_object(entry->at("metadata"),
                      "native Acb SCC target metadata"),
            std::move(source_data), work_.work_min);
        validate_block_result(target_native, target_block, false, true);
        blocks_[target_block].chart->record_native_local_success(
            target_native.diagnostics);
        accumulate_diagnostics(aggregate, target_native.diagnostics);
        auto diagnostic = block_diagnostic(
            target_block, "particular-tag", predecessors, &group,
            target_native);
        diagnostic["source_a"] = group_a.str();
        diagnostic["source_b"] = group_b.str();
        diagnostics.push_back(std::move(diagnostic));
        target_parts.push_back(std::move(target_native.solution));
      }
      auto target_state = target_parts.size() == 1
          ? std::move(target_parts.front())
          : combine_local_solutions(
                target_parts,
                checkpoint_identity + ":target-tag-sum:" +
                    std::to_string(target_block));
      require_work_local(target_state, "combined Acb target tag particulars");
      state[target_block] = std::move(target_state);
    }

    std::vector<LocalSolution<Scalar>> embedded;
    for (std::uint32_t block = 0; block < state.size(); ++block) {
      if (!state[block].has_value()) continue;
      auto physical = block_gauge_transform(
          *state[block], block, true,
          checkpoint_identity + ":final-gauge:" + std::to_string(block));
      embedded.push_back(local_algebra_detail::embedded_components(
          physical, blocks_[block].vertices, dimension_));
    }
    if (embedded.empty())
      throw std::logic_error("native Acb SCC column produced no block state");
    auto parent = embedded.size() == 1
        ? std::move(embedded.front())
        : combine_local_solutions(
              embedded, checkpoint_identity + ":work-parent");
    require_work_local(parent, "combined Acb parent work state");
    // A retained SCC column is a receiving-basis object, not a public
    // observable.  Preserve every coefficient that the recurrence actually
    // proved complete.  The public requirement remains requested_max and is
    // checked below; an arbitrary dimension-scaled cap would discard private
    // Laurent-match reservoir without reducing any upstream computation.
    const auto retained_match_max = static_cast<std::int32_t>(
        std::min<std::int64_t>(
            parent.epsilon.complete_max,
            work_.work_complete_max));
    if (aggregate.requires_parent_completeness_certificate)
      throw std::logic_error(
          "native Acb SCC column reached publication with deferred completeness instead of requesting its exact Rational shadow");
    parent = cap_composite_public_local(
        parent, retained_match_max, work_.public_t_order,
        retained_geometry_.chart, retained_geometry_.prescriptions,
        checkpoint_identity);
    validate_local_solution(parent, false);

    const auto scalar_execution = scalar_block_shape();
    json::object column_identity_record{
        {"schema", regular_singular_execution
             ? (scalar_execution
                   ? "diffexp2-native-scc-acb-regular-singular-scalar-column-v1"
                   : "diffexp2-native-scc-acb-regular-singular-jordan-column-v1")
             : (scalar_execution
                   ? "diffexp2-native-scc-acb-regular-scalar-column-v1"
                   : "diffexp2-native-scc-acb-regular-column-v2")},
        {"scc_exact_identity", exact_identity_},
        {"basis_index", basis_index},
        {"seed", seed_request},
        {"targets", target_requests},
        {"pseudo_compensation",
         aggregate.pseudo_hits == 0
             ? "none"
             : "exact-schedule-acb-ball-floor-certified-v1"}};
    if (!scalar_execution)
      column_identity_record["seed_local_component"] =
          seed_local_component;
    SCCColumnProvenance column_provenance{
        handle_, exact_identity_, seed_block, basis_index,
        json::serialize(canonical_json_value(column_identity_record))};
    auto local = retain_completed_parent_local(
        local_handle, seed_block, std::move(parent), aggregate,
        std::move(column_provenance), std::move(equation_owner));
    const auto elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    column_solves_.fetch_add(1);
    {
      std::lock_guard<std::mutex> lock(column_stats_mutex_);
      column_solve_ms_ += elapsed_ms;
    }
    return {std::move(local), std::move(diagnostics), elapsed_ms};
  }

  bool regular_block_column_ready() const {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      return false;
    } else {
      if (blocks_.size() < 2 || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return !block.regular || block.vertices.empty() ||
                   block.chart->dimension() != block.vertices.size() ||
                   !spectral_frame_ready(block) || !gauge_frame_ready(block) ||
                   !block.chart->has_regular_singleton_partition();
          }))
        return false;
      return std::all_of(
          couplings_.begin(), couplings_.end(), [&](const auto& coupling) {
            return regular_coupling_ready(coupling);
          });
    }
  }

  bool scalar_column_ready() const {
    return regular_block_column_ready() && scalar_block_shape();
  }

  const char* regular_column_capability() const {
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return scalar_block_shape()
          ? "acb-regular-scalar-block-dag-column-v1"
          : "acb-regular-block-dag-column-v2";
    else
      return scalar_block_shape()
          ? "exact-rational-regular-scalar-block-dag-column-v1"
          : "exact-rational-regular-block-dag-column-v2";
  }

  bool regular_singular_jordan_column_ready() const {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      return false;
    } else {
      if (blocks_.empty() || work_.work_min > 0 ||
          work_.requested_max < 0 || work_.work_complete_max < 0 ||
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return block.vertices.empty() ||
                   block.chart->dimension() != block.vertices.size() ||
                   !spectral_frame_ready(block) || !gauge_frame_ready(block) ||
                   !block.exact_jordan_indicial.has_value();
          }))
        return false;
      // Singularity is a property of the complete triangular system, not
      // only of its diagonal SCCs.  A polar (or center-nonvanishing theta)
      // cross coupling can make the parent regular-singular while every
      // diagonal block remains ordinary.  Such couplings are already closed
      // by the exact-tagged source algebra below; require at least one
      // genuinely non-ordinary block or coupling so the ordinary path keeps
      // its narrower capability and provenance.
      const auto singular_parent =
          std::any_of(blocks_.begin(), blocks_.end(), [](const auto& block) {
            return !block.regular;
          }) ||
          std::any_of(couplings_.begin(), couplings_.end(),
                      [&](const auto& coupling) {
                        return !regular_coupling_ready(coupling);
                      });
      if (!singular_parent) return false;
      // Acb admission is schedule-specific, not producer-global.  A block
      // may advertise possible family collisions even when the submitted
      // finite Taylor schedule contains no CASE-P step.  checked_column_run
      // reconstructs every affine offset exactly and rejects an actual
      // CASE-P before any Acb recurrence executes.
      return std::all_of(
          couplings_.begin(), couplings_.end(), [&](const auto& coupling) {
            return exact_tagged_coupling_ready(coupling);
          });
    }
  }

  const char* regular_singular_column_capability() const {
    const auto scalar = scalar_block_shape();
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return scalar
          ? "acb-regular-singular-scalar-block-dag-column-v1"
          : "acb-regular-singular-jordan-block-dag-column-v1";
    else
      return scalar
          ? "exact-rational-regular-singular-scalar-block-dag-column-v1"
          : "exact-rational-regular-singular-jordan-block-dag-column-v2";
  }

  bool scalar_block_shape() const {
    return std::all_of(
        blocks_.begin(), blocks_.end(), [](const auto& block) {
          return block.vertices.size() == 1;
        });
  }

  static bool spectral_frame_ready(
      const CompositeSCCBlock<Scalar>& block) {
    const auto& transform = block.source_transform;
    if (!block.chart->has_assembly() ||
        !transform.epsilon_unimodular ||
        transform.matrix.rows != block.chart->dimension() ||
        transform.matrix.columns != block.chart->dimension() ||
        transform.identity != block.chart->has_identity_assembly())
      return false;
    if (transform.identity)
      return true;  // Backward-compatible identity manifests need no payload.
    if (transform.producer_identity.empty() ||
        transform.v_exact_identity.empty() ||
        transform.vinv_exact_identity.empty() ||
        transform.det_exact_identity.empty() ||
        transform.v_exact_identity !=
            block.chart->assembly_exact_identity() ||
        transform.matrix.entries.empty())
      return false;
    return true;
  }

  static bool gauge_frame_ready(const CompositeSCCBlock<Scalar>& block) {
    const auto ready = [&](const auto& transform, const char* role) {
      if (transform.identity)
        return transform.role == role &&
            transform.matrix.rows == block.chart->dimension() &&
            transform.matrix.columns == block.chart->dimension();
      return transform.role == role && !transform.producer_identity.empty() &&
          !transform.gauge_exact_identity.empty() &&
          !transform.gauge_inverse_exact_identity.empty() &&
          !transform.gauge_det_exact_identity.empty() &&
          transform.matrix.rows == block.chart->dimension() &&
          transform.matrix.columns == block.chart->dimension() &&
          !transform.matrix.entries.empty();
    };
    return ready(block.to_physical, "to_physical") &&
        ready(block.to_reduced, "to_reduced") &&
        block.to_physical.identity == block.to_reduced.identity &&
        (block.to_physical.identity ||
         (block.to_physical.gauge_exact_identity ==
              block.to_reduced.gauge_exact_identity &&
          block.to_physical.gauge_inverse_exact_identity ==
              block.to_reduced.gauge_inverse_exact_identity &&
          block.to_physical.gauge_det_exact_identity ==
              block.to_reduced.gauge_det_exact_identity));
  }

  bool regular_coupling_ready(
      const CompositeSCCCoupling<Scalar>& coupling) const {
    if (!exact_tagged_coupling_ready(coupling)) return false;
    for (const auto& entry : coupling.matrix.entries) {
      if (entry.multiplier.proven_zero) continue;
      if (entry.multiplier.center_pole_order != 0 ||
          std::any_of(entry.multiplier.kernels.begin(),
                      entry.multiplier.kernels.end(),
                      [](const auto& kernel) {
                        return !ScalarTraits<Scalar>::is_zero(kernel.front());
                      }))
        return false;
    }
    return true;
  }

  // A prepared center pole is closed in the exact local algebra: t^-p
  // changes the exact sector tag from a to a-p before source grouping.
  // Rational regular-singular execution derives one schedule per resulting
  // tag; Acb either keeps the submitted tag or requests its Rational shadow.
  // Ordinary regular execution retains its stricter pole-free check above.
  bool exact_tagged_coupling_ready(
      const CompositeSCCCoupling<Scalar>& coupling) const {
    if (coupling.source_block >= blocks_.size() ||
        coupling.target_block >= blocks_.size() ||
        coupling.matrix.columns !=
            blocks_[coupling.source_block].vertices.size() ||
        coupling.matrix.rows !=
            blocks_[coupling.target_block].vertices.size())
      return false;
    bool active = false;
    for (const auto& entry : coupling.matrix.entries) {
      if (entry.row >= coupling.matrix.rows ||
          entry.column >= coupling.matrix.columns)
        return false;
      if (entry.multiplier.proven_zero) continue;
      active = true;
      if (entry.multiplier.kernels.empty() ||
          std::any_of(entry.multiplier.kernels.begin(),
                      entry.multiplier.kernels.end(),
                      [](const auto& kernel) {
                        return kernel.empty();
                      }))
        return false;
    }
    return active;
  }

  static bool scalar_identical(const Scalar& left, const Scalar& right) {
    if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return acb_equal(left.raw(), right.raw());
    else
      return left == right;
  }

  static Scalar scalar_from_exact_rational(const Rational& value) {
    if constexpr (std::is_same_v<Scalar, Rational>)
      return value;
    else if constexpr (std::is_same_v<Scalar, ComplexBall>)
      return ComplexBall::from_strings(value.str());
    else
      return SymbolicRational(value.str());
  }

  const json::object& checked_column_run(
      const json::object& entry, std::uint32_t block_index, bool seed,
      const LocalSolution<Scalar>* source,
      bool regular_singular_execution,
      std::vector<ValidatedScheduleEvidence>* validated_schedules =
          nullptr) const {
    if (block_index >= blocks_.size())
      throw std::invalid_argument("native SCC run block is out of range");
    const auto block_dimension = blocks_[block_index].chart->dimension();
    const auto frame_width = blocks_[block_index].chart->frame_width();
    const ExactJordanIndicialCertificate* retained_indicial = nullptr;
    if (regular_singular_execution) {
      const auto& certificate =
          blocks_[block_index].exact_jordan_indicial;
      if (!certificate.has_value())
        throw std::invalid_argument(
            "native regular-singular SCC chart has no retained exact affine Jordan indicial certificate");
      retained_indicial = &*certificate;
    }
    const auto& run = as_object(entry.at("run"), "native SCC recurrence run");
    const auto& metadata_object = as_object(
        entry.at("metadata"), "native SCC local metadata");
    validate_metadata_geometry(metadata_object);
    const auto* raw_audit = run.if_contains("cancellation_audit_base");
    const std::optional<std::int32_t> run_audit =
        raw_audit == nullptr || raw_audit->is_null()
            ? std::nullopt
            : std::optional<std::int32_t>(
                  as_i32(*raw_audit, "native SCC cancellation audit base"));
    if (run_audit != work_.cancellation_audit_base)
      throw std::invalid_argument(
          "native SCC run cancellation audit differs from its retained work contract");
    auto metadata = parse_local_metadata(metadata_object);
    if (!run.at("source").is_null())
      throw std::invalid_argument(
          "native SCC column rejects caller-supplied recurrence source data");
    if (run.at("return_u").as_bool())
      throw std::invalid_argument(
          "native SCC column requires retained assembly without U JSON");
    const auto nmax = as_u32(run.at("nmax"), "native SCC Taylor order");
    const auto log_max = as_u32(run.at("p"), "native SCC log maximum");
    if (nmax != work_.work_t_order)
      throw std::invalid_argument(
          "native SCC run Taylor order differs from its retained work contract");
    if (as_i32(run.at("a_shift_min"),
               "native SCC a-shift minimum") != 0)
      throw std::invalid_argument(
          "native SCC column requires a zero exact a-shift origin");
    const auto& a_shifts = as_array(
        run.at("a_shifts"), "native SCC exact a-shift schedule");
    if (a_shifts.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC regular block run has an incomplete a-shift schedule");
    std::optional<Rational> exact_a_target;
    std::optional<Rational> exact_b_target;
    if (regular_singular_execution) {
      if (metadata.a.domain != ExactDomain::Rational ||
          metadata.b.domain != ExactDomain::Rational)
        throw std::invalid_argument(
            "native regular-singular SCC execution requires exact Rational task tags");
      exact_a_target = Rational(metadata.a.canonical);
      exact_b_target = Rational(metadata.b.canonical);
      const auto encoded_a = parse_scalar<Scalar>(run.at("a_target"));
      const auto encoded_b = parse_scalar<Scalar>(run.at("b_target"));
      if (!scalar_identical(
              encoded_a, scalar_from_exact_rational(*exact_a_target)) ||
          !scalar_identical(
              encoded_b, scalar_from_exact_rational(*exact_b_target)))
        throw std::invalid_argument(
            "native regular-singular SCC numeric targets differ from their exact task tags");
      const auto& raw_tag = as_object(
          metadata_object.at("tag"), "native SCC exact task tag");
      if (const auto* raw_p = raw_tag.if_contains("p")) {
        const auto& p = as_object(*raw_p, "native SCC exact log tag");
        if (required_string(p, "domain") != "integer" ||
            Rational(required_string(p, "canonical")) !=
                Rational(std::to_string(log_max)))
          throw std::invalid_argument(
              "native regular-singular SCC exact log tag differs from its recurrence run");
      } else if constexpr (std::is_same_v<Scalar, ComplexBall>) {
        throw std::invalid_argument(
            "native Acb regular-singular SCC metadata must retain the exact log tag");
      }
      for (std::size_t n = 0; n < a_shifts.size(); ++n) {
        const auto expected =
            *exact_a_target + Rational(std::to_string(n));
        if (!scalar_identical(parse_scalar<Scalar>(a_shifts[n]),
                              scalar_from_exact_rational(expected)))
          throw std::invalid_argument(
              "native SCC a-shift schedule must equal a_target plus the exact Taylor index");
      }
    } else {
      const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
      const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
      if (log_max != 0 || !ScalarTraits<Scalar>::is_zero(a_target) ||
          !ScalarTraits<Scalar>::is_zero(b_target))
        throw std::invalid_argument(
            "native SCC regular block runs require p=0 and exact a=b=0");
      for (std::size_t n = 0; n < a_shifts.size(); ++n) {
        const auto shift = parse_scalar<Scalar>(a_shifts[n]);
        const auto expected = ScalarTraits<Scalar>::integer(
            static_cast<long>(n));
        if (!scalar_identical(shift, expected))
          throw std::invalid_argument(
              "native SCC regular a-shift schedule must equal the exact Taylor index");
      }
    }
    if (seed) {
      if (!run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC seed requires one initialized log-zero sector");
    } else {
      if (run.at("has_initial").as_bool())
        throw std::invalid_argument(
            "native SCC target rejects caller-supplied initial state");
      const auto& initial = as_array(
          run.at("initial"), "native SCC target initial tensor");
      const auto& validity = as_array(
          run.at("initial_validity"),
          "native SCC target initial validity");
      const auto log_points = static_cast<std::size_t>(log_max) + 1;
      if (log_points > std::numeric_limits<std::size_t>::max() /
                           block_dimension)
        throw std::overflow_error(
            "native SCC target initial validity size overflow");
      const auto expected_validity = log_points * block_dimension;
      if (expected_validity > std::numeric_limits<std::size_t>::max() /
                                  frame_width)
        throw std::overflow_error(
            "native SCC target initial tensor size overflow");
      const auto expected_initial_coefficients = expected_validity * frame_width;
      if (initial.size() != expected_initial_coefficients ||
          validity.size() != expected_validity ||
          std::any_of(initial.begin(), initial.end(), [](const auto& value) {
            return !ScalarTraits<Scalar>::is_zero(
                parse_scalar<Scalar>(value));
          }) ||
          std::any_of(validity.begin(), validity.end(), [](const auto& value) {
            return !value.is_null();
          }))
        throw std::invalid_argument(
            "native SCC target initial template must be explicit exact zero with unknown validity");
      if (source == nullptr)
        throw std::logic_error("native SCC target source was not constructed");
      if (source->dimension != block_dimension)
        throw std::invalid_argument(
            "native SCC coupling source dimension differs from its target block");
      for (const auto& sector : source->sectors)
        if (sector.log_power > log_max)
          throw std::invalid_argument(
              "native SCC target log template does not cover its exact source sectors");
    }
    const auto& schedule = as_array(
        run.at("schedule"), "native SCC recurrence schedule");
    if (schedule.size() != static_cast<std::size_t>(nmax) + 1)
      throw std::invalid_argument(
          "native SCC run has an inconsistent recurrence schedule height");
    if (regular_singular_execution) {
      std::vector<std::vector<BlockStep<Rational>>> exact_schedule;
      exact_schedule.reserve(schedule.size());
      for (std::size_t n = 0; n < schedule.size(); ++n) {
        const auto& raw_row = schedule[n];
        const auto& row = as_array(raw_row, "native SCC schedule row");
        if (row.size() != retained_indicial->blocks.size())
          throw std::invalid_argument(
              "native regular-singular SCC execution requires one step per retained exact Jordan block");
        std::vector<BlockStep<Rational>> exact_row;
        exact_row.reserve(row.size());
        for (std::size_t block = 0; block < row.size(); ++block) {
          const auto& raw_step = row[block];
          const auto& step = as_object(raw_step, "native SCC schedule step");
          const auto kind = required_string(step, "case");
          const auto step_case = kind == "T" ? StepCase::Taylor
              : kind == "P" ? StepCase::Pseudo
              : kind == "R" ? StepCase::Resonant
              : throw std::invalid_argument(
                    "native SCC schedule has an unknown recurrence case");
          if constexpr (std::is_same_v<Scalar, Rational>) {
            exact_row.push_back({step_case,
                parse_scalar<Rational>(step.at("da")),
                parse_scalar<Rational>(step.at("db"))});
          } else {
            const auto& exact_block = retained_indicial->blocks[block];
            const auto expected_da = *exact_a_target +
                Rational(std::to_string(n)) - exact_block.root.a;
            const auto expected_db =
                *exact_b_target - exact_block.root.b;
            if (!scalar_identical(
                    parse_scalar<Scalar>(step.at("da")),
                    scalar_from_exact_rational(expected_da)) ||
                !scalar_identical(
                    parse_scalar<Scalar>(step.at("db")),
                    scalar_from_exact_rational(expected_db)))
              throw std::invalid_argument(
                  "native Acb SCC schedule enclosure differs from its exact affine-Jordan offsets");
            exact_row.push_back(
                {step_case, expected_da, expected_db});
          }
        }
        exact_schedule.push_back(std::move(exact_row));
      }
      (void)certify_exact_affine_jordan_schedule(
          *retained_indicial, *exact_a_target, *exact_b_target,
          exact_schedule);
    } else {
      if (!blocks_[block_index].chart->has_regular_singleton_partition())
        throw std::invalid_argument(
            "native SCC regular block lost its retained singleton partition");
      for (std::size_t n = 0; n < schedule.size(); ++n) {
        const auto& row = as_array(schedule[n], "native SCC schedule row");
        if (row.size() != static_cast<std::size_t>(block_dimension))
          throw std::invalid_argument(
              "native SCC regular execution requires one step per retained Jordan singleton");
        const auto expected_kind = n == 0 ? "R" : "T";
        const auto expected_da = ScalarTraits<Scalar>::integer(
            static_cast<long>(n));
        for (const auto& raw_step : row) {
          const auto& step = as_object(raw_step, "native SCC schedule step");
          const auto kind = required_string(step, "case");
          const auto da = parse_scalar<Scalar>(step.at("da"));
          const auto db = parse_scalar<Scalar>(step.at("db"));
          if (kind != expected_kind ||
              !ScalarTraits<Scalar>::is_zero(db) ||
              !scalar_identical(da, expected_da))
            throw std::invalid_argument(
                "native SCC regular block schedule must be resonant at zero and Taylor by exact index for every Jordan singleton");
        }
      }
    }
    if (validated_schedules != nullptr) {
      const auto found = std::find_if(
          validated_schedules->begin(), validated_schedules->end(),
          [&](const auto& evidence) { return evidence.block == block_index; });
      if (found == validated_schedules->end()) {
        validated_schedules->push_back(ValidatedScheduleEvidence{
            block_index, static_cast<std::uint32_t>(schedule.size()),
            regular_singular_execution, 1});
      } else {
        if (found->taylor_rows != schedule.size() ||
            found->exact_affine_jordan != regular_singular_execution)
          throw std::logic_error(
              "validated SCC tag schedules disagree within one reachable block");
        ++found->validation_runs;
      }
    }
    return run;
  }

  std::uint32_t seed_component_from_run(
      const json::object& run, std::uint32_t block_index,
      bool regular_singular_execution) const {
    const auto dimension = blocks_[block_index].chart->dimension();
    const auto frame_width = blocks_[block_index].chart->frame_width();
    const auto& initial = as_array(
        run.at("initial"), "native SCC seed initial tensor");
    const auto& validity = as_array(
        run.at("initial_validity"), "native SCC seed initial validity");
    const auto unit_index_i64 = -static_cast<std::int64_t>(work_.work_min);
    if (regular_singular_execution) {
      const auto& indicial =
          blocks_[block_index].exact_jordan_indicial;
      if (!indicial.has_value())
        throw std::logic_error(
            "native regular-singular seed lost its retained indicial certificate");
      const auto log_max = as_u32(
          run.at("p"), "native regular-singular seed log maximum");
      const auto log_points = static_cast<std::size_t>(log_max) + 1;
      if (unit_index_i64 < 0 || unit_index_i64 >= frame_width ||
          log_points > std::numeric_limits<std::size_t>::max() / dimension ||
          log_points * dimension >
              std::numeric_limits<std::size_t>::max() / frame_width ||
          initial.size() != log_points * dimension * frame_width ||
          validity.size() != log_points * dimension)
        throw std::invalid_argument(
            "native regular-singular SCC seed has a malformed finite Jordan/log tensor");
      for (const auto& value : validity)
        if (value.is_null() ||
            as_i32(value, "native regular-singular seed validity") !=
                work_.work_complete_max)
          throw std::invalid_argument(
              "native regular-singular SCC seed requires finite validity through the retained work maximum");

      const auto unit_index = static_cast<std::size_t>(unit_index_i64);
      const auto initial_index = [&](std::uint32_t log,
                                     std::uint32_t component,
                                     std::size_t epsilon) {
        return ((static_cast<std::size_t>(log) * dimension + component) *
                frame_width) + epsilon;
      };
      const auto one = ScalarTraits<Scalar>::one();
      std::optional<std::uint32_t> selected;
      for (std::uint32_t component = 0; component < dimension; ++component) {
        for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
          const auto coefficient = parse_scalar<Scalar>(
              initial[initial_index(0, component, epsilon)]);
          if (epsilon == unit_index && scalar_identical(coefficient, one)) {
            if (selected.has_value())
              throw std::invalid_argument(
                  "native regular-singular SCC seed contains more than one log-zero eps^0 unit component");
            selected = component;
          } else if (!ScalarTraits<Scalar>::is_zero(coefficient)) {
            throw std::invalid_argument(
                "native regular-singular SCC seed log-zero frame is not one exact eps^0 unit component");
          }
        }
      }
      if (!selected.has_value())
        throw std::invalid_argument(
            "native regular-singular SCC seed contains no log-zero eps^0 unit component");

      const auto spectral_block_index = indicial->block_of_column[*selected];
      const auto position = indicial->position_in_block[*selected];
      const auto& spectral_block =
          indicial->blocks[spectral_block_index];
      const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
      const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
      if (!scalar_identical(
              a_target, scalar_from_exact_rational(spectral_block.root.a)) ||
          !scalar_identical(
              b_target, scalar_from_exact_rational(spectral_block.root.b)))
        throw std::invalid_argument(
            "native regular-singular SCC seed tag is not the exact affine root of its selected Jordan chain");
      if (log_max < position)
        throw std::invalid_argument(
            "native regular-singular SCC seed log ceiling truncates its exact Jordan chain");

      for (std::uint32_t log = 0; log <= log_max; ++log) {
        std::optional<std::uint32_t> expected_component;
        std::optional<std::size_t> expected_epsilon;
        if (log <= position) {
          expected_component = spectral_block.columns[position - log];
          const auto epsilon_i64 = -static_cast<std::int64_t>(log) -
                                   work_.work_min;
          if (epsilon_i64 < 0 || epsilon_i64 >= frame_width)
            throw RecurrenceError(
                "E4",
                "canonical Jordan seed exceeds the retained lower epsilon frame",
                work_.work_min, -static_cast<std::int32_t>(log));
          expected_epsilon = static_cast<std::size_t>(epsilon_i64);
        }
        for (std::uint32_t component = 0; component < dimension;
             ++component) {
          for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
            const auto coefficient = parse_scalar<Scalar>(
                initial[initial_index(log, component, epsilon)]);
            const auto expected = expected_component.has_value() &&
                                  component == *expected_component &&
                                  epsilon == *expected_epsilon;
            if ((expected && !scalar_identical(coefficient, one)) ||
                (!expected && !ScalarTraits<Scalar>::is_zero(coefficient)))
              throw std::invalid_argument(
                  "native regular-singular SCC seed differs from the captured canonical Jordan/log normalization");
          }
        }
      }
      return *selected;
    }
    if (unit_index_i64 < 0 || unit_index_i64 >= frame_width ||
        dimension > std::numeric_limits<std::size_t>::max() / frame_width ||
        initial.size() != static_cast<std::size_t>(dimension) * frame_width ||
        validity.size() != dimension)
      throw std::invalid_argument(
          "native SCC regular block seed requires one honest finite eps^0 unit frame");
    for (const auto& value : validity)
      if (value.is_null() ||
          as_i32(value, "native SCC seed validity") !=
              work_.work_complete_max)
        throw std::invalid_argument(
            "native SCC regular block seed requires finite validity through the retained work maximum");

    const auto unit_index = static_cast<std::size_t>(unit_index_i64);
    const auto one = ScalarTraits<Scalar>::one();
    std::optional<std::uint32_t> selected;
    for (std::uint32_t component = 0; component < dimension; ++component) {
      for (std::size_t epsilon = 0; epsilon < frame_width; ++epsilon) {
        const auto coefficient = parse_scalar<Scalar>(
            initial[static_cast<std::size_t>(component) * frame_width +
                    epsilon]);
        if (epsilon == unit_index && scalar_identical(coefficient, one)) {
          if (selected.has_value())
            throw std::invalid_argument(
                "native SCC seed contains more than one eps^0 unit component");
          selected = component;
        } else if (!ScalarTraits<Scalar>::is_zero(coefficient)) {
          throw std::invalid_argument(
              "native SCC regular block seed must be exactly one eps^0 unit column");
        }
      }
    }
    if (!selected.has_value())
      throw std::invalid_argument(
          "native SCC regular block seed contains no eps^0 unit component");
    return *selected;
  }

  void validate_metadata_geometry(const json::object& metadata) const {
    const auto& chart = as_object(
        metadata.at("chart"), "native SCC local chart metadata");
    const auto infinite = chart.at("infinite_radius").as_bool();
    if (required_string(chart, "center_exact") !=
            retained_geometry_.chart.center_exact ||
        required_string(chart, "scale_exact") !=
            retained_geometry_.chart.scale_exact ||
        infinite != retained_geometry_.chart.infinite_radius)
      throw std::invalid_argument(
          "native SCC local metadata differs from retained parent chart geometry");
    if (!infinite) {
      const auto* radius = chart.if_contains("radius");
      const auto* radius_exact = chart.if_contains("radius_exact");
      const bool exact_identity_matches =
          radius_exact != nullptr && radius_exact->is_string() &&
          radius_exact->as_string() == retained_geometry_.radius_exact;
      // Schema-v1 local metadata carried only a Rational radius string.
      // Preserve that path while allowing the exact algebraic identity plus
      // Acb specialization emitted by current Wolfram producers.
      const bool legacy_rational_matches =
          radius_exact == nullptr && radius != nullptr &&
          radius->is_string() &&
          Rational(std::string(radius->as_string())).str() ==
              retained_geometry_.radius_exact;
      if (radius == nullptr ||
          !(exact_identity_matches || legacy_rational_matches))
        throw std::invalid_argument(
            "native SCC local radius differs from retained exact parent radius");
    }

    const auto& prescriptions = as_array(
        metadata.at("prescriptions"),
        "native SCC local prescriptions");
    if (prescriptions.size() != retained_geometry_.prescriptions.size())
      throw std::invalid_argument(
          "native SCC local prescriptions differ from retained parent prescriptions");
    for (std::size_t index = 0; index < prescriptions.size(); ++index) {
      const auto& raw = as_object(
          prescriptions[index], "native SCC local prescription");
      const auto& retained = retained_geometry_.prescriptions[index];
      if (required_string(raw, "factor_exact") != retained.factor_exact ||
          as_i32(raw.at("sign"), "prescription sign") != retained.sign ||
          as_u32(raw.at("multiplicity"), "prescription multiplicity") !=
              retained.multiplicity ||
          as_i32(raw.at("leading_coefficient_sign"),
                 "prescription leading coefficient sign") !=
              retained.leading_coefficient_sign)
        throw std::invalid_argument(
            "native SCC local prescriptions differ from retained parent prescriptions");
    }
  }

  void require_work_local(const LocalSolution<Scalar>& solution,
                          const char* label) const {
    validate_local_solution(solution, false);
    if (!solution.error.empty() ||
        !local_algebra_detail::same_chart(
            solution.chart, retained_geometry_.chart) ||
        !local_algebra_detail::same_prescriptions(
            solution.prescriptions, retained_geometry_.prescriptions))
      throw std::invalid_argument(
          std::string(label) +
          " differs from retained composite geometry or carries an unsupported error envelope");
    if (solution.epsilon.min_power < work_.work_min ||
        solution.epsilon.complete_max > work_.work_complete_max ||
        solution.taylor_complete_max != work_.work_t_order)
      throw std::invalid_argument(
          std::string(label) + " lies outside the retained SCC work rectangle");
  }

  LocalSolution<Scalar> target_recurrence_source(
      LocalSolution<Scalar> physical_source, std::uint32_t target_block,
      std::string checkpoint_identity) const {
    if (target_block >= blocks_.size())
      throw std::invalid_argument(
          "native SCC target spectral transform block is out of range");
    const auto& block = blocks_[target_block];
    if (!spectral_frame_ready(block))
      throw std::invalid_argument(
          "native SCC target lacks a certified epsilon-unimodular spectral frame");
    if (physical_source.dimension != block.vertices.size())
      throw std::invalid_argument(
          "physical SCC source dimension differs from its target spectral frame");
    physical_source = block_gauge_transform(
        physical_source, target_block, false,
        checkpoint_identity + ":gauge-inverse");
    // The epsilon-regular recurrence accepts reduced physical source data
    // and performs V^-1 only inside its bounded resonant-layer transaction.
    // Retain the prepared V^-1 nevertheless: matching uses it as the exact
    // local normal frame for both sides of the handoff.
    if (block.chart->uses_epsilon_regular_principal())
      return physical_source;
    if (block.source_transform.identity)
      return physical_source;
    auto transformed = apply_prepared_sparse_local_matrix(
        block.source_transform.matrix, physical_source,
        std::move(checkpoint_identity));
    if (!transformed.has_value())
      throw std::logic_error(
          "certified nonidentity SCC VInv produced no structural source");
    return std::move(*transformed);
  }

  LocalSolution<Scalar> block_gauge_transform(
      const LocalSolution<Scalar>& source, std::uint32_t block_index,
      bool to_physical, std::string checkpoint_identity) const {
    if (block_index >= blocks_.size() ||
        source.dimension != blocks_[block_index].vertices.size() ||
        !gauge_frame_ready(blocks_[block_index]))
      throw std::invalid_argument(
          "native SCC block lacks a certified exact gauge frame");
    const auto& transform = to_physical
        ? blocks_[block_index].to_physical
        : blocks_[block_index].to_reduced;
    if (transform.identity) return source;
    auto result = apply_prepared_sparse_local_matrix(
        transform.matrix, source, std::move(checkpoint_identity));
    if (!result.has_value())
      throw std::logic_error(
          "certified nonidentity SCC gauge produced no structural state");
    return restrict_local_epsilon_frame_strict_lower(
        *result, work_.work_min, work_.work_complete_max,
        transform.role + ":work-frame");
  }

  void validate_block_result(const NativeLocalRun<Scalar>& native,
                             std::uint32_t block_index, bool seed,
                             bool regular_singular_execution) const {
    if (!native.pseudo_hits.empty())
      throw std::invalid_argument(
          "native SCC execution encountered uncompensated pseudo hits after its exact schedule certificate");
    if (block_index >= blocks_.size() ||
        native.solution.dimension != blocks_[block_index].vertices.size())
      throw std::invalid_argument(
          "native SCC block solve returned the wrong retained dimension");
    require_work_local(native.solution, "native SCC block result");
    if (native.solution.sectors.empty())
      throw std::invalid_argument("native SCC block solve returned no sectors");
    if (seed && !regular_singular_execution &&
        (native.solution.sectors.size() != 1 ||
         native.solution.sectors.front().log_power != 0))
      throw std::invalid_argument(
          "native SCC seed must assemble exactly one log-zero sector");
    for (const auto& sector : native.solution.sectors)
      if (sector.a.domain != ExactDomain::Rational ||
          sector.b.domain != ExactDomain::Rational)
        throw std::invalid_argument(
            "native SCC column execution requires exact rational sector tags");
  }

  void validate_column_coupling(
      const CompositeSCCCoupling<Scalar>& coupling,
      bool regular_singular_execution) const {
    const bool ready = regular_singular_execution
        ? exact_tagged_coupling_ready(coupling)
        : regular_coupling_ready(coupling);
    if (!ready)
      throw std::invalid_argument(
          regular_singular_execution
              ? "native regular-singular SCC column requires a dimension-matched exact coupling matrix closed under exact sector-tag shifts"
              : "native SCC column requires a dimension-matched exact pole-free coupling matrix whose active entries vanish at chart center");
  }

  void require_source_tag_matches_run(
      const LocalSolution<Scalar>& source, const json::object& run) const {
    const auto a_target = parse_scalar<Scalar>(run.at("a_target"));
    const auto b_target = parse_scalar<Scalar>(run.at("b_target"));
    const auto log_max = as_u32(run.at("p"), "target log maximum");
    for (const auto& sector : source.sectors) {
      if (sector.a.domain != ExactDomain::Rational ||
          sector.b.domain != ExactDomain::Rational ||
          sector.log_power > log_max)
        throw std::invalid_argument(
            "native SCC coupling source tag/log sector differs from its target run");
      verify_tag_binding<Scalar>(sector.a, a_target,
                                 "native SCC coupling source a");
      verify_tag_binding<Scalar>(sector.b, b_target,
                                 "native SCC coupling source b");
    }
  }

  static void accumulate_diagnostics(
      NativeLocalDiagnostics& total,
      const NativeLocalDiagnostics& current) {
    total.top_valid = std::min(total.top_valid, current.top_valid);
    total.parse_ms += current.parse_ms;
    total.kernel_ms += current.kernel_ms;
    total.requires_parent_completeness_certificate =
        total.requires_parent_completeness_certificate ||
        current.requires_parent_completeness_certificate;
  }

  static json::object block_diagnostic(
      std::uint32_t block, const char* role,
      const std::vector<std::uint32_t>& predecessors,
      const LocalSolution<Scalar>* source,
      const NativeLocalRun<Scalar>& result) {
    json::object diagnostic{
        {"block", block}, {"role", role},
        {"predecessors", encode_indices(predecessors)},
        {"result_epsilon_min", result.solution.epsilon.min_power},
        {"result_epsilon_max", result.solution.epsilon.complete_max},
        {"result_taylor_max", result.solution.taylor_complete_max},
        {"result_sectors", result.solution.sectors.size()},
        {"pseudo_hit_count", result.diagnostics.pseudo_hits},
        {"pseudo_compensation_count",
         result.diagnostics.pseudo_compensations},
        {"max_pseudo_depth", result.diagnostics.max_pseudo_depth},
        {"pseudo_value_certified",
         result.diagnostics.pseudo_value_certified},
        {"parent_completeness_certificate_required",
         result.diagnostics.requires_parent_completeness_certificate},
        {"uncompensated_pseudo_hit_count", result.pseudo_hits.size()},
        {"top_valid", encode_validity(result.diagnostics.top_valid)},
        {"parse_ms", result.diagnostics.parse_ms},
        {"kernel_ms", result.diagnostics.kernel_ms}};
    if (source != nullptr) {
      diagnostic["source_epsilon_min"] = source->epsilon.min_power;
      diagnostic["source_epsilon_max"] = source->epsilon.complete_max;
      diagnostic["source_taylor_max"] = source->taylor_complete_max;
      diagnostic["source_sectors"] = source->sectors.size();
    }
    return diagnostic;
  }

  std::uint32_t dimension_ = 0;
  SCCCertificate graph_;
  std::string exact_system_record_;
  std::string exact_theta_record_;
  std::string geometry_record_;
  RetainedCompositeGeometry retained_geometry_;
  CompositeWorkContract work_;
  std::vector<CompositeSCCBlock<Scalar>> blocks_;
  std::vector<CompositeSCCCoupling<Scalar>> couplings_;
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
      physical_equation_;
  std::unique_ptr<PseudoTargetCache<Scalar>> pseudo_target_cache_;
  std::atomic<std::uint64_t> column_solves_{0};
  mutable std::mutex column_stats_mutex_;
  double column_solve_ms_ = 0.0;
};
