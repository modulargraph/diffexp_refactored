// Certify the accuracy of a value-handoff coefficient without reducing an
// Arb enclosure to midpoint-relative bits.  The component radii and Acb
// upper magnitude are exact dyadic mag values; a zero-crossing coefficient
// therefore still passes when its absolute enclosure is sufficiently tight.
bool value_handoff_accurate(const ComplexBall& value,
                            const Rational& relative_error_max) {
  if (!value.is_finite() || relative_error_max.sign() <= 0 ||
      !(relative_error_max < Rational(1)))
    return false;
  mag_t scale;
  mag_init(scale);
  acb_get_mag(scale, value.raw());
  if (mag_is_inf(scale)) {
    mag_clear(scale);
    return false;
  }
  if (mag_cmp_2exp_si(scale, 0) < 0) mag_one(scale);

  fmpq_t threshold, scale_exact, allowed, radius_exact;
  fmpq_init(threshold);
  fmpq_init(scale_exact);
  fmpq_init(allowed);
  fmpq_init(radius_exact);
  const auto threshold_string = relative_error_max.str();
  const bool parsed =
      fmpq_set_str(threshold, threshold_string.c_str(), 10) == 0;
  if (parsed) {
    fmpq_canonicalise(threshold);
    mag_get_fmpq(scale_exact, scale);
    fmpq_mul(allowed, threshold, scale_exact);
  }
  const auto radius_passes = [&](const mag_t radius) {
    if (!parsed || mag_is_inf(radius)) return false;
    mag_get_fmpq(radius_exact, radius);
    return fmpq_cmp(radius_exact, allowed) <= 0;
  };
  const bool accurate =
      radius_passes(arb_radref(acb_realref(value.raw()))) &&
      radius_passes(arb_radref(acb_imagref(value.raw())));
  fmpq_clear(radius_exact);
  fmpq_clear(allowed);
  fmpq_clear(scale_exact);
  fmpq_clear(threshold);
  mag_clear(scale);
  return accurate;
}

bool is_retained_plan_value_handoff_schema(std::string_view schema) {
  return schema == "diffexp2-retained-plan-value-handoff-v1" ||
         schema == "diffexp2-retained-plan-value-handoff-v2";
}

std::string public_provenance_fingerprint(std::string_view identity) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto byte : identity) {
    hash ^= static_cast<std::uint8_t>(byte);
    hash *= UINT64_C(1099511628211);
  }
  constexpr char digits[] = "0123456789abcdef";
  std::string result = "fnv1a64:";
  for (int shift = 60; shift >= 0; shift -= 4)
    result.push_back(digits[(hash >> shift) & UINT64_C(0xf)]);
  return result;
}

std::vector<std::string> parse_symbols(const json::object& object) {
  std::vector<std::string> symbols;
  if (const auto* raw_symbols = object.if_contains("symbols")) {
    for (const auto& value : as_array(*raw_symbols, "symbolic variables")) {
      if (!value.is_string())
        throw std::invalid_argument("symbolic variable names must be strings");
      symbols.emplace_back(value.as_string());
    }
  }
  return symbols;
}

template <typename Scalar>
ExactEpsilonRational<Scalar> parse_physical_epsilon_rational(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  if (object.if_contains("zero") == nullptr ||
      !object.at("zero").is_bool())
    throw std::invalid_argument(std::string(label) +
                                " zero marker must be Boolean");
  ExactEpsilonRational<Scalar> result;
  result.zero = object.at("zero").as_bool();
  if (result.zero) {
    require_exact_keys(object, {"zero"}, label);
  } else {
    require_exact_keys(
        object, {"zero", "valuation", "numerator", "denominator"}, label);
    result.valuation = as_i32(object.at("valuation"),
                              "physical epsilon valuation");
    for (const auto& coefficient : as_array(
             object.at("numerator"), "physical epsilon numerator"))
      result.numerator.push_back(parse_scalar<Scalar>(coefficient));
    for (const auto& coefficient : as_array(
             object.at("denominator"), "physical epsilon denominator"))
      result.denominator.push_back(parse_scalar<Scalar>(coefficient));
  }
  physical_ode_detail::validate_rational(result, label);
  return result;
}

template <typename Scalar>
std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
parse_prepared_physical_ode(const json::value& raw,
                            std::uint32_t expected_dimension,
                            bool require_operator_token = true) {
  const auto& object = as_object(raw, "prepared physical cleared ODE");
  require_exact_keys(
      object,
      {"schema", "basis", "theta_coordinate", "owner_signature_identity",
       "payload_identity", "q", "c"},
      "prepared physical cleared ODE");
  if (required_string(object, "schema") !=
          "diffexp2-physical-cleared-ode-v1" ||
      required_string(object, "basis") != "physical-original-master" ||
      required_string(object, "theta_coordinate") != "local-t")
    throw std::invalid_argument(
        "prepared physical cleared ODE changes its equation basis or coordinate");
  auto result = std::make_shared<PreparedPhysicalClearedODE<Scalar>>();
  result->dimension = expected_dimension;
  result->owner_signature_identity = required_string(
      object, "owner_signature_identity");
  result->payload_identity = required_string(object, "payload_identity");
  if ((require_operator_token &&
       !result->owner_signature_identity.starts_with("de2-operator-")) ||
      !result->payload_identity.starts_with("de2-physical-ode-"))
    throw std::invalid_argument(
        "prepared physical cleared ODE has malformed owner/payload identity tokens");
  result->exact_payload_record = json::serialize(
      canonical_json_value(object));
  for (const auto& coefficient : as_array(object.at("q"), "physical q lags"))
    result->q_lags.push_back(parse_physical_epsilon_rational<Scalar>(
        coefficient, "physical q lag"));
  for (const auto& raw_lag : as_array(object.at("c"), "physical C lags")) {
    std::vector<PhysicalODEMatrixEntry<Scalar>> lag;
    for (const auto& raw_entry : as_array(raw_lag, "physical C lag")) {
      const auto& entry = as_object(raw_entry, "physical C entry");
      require_exact_keys(entry, {"r", "c", "v"}, "physical C entry");
      lag.push_back({as_u32(entry.at("r"), "physical C row"),
                     as_u32(entry.at("c"), "physical C column"),
                     parse_physical_epsilon_rational<Scalar>(
                         entry.at("v"), "physical C value")});
    }
    result->c_lags.push_back(std::move(lag));
  }
  physical_ode_detail::validate_ode(*result);
  return result;
}

std::string canonical_chart_geometry_record(const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  const auto center = required_string(geometry, "center_exact");
  const auto scale = required_string(geometry, "scale_exact");
  if (scale == "0")
    throw std::invalid_argument(
        "exact chart geometry requires a nonzero scale");
  const auto infinite = geometry.at("infinite_radius").as_bool();

  json::array prescriptions;
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "exact chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "exact chart prescription");
    const auto factor = required_string(prescription, "factor_exact");
    const auto sign = as_i32(prescription.at("sign"), "prescription sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"), "prescription multiplicity");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "prescription leading coefficient sign");
    if ((sign != -1 && sign != 1) || multiplicity == 0 ||
        (leading != -1 && leading != 1))
      throw std::invalid_argument(
          "malformed exact chart prescription");
    prescriptions.push_back(json::object{
        {"factor_exact", factor}, {"sign", sign},
        {"multiplicity", multiplicity},
        {"leading_coefficient_sign", leading}});
  }

  json::object canonical{{"center_exact", center},
                         {"scale_exact", scale},
                         {"infinite_radius", infinite}};
  if (const auto* numeric = geometry.if_contains("center_numeric"))
    canonical["center_numeric"] = *numeric;
  if (const auto* numeric = geometry.if_contains("scale_numeric"))
    canonical["scale_numeric"] = *numeric;
  if (infinite) {
    if (const auto* radius = geometry.if_contains("radius_exact");
        radius != nullptr && !radius->is_null() &&
        (!radius->is_string() || radius->as_string() != "Infinity"))
      throw std::invalid_argument(
          "infinite exact chart geometry radius identity must be Infinity");
    // The Wolfram producer historically emitted radius_exact even for an
    // infinite chart.  Canonicalize both its explicit form and an omitted
    // field to the same exact record so composite preparation remains
    // backward compatible without accepting a finite-radius mismatch.
    canonical["radius_exact"] = "Infinity";
  } else {
    canonical["radius_exact"] = required_string(geometry, "radius_exact");
    if (const auto* numeric = geometry.if_contains("radius_numeric"))
      canonical["radius_numeric"] = *numeric;
  }
  canonical["prescriptions"] = std::move(prescriptions);
  return json::serialize(canonical);
}

void validate_first_slice_rational_geometry(const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  const auto exact_specialization = [&](const char* exact_key,
                                        const char* numeric_key) {
    if (const auto* numeric = geometry.if_contains(numeric_key))
      return parse_scalar<ComplexBall>(*numeric);
    return ComplexBall::from_strings(
        Rational(required_string(geometry, exact_key)).str());
  };
  const auto center = exact_specialization("center_exact", "center_numeric");
  const auto scale = exact_specialization("scale_exact", "scale_numeric");
  if (!local_detail::exactly_real(center) ||
      !local_detail::exactly_real(scale) || scale.contains_zero())
    throw std::invalid_argument(
        "native SCC first slice requires certified real center/scale specializations and a scale excluding zero");
  if (geometry.at("infinite_radius").as_bool())
    throw std::invalid_argument(
        "native SCC first slice requires a positive finite rational radius");

  const auto radius = exact_specialization("radius_exact", "radius_numeric");
  if (!local_detail::exactly_real(radius) ||
      !arb_is_positive(acb_realref(radius.raw())))
    throw std::invalid_argument(
        "native SCC first slice requires a provably positive finite radius specialization");
}

struct RetainedCompositeGeometry {
  ChartGeometry chart;
  std::string radius_exact;
  ComplexBall center_numeric;
  ComplexBall scale_numeric;
  std::vector<Prescription> prescriptions;
};

RetainedCompositeGeometry parse_retained_composite_geometry(
    const json::value& raw) {
  const auto& geometry = as_object(raw, "exact chart geometry");
  const auto exact_specialization = [&](const char* exact_key,
                                        const char* numeric_key) {
    if (const auto* numeric = geometry.if_contains(numeric_key))
      return parse_scalar<ComplexBall>(*numeric);
    return ComplexBall::from_strings(
        Rational(required_string(geometry, exact_key)).str());
  };
  RetainedCompositeGeometry retained;
  retained.chart.center_exact = required_string(geometry, "center_exact");
  retained.chart.scale_exact = required_string(geometry, "scale_exact");
  retained.center_numeric = exact_specialization(
      "center_exact", "center_numeric");
  retained.scale_numeric = exact_specialization(
      "scale_exact", "scale_numeric");
  retained.chart.infinite_radius = geometry.at("infinite_radius").as_bool();
  retained.radius_exact = retained.chart.infinite_radius
      ? "Infinity" : required_string(geometry, "radius_exact");
  if (!retained.chart.infinite_radius)
    retained.chart.radius = exact_specialization(
        "radius_exact", "radius_numeric");
  for (const auto& raw_prescription : as_array(
           geometry.at("prescriptions"), "exact chart prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "exact chart prescription");
    retained.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "prescription sign"),
        as_u32(prescription.at("multiplicity"),
               "prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "prescription leading coefficient sign")});
  }
  return retained;
}

std::string canonical_native_scc_capabilities(const json::value& raw) {
  const auto& capabilities = as_object(
      raw, "native SCC chart capabilities");
  const bool identity_v = capabilities.at("identity_v").as_bool();
  const bool identity_gauge = capabilities.at("identity_gauge").as_bool();
  const bool exact_gauge = capabilities.if_contains("exact_gauge") != nullptr
      ? capabilities.at("exact_gauge").as_bool()
      : identity_gauge;
  const bool epsilon_unimodular_v =
      capabilities.if_contains("epsilon_unimodular_v") != nullptr
          ? capabilities.at("epsilon_unimodular_v").as_bool()
          : identity_v;
  return json::serialize(json::object{
      {"regular", capabilities.at("regular").as_bool()},
      {"identity_gauge", identity_gauge},
      {"exact_gauge", exact_gauge},
      {"identity_v", identity_v},
      {"epsilon_unimodular_v", epsilon_unimodular_v},
      {"no_pseudo", capabilities.at("no_pseudo").as_bool()}});
}

TruthValue parse_truth_value(const json::value& value, const char* label) {
  if (!value.is_string())
    throw std::invalid_argument(std::string(label) + " must be yes, no or unknown");
  const auto text = std::string(value.as_string());
  if (text == "yes") return TruthValue::Yes;
  if (text == "no") return TruthValue::No;
  if (text == "unknown") return TruthValue::Unknown;
  throw std::invalid_argument(std::string(label) + " must be yes, no or unknown");
}

ExactSign parse_exact_sign(const json::value& value, const char* label) {
  if (!value.is_string())
    throw std::invalid_argument(
        std::string(label) + " must be negative, zero, positive or unknown");
  const auto text = std::string(value.as_string());
  if (text == "negative") return ExactSign::Negative;
  if (text == "zero") return ExactSign::Zero;
  if (text == "positive") return ExactSign::Positive;
  if (text == "unknown") return ExactSign::Unknown;
  throw std::invalid_argument(
      std::string(label) + " must be negative, zero, positive or unknown");
}

void verify_optional_truth(const json::object& object, const char* key,
                           TruthValue expected) {
  if (const auto* raw = object.if_contains(key);
      raw != nullptr && parse_truth_value(*raw, key) != expected)
    throw std::invalid_argument(
        std::string("contradictory rational exact-tag fact: ") + key);
}

void verify_optional_sign(const json::object& object, const char* key,
                          ExactSign expected) {
  if (const auto* raw = object.if_contains(key);
      raw != nullptr && parse_exact_sign(*raw, key) != expected)
    throw std::invalid_argument(
        std::string("contradictory rational exact-tag fact: ") + key);
}

ExactScalarDescriptor parse_exact_descriptor(const json::value& raw,
                                              const char* label) {
  const auto& object = as_object(raw, label);
  const auto domain = required_string(object, "domain");
  const auto canonical = required_string(object, "canonical");
  if (domain == "rational") {
    auto descriptor = ExactScalarDescriptor::rational(canonical);
    verify_optional_truth(object, "is_zero", descriptor.is_zero);
    verify_optional_truth(object, "is_integer", descriptor.is_integer);
    verify_optional_sign(object, "sign", descriptor.sign);
    if (const auto* raw_specialization = object.if_contains("specialization")) {
      const auto specialization = parse_scalar<ComplexBall>(*raw_specialization);
      if (!acb_equal(specialization.raw(), descriptor.numeric().raw()))
        throw std::invalid_argument(
            std::string(label) + " specialization contradicts exact rational tag");
    }
    return descriptor;
  }

  const auto zero = parse_truth_value(object.at("is_zero"), "is_zero");
  const auto integer = parse_truth_value(object.at("is_integer"), "is_integer");
  const auto sign = parse_exact_sign(object.at("sign"), "sign");
  std::optional<ComplexBall> specialization;
  if (const auto* raw_specialization = object.if_contains("specialization"))
    specialization = parse_scalar<ComplexBall>(*raw_specialization);
  if (domain == "symbolic-rational") {
    std::vector<std::string> symbols;
    for (const auto& symbol : as_array(object.at("symbols"), "tag symbols")) {
      if (!symbol.is_string() || symbol.as_string().empty())
        throw std::invalid_argument("tag symbols must be nonempty strings");
      symbols.emplace_back(symbol.as_string());
    }
    return ExactScalarDescriptor::symbolic(
        canonical, std::move(symbols), zero, integer, sign,
        std::move(specialization));
  }
  if (domain == "algebraic") {
    if (!specialization.has_value())
      throw std::invalid_argument(
          std::string(label) + " algebraic tag requires a specialization");
    return ExactScalarDescriptor::algebraic(
        canonical, zero, integer, sign, std::move(*specialization));
  }
  throw std::invalid_argument(
      std::string(label) + " has unsupported exact domain: " + domain);
}

const char* encode_truth_value(TruthValue value) {
  if (value == TruthValue::Yes) return "yes";
  if (value == TruthValue::No) return "no";
  return "unknown";
}

const char* encode_exact_sign(ExactSign value) {
  if (value == ExactSign::Negative) return "negative";
  if (value == ExactSign::Zero) return "zero";
  if (value == ExactSign::Positive) return "positive";
  return "unknown";
}

json::object encode_exact_descriptor(const ExactScalarDescriptor& descriptor) {
  const char* domain = descriptor.domain == ExactDomain::Rational
      ? "rational"
      : descriptor.domain == ExactDomain::SymbolicRational
          ? "symbolic-rational" : "algebraic";
  json::array symbols;
  for (const auto& symbol : descriptor.symbols) symbols.emplace_back(symbol);
  return json::object{{"domain", domain}, {"canonical", descriptor.canonical},
                      {"symbols", std::move(symbols)},
                      {"is_zero", encode_truth_value(descriptor.is_zero)},
                      {"is_integer", encode_truth_value(descriptor.is_integer)},
                      {"sign", encode_exact_sign(descriptor.sign)},
                      {"has_specialization",
                       descriptor.specialization.has_value()}};
}

struct LocalMetadata {
  ChartGeometry chart;
  ExactScalarDescriptor a;
  ExactScalarDescriptor b;
  std::vector<Prescription> prescriptions;
  std::string checkpoint_identity;
};

LocalMetadata parse_local_metadata(const json::object& metadata) {
  LocalMetadata out;
  const auto& chart = as_object(metadata.at("chart"), "local chart metadata");
  out.chart.center_exact = required_string(chart, "center_exact");
  out.chart.scale_exact = required_string(chart, "scale_exact");
  out.chart.infinite_radius = chart.at("infinite_radius").as_bool();
  if (const auto* radius = chart.if_contains("radius"))
    out.chart.radius = parse_scalar<ComplexBall>(*radius);
  else if (!out.chart.infinite_radius)
    throw std::invalid_argument("finite local chart requires a radius");
  if (!out.chart.infinite_radius &&
      (!local_detail::exactly_real(out.chart.radius) ||
       !arb_is_positive(acb_realref(out.chart.radius.raw()))))
    throw std::invalid_argument(
        "finite local chart radius must be a provably positive real ball");

  const auto& tag = as_object(metadata.at("tag"), "local exact tag");
  out.a = parse_exact_descriptor(tag.at("a"), "local a tag");
  out.b = parse_exact_descriptor(tag.at("b"), "local b tag");
  for (const auto* descriptor : {&out.a, &out.b})
    if (descriptor->specialization.has_value() &&
        !local_detail::exactly_real(*descriptor->specialization))
      throw std::invalid_argument(
          "local exact-tag specialization must be real");
  for (const auto& raw_prescription : as_array(
           metadata.at("prescriptions"), "local prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "local prescription");
    Prescription parsed{
        required_string(prescription, "factor_exact"),
        as_i32(prescription.at("sign"), "prescription sign"),
        as_u32(prescription.at("multiplicity"), "prescription multiplicity"),
        as_i32(prescription.at("leading_coefficient_sign"),
               "prescription leading coefficient sign")};
    if ((parsed.sign != -1 && parsed.sign != 1) ||
        (parsed.leading_coefficient_sign != -1 &&
         parsed.leading_coefficient_sign != 1) ||
        parsed.multiplicity == 0)
      throw std::invalid_argument(
          "malformed analytic-continuation prescription");
    out.prescriptions.push_back(std::move(parsed));
  }
  out.checkpoint_identity = required_string(metadata, "checkpoint_identity");
  return out;
}

template <typename Scalar>
void verify_tag_binding(const ExactScalarDescriptor& descriptor,
                        const Scalar& target, const char* label);

template <>
void verify_tag_binding<Rational>(const ExactScalarDescriptor& descriptor,
                                  const Rational& target,
                                  const char* label) {
  if (descriptor.domain != ExactDomain::Rational ||
      !(Rational(descriptor.canonical) == target))
    throw std::invalid_argument(
        std::string(label) + " exact tag does not equal recurrence target");
}

template <>
void verify_tag_binding<SymbolicRational>(
    const ExactScalarDescriptor& descriptor,
    const SymbolicRational& target, const char* label) {
  if (descriptor.domain == ExactDomain::Algebraic ||
      !(SymbolicRational(descriptor.canonical) == target))
    throw std::invalid_argument(
        std::string(label) + " exact tag does not equal recurrence target");
}

template <>
void verify_tag_binding<ComplexBall>(const ExactScalarDescriptor& descriptor,
                                     const ComplexBall& target,
                                     const char* label) {
  if (!descriptor.specialization.has_value())
    throw std::invalid_argument(
        std::string(label) + " exact tag needs a numeric specialization");
  if (!acb_equal(descriptor.specialization->raw(), target.raw()))
    throw std::invalid_argument(
        std::string(label) + " specialization does not equal recurrence target");
}

struct ChartStats {
  std::uint64_t runs = 0;
  std::uint64_t local_runs = 0;
  double prepare_parse_ms = 0.0;
  double run_parse_ms = 0.0;
  double kernel_ms = 0.0;
  double local_run_parse_ms = 0.0;
  double local_kernel_ms = 0.0;
};

class StoredLocalBase;
class StoredTilePlan;
class StoredTransportArmState;

// Equation ownership is deliberately independent of the current primitive
// chart implementation.  Composite SCC owners can implement this interface
// once they retain one full original-system q/C certificate, without changing
// StoredLocal or its checkpoint provenance shape.
class PhysicalEquationOwnerBase {
 public:
  virtual ~PhysicalEquationOwnerBase() = default;
  virtual const std::string& equation_owner_handle() const = 0;
  virtual const std::string& equation_operator_identity() const = 0;
  virtual const char* equation_owner_kind() const = 0;
  virtual const char* equation_scalar_domain() const = 0;
  virtual std::shared_ptr<const void> physical_ode_erased() const = 0;
  virtual const std::string& physical_payload_identity() const = 0;
  virtual const std::string& physical_payload_record() const = 0;
  virtual const std::string& owner_signature_identity() const = 0;
  virtual std::optional<std::pair<
      FiniteLaurentVector<ComplexBall>, std::string>>
  normalize_acb_matching_vector(
      const FiniteLaurentVector<ComplexBall>& physical) const {
    (void)physical;
    return std::nullopt;
  }
  virtual std::optional<FiniteLaurentMatrix<ComplexBall>>
  right_normalize_acb_matching_basis(
      const FiniteLaurentMatrix<ComplexBall>& left_normalized) const {
    (void)left_normalized;
    return std::nullopt;
  }
  virtual std::optional<FiniteLaurentVector<ComplexBall>>
  denormalize_acb_matching_weights(
      const FiniteLaurentVector<ComplexBall>& normalized) const {
    (void)normalized;
    return std::nullopt;
  }
};

class PreparedChartBase : public PhysicalEquationOwnerBase {
 public:
  PreparedChartBase(std::string handle, std::string key,
                    std::string exact_identity, std::string signature,
                    std::optional<std::string> geometry_record,
                    std::optional<std::string> principal_matrix_record,
                    std::optional<std::string> native_scc_capabilities,
                    std::optional<std::string>
                        regular_value_relative_accuracy_max_exact,
                    SCCCertificate scc, double prepare_parse_ms)
      : handle_(std::move(handle)), key_(std::move(key)),
        exact_identity_(std::move(exact_identity)),
        signature_(std::move(signature)),
        geometry_record_(std::move(geometry_record)),
        principal_matrix_record_(std::move(principal_matrix_record)),
        native_scc_capabilities_(std::move(native_scc_capabilities)),
        regular_value_relative_accuracy_max_exact_(
            std::move(regular_value_relative_accuracy_max_exact)),
        scc_(std::move(scc)),
        prepare_parse_ms_(prepare_parse_ms) {}
  virtual ~PreparedChartBase() = default;

  virtual json::object solve(const json::object& run, int output_digits) = 0;
  virtual std::shared_ptr<StoredLocalBase> solve_local(
      const std::string& local_handle, const json::object& run,
      const json::object& metadata,
      std::shared_ptr<PhysicalEquationOwnerBase> equation_owner) = 0;
  virtual std::uint32_t dimension() const = 0;
  virtual std::int32_t frame_base() const = 0;
  virtual std::uint32_t frame_width() const = 0;
  virtual const char* d0_inverse_mode() const = 0;
  virtual std::string regular_value_tail_proxy_max_exact() const = 0;
  virtual ChartStats stats() const = 0;

  const std::string& equation_owner_handle() const override {
    return handle_;
  }
  const std::string& equation_operator_identity() const override {
    return exact_identity_;
  }
  const char* equation_owner_kind() const override {
    return "prepared-chart";
  }

  const std::string& handle() const { return handle_; }
  const std::string& key() const { return key_; }
  const std::string& exact_identity() const { return exact_identity_; }
  const std::string& signature() const { return signature_; }
  const std::optional<std::string>& geometry_record() const {
    return geometry_record_;
  }
  const std::optional<std::string>& principal_matrix_record() const {
    return principal_matrix_record_;
  }
  const std::optional<std::string>& native_scc_capabilities() const {
    return native_scc_capabilities_;
  }
  const std::optional<std::string>&
  regular_value_relative_accuracy_max_exact() const {
    return regular_value_relative_accuracy_max_exact_;
  }
  const SCCCertificate& scc() const { return scc_; }

 protected:
  std::string handle_;
  std::string key_;
  std::string exact_identity_;
  std::string signature_;
  std::optional<std::string> geometry_record_;
  std::optional<std::string> principal_matrix_record_;
  std::optional<std::string> native_scc_capabilities_;
  std::optional<std::string> regular_value_relative_accuracy_max_exact_;
  SCCCertificate scc_;
  double prepare_parse_ms_ = 0.0;
};

std::recursive_mutex& symbolic_run_mutex() {
  static std::recursive_mutex mutex;
  return mutex;
}

template <typename Scalar, typename Object, typename... Arguments>
std::shared_ptr<Object> make_retained_typed_shared(Arguments&&... arguments) {
  if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
    // SymbolicRational destruction touches the process-global FLINT context
    // and live-object count.  The deleter travels with the shared_ptr control
    // block through base conversions and dynamic casts, so releases, session
    // close, library reset, and delayed in-flight holders all serialize.  A
    // recursive mutex is required because destroying a composite can release
    // its last typed chart pointers under this same deleter.
    std::lock_guard<std::recursive_mutex> construction_lock(
        symbolic_run_mutex());
    return std::shared_ptr<Object>(
        new Object(std::forward<Arguments>(arguments)...),
        [](Object* object) {
          std::lock_guard<std::recursive_mutex> lock(symbolic_run_mutex());
          delete object;
        });
  } else {
    return std::make_shared<Object>(
        std::forward<Arguments>(arguments)...);
  }
}

struct AcbExecutionState {
  std::mutex mutex;
  std::condition_variable changed;
  std::optional<slong> precision_bits;
  std::size_t active = 0;
};

AcbExecutionState& acb_execution_state() {
  static AcbExecutionState state;
  return state;
}

class AcbPrecisionLease {
 public:
  explicit AcbPrecisionLease(slong precision_bits)
      : state_(acb_execution_state()) {
    std::unique_lock<std::mutex> lock(state_.mutex);
    state_.changed.wait(lock, [&] {
      return state_.active == 0 || state_.precision_bits == precision_bits;
    });
    if (state_.active == 0) state_.precision_bits = precision_bits;
    ++state_.active;
  }
  AcbPrecisionLease(const AcbPrecisionLease&) = delete;
  AcbPrecisionLease& operator=(const AcbPrecisionLease&) = delete;
  ~AcbPrecisionLease() {
    std::lock_guard<std::mutex> lock(state_.mutex);
    --state_.active;
    if (state_.active == 0) state_.changed.notify_all();
  }

 private:
  AcbExecutionState& state_;
};

json::object encode_epsilon_vector(const EpsilonVector& vector, int digits) {
  json::array coefficients;
  coefficients.reserve(vector.coefficients.size());
  for (const auto& coefficient : vector.coefficients)
    coefficients.push_back(encode_scalar(coefficient, digits));
  json::object encoded{{"min", vector.epsilon.min_power},
                       {"max", vector.epsilon.complete_max},
                       {"dimension", vector.dimension},
                       {"coefficients", std::move(coefficients)}};
  if (!vector.error.empty()) {
    json::array upper;
    upper.reserve(vector.error.absolute.size());
    for (const auto& bound : vector.error.absolute)
      upper.push_back(bound.approximate_upper());
    const char* guarantee = vector.error.guarantee == ErrorGuarantee::Certified
        ? "certified"
        : vector.error.guarantee == ErrorGuarantee::Advisory
            ? "advisory" : "none";
    encoded["error"] = json::object{
        {"min", vector.error.frame.min_power},
        {"max", vector.error.frame.complete_max},
        {"guarantee", guarantee},
        {"absolute_upper_approx", std::move(upper)},
        {"bound_encoding", "approximate-double"},
        {"provenance", vector.error.provenance}};
  }
  return encoded;
}

const char* tail_majorant_status_name(TailMajorantStatus status) {
  switch (status) {
    case TailMajorantStatus::Certified:
      return "certified";
    case TailMajorantStatus::Inconclusive:
      return "inconclusive";
    case TailMajorantStatus::Unsupported:
      return "unsupported";
  }
  throw std::logic_error("unknown tail-majorant status");
}

RegularTaylorTailModelResult unavailable_tail_model(std::string detail) {
  return {TailMajorantStatus::Unsupported, std::nullopt,
          std::move(detail)};
}

struct TailModelCheckpointMarker {
  std::string saved_status;
  bool attached_before_save = false;
};

json::object encode_tail_model_status(
    const RegularTaylorTailModelResult& result) {
  json::object encoded{
      {"capability", kRegularTailMajorantCapability},
      {"status", tail_majorant_status_name(result.status)},
      {"attached", result.model.has_value()},
      {"detail", result.detail},
      {"checkpoint_serialized", result.model.has_value()}};
  if (result.model.has_value()) {
    encoded["operator_identity"] = result.model->operator_identity;
    encoded["local_checkpoint_identity"] =
        result.model->local_checkpoint_identity;
    encoded["epsilon"] = json::object{
        {"min", result.model->epsilon.min_power},
        {"max", result.model->epsilon.complete_max}};
    encoded["taylor_complete_max"] =
        result.model->taylor_complete_max;
    encoded["provenance"] = result.model->provenance;
  }
  return encoded;
}

std::optional<std::string> parse_certified_tail_witness(
    const json::object& request) {
  const auto* raw_options = request.if_contains("options");
  if (raw_options == nullptr) return std::nullopt;
  const auto& options = as_object(*raw_options, "local evaluation options");
  const auto* raw_radius =
      options.if_contains("certified_tail_radius_exact");
  if (raw_radius == nullptr || raw_radius->is_null()) return std::nullopt;
  if (!raw_radius->is_string() || raw_radius->as_string().empty())
    throw std::invalid_argument(
        "certified tail radius must be a nonempty exact rational string");
  const std::string radius(raw_radius->as_string());
  const Rational parsed(radius);
  if (parsed.sign() <= 0)
    throw std::invalid_argument("certified tail radius must be positive");
  return parsed.str();
}

json::object encode_point_tail_certificate(
    const RegularTaylorPointTailCertificate& certificate,
    const std::string& witness_radius) {
  json::object result{
      {"capability", kRegularTailMajorantCapability},
      {"requested", true},
      {"status", tail_majorant_status_name(certificate.status)},
      {"witness_radius_exact", witness_radius},
      {"detail", certificate.detail}};
  if (certificate.disk.status == TailMajorantStatus::Certified) {
    result["q_lower_approx"] =
        certificate.disk.q_lower.approximate_upper();
    result["ode_norm_upper_approx"] =
        certificate.disk.ode_norm_upper.approximate_upper();
    result["bound_encoding"] = "approximate-double-diagnostics";
  }
  return result;
}

EvaluationOptions parse_local_evaluation_options(
    const json::object& request, bool default_tail_estimate) {
  EvaluationOptions options;
  options.compute_tail_estimate = default_tail_estimate;
  if (const auto* raw_options = request.if_contains("options")) {
    const auto& object = as_object(*raw_options, "local evaluation options");
    if (const auto* raw_sign = object.if_contains("imaginary_sign");
        raw_sign != nullptr && !raw_sign->is_null()) {
      const auto sign = as_i32(*raw_sign, "imaginary sign");
      if (sign != -1 && sign != 1)
        throw std::invalid_argument("imaginary sign must be +1 or -1");
      options.imaginary_sign = sign;
    }
    if (const auto* reduction = object.if_contains("t_order_reduction"))
      options.t_order_reduction = as_u32(
          *reduction, "Taylor-order reduction");
    if (const auto* tail = object.if_contains("tail_estimate"))
      options.compute_tail_estimate = tail->as_bool();
  }
  return options;
}

RealEvaluationPoint parse_local_evaluation_point(const json::object& request) {
  const auto& point_object = as_object(
      request.at("point"), "local evaluation point");
  return RealEvaluationPoint::rational(
      required_string(point_object, "exact"));
}

EpsilonWindow parse_epsilon_window(const json::object& object,
                                   const char* label) {
  EpsilonWindow window{as_i32(object.at("min"), label),
                       as_i32(object.at("max"), label)};
  (void)window.width();
  return window;
}

std::size_t checked_flat_count(std::size_t left, std::size_t right,
                               const char* label) {
  if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right)
    throw std::overflow_error(std::string(label) + " size overflow");
  return left * right;
}

EpsilonMatrix parse_epsilon_matrix(const json::value& raw,
                                   const char* label) {
  const auto& object = as_object(raw, label);
  EpsilonMatrix matrix;
  matrix.epsilon = parse_epsilon_window(object, label);
  matrix.dimension = as_u32(object.at("dimension"), label);
  if (matrix.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " dimension must be positive");
  const auto matrix_size = checked_flat_count(
      matrix.dimension, matrix.dimension, label);
  const auto expected = checked_flat_count(
      matrix.epsilon.width(), matrix_size, label);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient tensor has the wrong size");
  matrix.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    matrix.coefficients.push_back(parse_scalar<ComplexBall>(coefficient));
  return matrix;
}

EpsilonVector parse_epsilon_vector(const json::value& raw,
                                   const char* label) {
  const auto& object = as_object(raw, label);
  EpsilonVector vector;
  vector.epsilon = parse_epsilon_window(object, label);
  vector.dimension = as_u32(object.at("dimension"), label);
  if (vector.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " dimension must be positive");
  const auto expected = checked_flat_count(
      vector.epsilon.width(), vector.dimension, label);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient tensor has the wrong size");
  vector.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    vector.coefficients.push_back(parse_scalar<ComplexBall>(coefficient));
  return vector;
}

json::array encode_magnitude_diagnostics(
    const std::vector<Magnitude>& magnitudes) {
  json::array encoded;
  encoded.reserve(magnitudes.size());
  for (const auto& magnitude : magnitudes)
    encoded.push_back(magnitude.is_finite()
                          ? json::value(magnitude.approximate_upper())
                          : json::value(nullptr));
  return encoded;
}

struct StoredLocalStats {
  std::uint64_t evaluations = 0;
  std::uint64_t residual_certifications = 0;
  std::uint64_t endpoint_limits = 0;
  std::uint64_t line_integrations = 0;
  double evaluate_ms = 0.0;
  double residual_certify_ms = 0.0;
  double endpoint_limit_ms = 0.0;
  double line_integration_ms = 0.0;
  double create_parse_ms = 0.0;
  double create_kernel_ms = 0.0;
  std::size_t coefficient_count = 0;
  std::uint64_t tail_certificate_requests = 0;
  std::uint64_t tail_certificate_certified = 0;
  std::uint64_t tail_certificate_inconclusive = 0;
  std::uint64_t tail_certificate_unsupported = 0;
};

struct NativeLocalDiagnostics {
  std::int32_t top_valid = kCompleteInfinity;
  double parse_ms = 0.0;
  double kernel_ms = 0.0;
  std::uint64_t pseudo_hits = 0;
  std::uint64_t pseudo_compensations = 0;
  std::uint32_t max_pseudo_depth = 0;
  bool pseudo_value_certified = true;
};

json::object checkpoint_ball_record(const ComplexBall& value) {
  const auto dump = checkpoint::dump_complex_ball_exact(value);
  return json::object{{"real", dump.real}, {"imaginary", dump.imaginary}};
}

ComplexBall parse_checkpoint_ball(const json::value& raw,
                                  const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"real", "imaginary"}, label);
  return checkpoint::load_complex_ball_exact(
      {required_string(object, "real"),
       required_string(object, "imaginary")});
}

Magnitude parse_checkpoint_magnitude(const json::value& raw,
                                     const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) +
                                " must be an exact dump string");
  auto result = Magnitude::from_exact_dump(std::string(raw.as_string()));
  if (!result.is_finite())
    throw std::invalid_argument(std::string(label) +
                                " must be finite");
  return result;
}

template <typename Scalar>
json::value checkpoint_scalar_record(const Scalar& value);

template <>
json::value checkpoint_scalar_record<Rational>(const Rational& value) {
  return json::string(value.str());
}

template <>
json::value checkpoint_scalar_record<ComplexBall>(const ComplexBall& value) {
  return checkpoint_ball_record(value);
}

template <typename Scalar>
Scalar parse_checkpoint_scalar(const json::value& raw, const char* label);

template <>
Rational parse_checkpoint_scalar<Rational>(const json::value& raw,
                                            const char* label) {
  if (!raw.is_string())
    throw std::invalid_argument(std::string(label) +
                                " must be an exact rational string");
  return Rational(std::string(raw.as_string()));
}

template <>
ComplexBall parse_checkpoint_scalar<ComplexBall>(const json::value& raw,
                                                  const char* label) {
  return parse_checkpoint_ball(raw, label);
}

json::object checkpoint_exact_descriptor_record(
    const ExactScalarDescriptor& descriptor) {
  auto output = encode_exact_descriptor(descriptor);
  output.erase("has_specialization");
  output["specialization"] = descriptor.specialization.has_value()
      ? json::value(checkpoint_ball_record(*descriptor.specialization))
      : json::value(nullptr);
  return output;
}

ExactScalarDescriptor parse_checkpoint_exact_descriptor(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object,
      {"domain", "canonical", "symbols", "is_zero", "is_integer",
       "sign", "specialization"}, label);
  const auto domain = required_string(object, "domain");
  const auto canonical = required_string(object, "canonical");
  if (canonical.empty())
    throw std::invalid_argument(std::string(label) +
                                " canonical form is empty");
  std::vector<std::string> symbols;
  for (const auto& raw_symbol : as_array(object.at("symbols"), label)) {
    if (!raw_symbol.is_string() || raw_symbol.as_string().empty())
      throw std::invalid_argument(std::string(label) +
                                  " symbols must be nonempty strings");
    symbols.emplace_back(raw_symbol.as_string());
  }
  const auto zero = parse_truth_value(object.at("is_zero"), "is_zero");
  const auto integer = parse_truth_value(object.at("is_integer"),
                                         "is_integer");
  const auto sign = parse_exact_sign(object.at("sign"), "sign");
  std::optional<ComplexBall> specialization;
  if (!object.at("specialization").is_null())
    specialization = parse_checkpoint_ball(object.at("specialization"),
                                           "exact tag specialization");

  if (domain == "rational") {
    if (!symbols.empty())
      throw std::invalid_argument(std::string(label) +
                                  " rational tag cannot name symbols");
    auto result = ExactScalarDescriptor::rational(canonical);
    if (result.is_zero != zero || result.is_integer != integer ||
        result.sign != sign)
      throw std::invalid_argument(std::string(label) +
                                  " rational facts contradict its value");
    fmpq_t exact;
    fmpq_init(exact);
    const auto parse_status = fmpq_set_str(exact, canonical.c_str(), 10);
    if (parse_status == 0) fmpq_canonicalise(exact);
    const bool consistent = parse_status == 0 && specialization.has_value() &&
        acb_contains_fmpq(specialization->raw(), exact);
    fmpq_clear(exact);
    if (!consistent)
      throw std::invalid_argument(std::string(label) +
                                  " rational specialization is missing or inconsistent");
    result.specialization = std::move(specialization);
    return result;
  }
  if (domain == "symbolic-rational")
    return ExactScalarDescriptor::symbolic(
        canonical, std::move(symbols), zero, integer, sign,
        std::move(specialization));
  if (domain == "algebraic") {
    if (!symbols.empty())
      throw std::invalid_argument(std::string(label) +
                                  " algebraic tag cannot name regulator symbols");
    if (!specialization.has_value())
      throw std::invalid_argument(std::string(label) +
                                  " algebraic tag lost its specialization");
    return ExactScalarDescriptor::algebraic(
        canonical, zero, integer, sign, std::move(*specialization));
  }
  throw std::invalid_argument(std::string(label) +
                              " has an unsupported exact domain");
}

const char* checkpoint_error_guarantee_name(ErrorGuarantee guarantee) {
  if (guarantee == ErrorGuarantee::Certified) return "certified";
  if (guarantee == ErrorGuarantee::Advisory) return "advisory";
  return "none";
}

json::object checkpoint_error_envelope_record(
    const ErrorEnvelope& envelope) {
  json::array absolute;
  absolute.reserve(envelope.absolute.size());
  for (const auto& magnitude : envelope.absolute)
    absolute.emplace_back(magnitude.dump_exact());
  return json::object{
      {"frame", json::object{{"min", envelope.frame.min_power},
                              {"max", envelope.frame.complete_max}}},
      {"guarantee", checkpoint_error_guarantee_name(envelope.guarantee)},
      {"absolute_exact", std::move(absolute)},
      {"provenance", envelope.provenance}};
}

ErrorEnvelope parse_checkpoint_error_envelope(const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint error envelope");
  require_exact_keys(object,
      {"frame", "guarantee", "absolute_exact", "provenance"},
      "checkpoint error envelope");
  const auto& frame = as_object(object.at("frame"),
                                "checkpoint error frame");
  require_exact_keys(frame, {"min", "max"}, "checkpoint error frame");
  ErrorEnvelope result;
  result.frame = {as_i32(frame.at("min"), "checkpoint error minimum"),
                  as_i32(frame.at("max"), "checkpoint error maximum")};
  (void)result.frame.width();
  const auto guarantee = required_string(object, "guarantee");
  result.guarantee = guarantee == "none" ? ErrorGuarantee::None
      : guarantee == "advisory" ? ErrorGuarantee::Advisory
      : guarantee == "certified" ? ErrorGuarantee::Certified
      : throw std::invalid_argument(
            "checkpoint error guarantee is unsupported");
  for (const auto& magnitude : as_array(
           object.at("absolute_exact"), "checkpoint error magnitudes")) {
    if (!magnitude.is_string())
      throw std::invalid_argument(
          "checkpoint error magnitudes must be exact dump strings");
    result.absolute.push_back(Magnitude::from_exact_dump(
        std::string(magnitude.as_string())));
  }
  if (!object.at("provenance").is_string())
    throw std::invalid_argument(
        "checkpoint error provenance must be a string");
  result.provenance = std::string(object.at("provenance").as_string());
  if (!result.absolute.empty() &&
      result.absolute.size() != result.frame.width())
    throw std::invalid_argument(
        "checkpoint error envelope width is inconsistent");
  if (result.absolute.empty() && result.guarantee != ErrorGuarantee::None)
    throw std::invalid_argument(
        "empty checkpoint error envelope cannot claim a guarantee");
  return result;
}

template <typename Scalar>
json::object checkpoint_epsilon_frame_record(
    const EpsilonFrame<Scalar>& frame) {
  json::array coefficients;
  coefficients.reserve(frame.coefficients().size());
  for (const auto& coefficient : frame.coefficients())
    coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
  return json::object{{"min", frame.min_power()},
                      {"max", frame.complete_max()},
                      {"coefficients", std::move(coefficients)}};
}

template <typename Scalar>
EpsilonFrame<Scalar> parse_checkpoint_epsilon_frame(
    const json::value& raw, const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"min", "max", "coefficients"}, label);
  EpsilonWindow window{as_i32(object.at("min"), label),
                       as_i32(object.at("max"), label)};
  const auto& raw_coefficients = as_array(object.at("coefficients"), label);
  if (raw_coefficients.size() != window.width())
    throw std::invalid_argument(std::string(label) +
                                " coefficient count is inconsistent");
  std::vector<Scalar> coefficients;
  coefficients.reserve(raw_coefficients.size());
  for (const auto& coefficient : raw_coefficients)
    coefficients.push_back(parse_checkpoint_scalar<Scalar>(coefficient,
                                                             label));
  return EpsilonFrame<Scalar>(window, std::move(coefficients));
}

template <typename Scalar>
json::array checkpoint_frame_vector_record(
    const FiniteLaurentVector<Scalar>& frames) {
  json::array output;
  output.reserve(frames.size());
  for (const auto& frame : frames)
    output.push_back(checkpoint_epsilon_frame_record(frame));
  return output;
}

template <typename Scalar>
FiniteLaurentVector<Scalar> parse_checkpoint_frame_vector(
    const json::value& raw, std::size_t expected_size, const char* label) {
  const auto& values = as_array(raw, label);
  if (values.size() != expected_size)
    throw std::invalid_argument(std::string(label) +
                                " dimension is inconsistent");
  FiniteLaurentVector<Scalar> output;
  output.reserve(values.size());
  for (const auto& value : values)
    output.push_back(parse_checkpoint_epsilon_frame<Scalar>(value, label));
  return output;
}

json::array checkpoint_exact_laurent_matrix_record(
    const ExactLaurentMatrix<Rational>& matrix) {
  json::array rows;
  rows.reserve(matrix.size());
  for (const auto& row : matrix) {
    json::array encoded_row;
    encoded_row.reserve(row.size());
    for (const auto& polynomial : row) {
      json::array terms;
      terms.reserve(polynomial.terms().size());
      for (const auto& [power, coefficient] : polynomial.terms())
        terms.push_back(json::object{{"power", power},
                                     {"coefficient", coefficient.str()}});
      encoded_row.push_back(std::move(terms));
    }
    rows.push_back(std::move(encoded_row));
  }
  return rows;
}

ExactLaurentMatrix<Rational> parse_checkpoint_exact_laurent_matrix(
    const json::value& raw, std::uint32_t dimension, const char* label) {
  const auto& rows = as_array(raw, label);
  if (rows.size() != dimension)
    throw std::invalid_argument(std::string(label) +
                                " row count differs from its dimension");
  ExactLaurentMatrix<Rational> matrix;
  matrix.reserve(rows.size());
  for (const auto& raw_row : rows) {
    const auto& row = as_array(raw_row, label);
    if (row.size() != dimension)
      throw std::invalid_argument(std::string(label) +
                                  " is not a square dimension-by-dimension matrix");
    std::vector<ExactLaurentPolynomial<Rational>> parsed_row;
    parsed_row.reserve(row.size());
    for (const auto& raw_entry : row) {
      ExactLaurentPolynomial<Rational> polynomial;
      std::optional<std::int32_t> previous_power;
      for (const auto& raw_term : as_array(raw_entry, label)) {
        const auto& term = as_object(raw_term, label);
        require_exact_keys(term, {"power", "coefficient"}, label);
        const auto power = as_i32(term.at("power"), label);
        if (previous_power.has_value() && power <= *previous_power)
          throw std::invalid_argument(std::string(label) +
                                      " terms are not in strict power order");
        previous_power = power;
        auto coefficient = parse_checkpoint_scalar<Rational>(
            term.at("coefficient"), label);
        if (coefficient.is_zero())
          throw std::invalid_argument(std::string(label) +
                                      " explicitly stores a structural-zero term");
        polynomial.add_term(power, std::move(coefficient));
      }
      parsed_row.push_back(std::move(polynomial));
    }
    matrix.push_back(std::move(parsed_row));
  }
  return matrix;
}

json::object checkpoint_saturation_diagnostics_record(
    const EpsilonLatticeSaturationDiagnostics<Rational>& diagnostics) {
  json::array valuations;
  for (const auto value : diagnostics.initial_column_valuations)
    valuations.push_back(value);
  json::array shifts;
  for (const auto value : diagnostics.initial_column_shifts)
    shifts.push_back(value);
  json::array actions;
  for (const auto& action : diagnostics.actions) {
    json::array relation;
    for (const auto& value : action.null_relation)
      relation.emplace_back(value.str());
    actions.push_back(json::object{
        {"leading_rank_before", action.leading_rank_before},
        {"target_column", action.target_column},
        {"null_relation", std::move(relation)}});
  }
  return json::object{
      {"initial_column_valuations", std::move(valuations)},
      {"initial_column_shifts", std::move(shifts)},
      {"normalized_determinant",
       checkpoint_epsilon_frame_record(diagnostics.normalized_determinant)},
      {"normalized_determinant_valuation",
       diagnostics.normalized_determinant_valuation},
      {"initial_leading_rank", diagnostics.initial_leading_rank},
      {"final_leading_rank", diagnostics.final_leading_rank},
      {"actions", std::move(actions)}};
}

EpsilonLatticeSaturationDiagnostics<Rational>
parse_checkpoint_saturation_diagnostics(const json::value& raw,
                                        std::uint32_t dimension) {
  const auto& object = as_object(raw, "checkpoint saturation diagnostics");
  require_exact_keys(
      object,
      {"initial_column_valuations", "initial_column_shifts",
       "normalized_determinant", "normalized_determinant_valuation",
       "initial_leading_rank", "final_leading_rank", "actions"},
      "checkpoint saturation diagnostics");
  std::vector<std::int32_t> valuations;
  std::vector<std::int32_t> shifts;
  for (const auto& value : as_array(object.at("initial_column_valuations"),
                                     "checkpoint saturation valuations"))
    valuations.push_back(
        as_i32(value, "checkpoint saturation valuation"));
  for (const auto& value : as_array(object.at("initial_column_shifts"),
                                     "checkpoint saturation shifts"))
    shifts.push_back(
        as_i32(value, "checkpoint saturation shift"));
  if (valuations.size() != dimension || shifts.size() != dimension)
    throw std::invalid_argument(
        "checkpoint saturation column metadata differs from its dimension");
  auto normalized_determinant = parse_checkpoint_epsilon_frame<Rational>(
      object.at("normalized_determinant"),
      "checkpoint normalized determinant");
  const auto normalized_determinant_valuation = as_i32(
      object.at("normalized_determinant_valuation"),
      "checkpoint normalized determinant valuation");
  const auto leading = finite_laurent_leading_power(
      normalized_determinant,
      "checkpoint normalized determinant valuation");
  if (!leading.has_value() ||
      *leading != normalized_determinant_valuation)
    throw std::invalid_argument(
        "checkpoint normalized determinant valuation is inconsistent");
  const auto initial_leading_rank = static_cast<std::size_t>(as_u64(
      object.at("initial_leading_rank"),
      "checkpoint initial leading rank"));
  const auto final_leading_rank = static_cast<std::size_t>(as_u64(
      object.at("final_leading_rank"), "checkpoint final leading rank"));
  if (initial_leading_rank > dimension || final_leading_rank != dimension ||
      initial_leading_rank > final_leading_rank)
    throw std::invalid_argument(
        "checkpoint saturation leading ranks are inconsistent");
  std::vector<EpsilonLatticeSaturationAction<Rational>> actions;
  for (const auto& raw_action : as_array(object.at("actions"),
                                          "checkpoint saturation actions")) {
    const auto& action = as_object(raw_action,
                                   "checkpoint saturation action");
    require_exact_keys(action,
        {"leading_rank_before", "target_column", "null_relation"},
        "checkpoint saturation action");
    EpsilonLatticeSaturationAction<Rational> parsed;
    parsed.leading_rank_before = static_cast<std::size_t>(as_u64(
        action.at("leading_rank_before"),
        "checkpoint saturation action rank"));
    parsed.target_column = static_cast<std::size_t>(as_u64(
        action.at("target_column"),
        "checkpoint saturation target column"));
    for (const auto& value : as_array(action.at("null_relation"),
                                      "checkpoint saturation null relation"))
      parsed.null_relation.push_back(parse_checkpoint_scalar<Rational>(
          value, "checkpoint saturation null relation coefficient"));
    if (parsed.leading_rank_before >= dimension ||
        parsed.target_column >= dimension ||
        parsed.null_relation.size() != dimension)
      throw std::invalid_argument(
          "checkpoint saturation action is outside its exact dimension");
    actions.push_back(std::move(parsed));
  }
  if (normalized_determinant_valuation < 0 ||
      actions.size() !=
          static_cast<std::size_t>(normalized_determinant_valuation))
    throw std::invalid_argument(
        "checkpoint saturation action count does not reproduce its determinant valuation");
  return {std::move(valuations), std::move(shifts),
          std::move(normalized_determinant),
          normalized_determinant_valuation, initial_leading_rank,
          final_leading_rank, std::move(actions)};
}

json::object checkpoint_epsilon_vector_record(const EpsilonVector& value) {
  json::array coefficients;
  coefficients.reserve(value.coefficients.size());
  for (const auto& coefficient : value.coefficients)
    coefficients.push_back(checkpoint_ball_record(coefficient));
  return json::object{
      {"epsilon", json::object{{"min", value.epsilon.min_power},
                                {"max", value.epsilon.complete_max}}},
      {"dimension", value.dimension},
      {"coefficients", std::move(coefficients)},
      {"error", checkpoint_error_envelope_record(value.error)}};
}

EpsilonVector parse_checkpoint_epsilon_vector(const json::value& raw,
                                               const char* label) {
  const auto& object = as_object(raw, label);
  require_exact_keys(object, {"epsilon", "dimension", "coefficients", "error"},
                     label);
  const auto& epsilon = as_object(object.at("epsilon"), label);
  require_exact_keys(epsilon, {"min", "max"}, label);
  EpsilonVector result;
  result.epsilon = {as_i32(epsilon.at("min"), label),
                    as_i32(epsilon.at("max"), label)};
  const auto width = result.epsilon.width();
  result.dimension = as_u32(object.at("dimension"), label);
  if (result.dimension == 0)
    throw std::invalid_argument(std::string(label) +
                                " has zero dimension");
  if (width > std::numeric_limits<std::size_t>::max() /
                  static_cast<std::size_t>(result.dimension))
    throw std::invalid_argument(std::string(label) +
                                " tensor size overflows size_t");
  const auto expected =
      width * static_cast<std::size_t>(result.dimension);
  const auto& coefficients = as_array(object.at("coefficients"), label);
  if (coefficients.size() != expected)
    throw std::invalid_argument(std::string(label) +
                                " coefficient count differs from its tensor");
  result.coefficients.reserve(expected);
  for (const auto& coefficient : coefficients)
    result.coefficients.push_back(parse_checkpoint_ball(coefficient, label));
  result.error = parse_checkpoint_error_envelope(object.at("error"));
  return result;
}

template <typename Scalar>
json::object checkpoint_local_analytic_metadata_record(
    const LocalSolution<Scalar>& solution) {
  json::array sectors;
  sectors.reserve(solution.sectors.size());
  for (const auto& sector : solution.sectors)
    sectors.push_back(json::object{
        {"a", checkpoint_exact_descriptor_record(sector.a)},
        {"b", checkpoint_exact_descriptor_record(sector.b)},
        {"log_power", sector.log_power}});
  json::array prescriptions;
  prescriptions.reserve(solution.prescriptions.size());
  for (const auto& prescription : solution.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"schema", "diffexp2-exact-local-analytic-metadata-v2"},
      {"chart", json::object{
          {"center_exact", solution.chart.center_exact},
          {"scale_exact", solution.chart.scale_exact},
          {"radius_exact_ball", checkpoint_ball_record(solution.chart.radius)},
          {"infinite_radius", solution.chart.infinite_radius}}},
      {"sectors", std::move(sectors)},
      {"prescriptions", std::move(prescriptions)}};
}

template <typename Scalar>
json::object checkpoint_local_solution_record(
    const LocalSolution<Scalar>& solution) {
  json::array sectors;
  sectors.reserve(solution.sectors.size());
  for (const auto& sector : solution.sectors) {
    json::array coefficients;
    coefficients.reserve(sector.coefficients.size());
    for (const auto& coefficient : sector.coefficients)
      coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
    sectors.push_back(json::object{
        {"a", checkpoint_exact_descriptor_record(sector.a)},
        {"b", checkpoint_exact_descriptor_record(sector.b)},
        {"log_power", sector.log_power},
        {"coefficients", std::move(coefficients)}});
  }
  json::array prescriptions;
  prescriptions.reserve(solution.prescriptions.size());
  for (const auto& prescription : solution.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"chart", json::object{
          {"center_exact", solution.chart.center_exact},
          {"scale_exact", solution.chart.scale_exact},
          {"radius_exact_ball", checkpoint_ball_record(solution.chart.radius)},
          {"infinite_radius", solution.chart.infinite_radius}}},
      {"epsilon", json::object{{"min", solution.epsilon.min_power},
                                {"max", solution.epsilon.complete_max}}},
      {"taylor_complete_max", solution.taylor_complete_max},
      {"dimension", solution.dimension},
      {"sectors", std::move(sectors)},
      {"prescriptions", std::move(prescriptions)},
      {"error", checkpoint_error_envelope_record(solution.error)},
      {"checkpoint_identity", solution.checkpoint_identity}};
}

json::object checkpoint_regular_tail_model_record(
    const RegularTaylorTailModel& model) {
  json::array q_coefficients;
  q_coefficients.reserve(model.q_coefficients.size());
  for (const auto& coefficient : model.q_coefficients)
    q_coefficients.push_back(checkpoint_ball_record(coefficient));
  json::array n_coefficients;
  n_coefficients.reserve(model.n_coefficients.size());
  for (const auto& matrix : model.n_coefficients) {
    json::array encoded_matrix;
    encoded_matrix.reserve(matrix.size());
    for (const auto& coefficient : matrix)
      encoded_matrix.push_back(checkpoint_ball_record(coefficient));
    n_coefficients.push_back(std::move(encoded_matrix));
  }
  json::array n_row_sum_upper;
  n_row_sum_upper.reserve(model.n_row_sum_upper.size());
  for (const auto& magnitude : model.n_row_sum_upper)
    n_row_sum_upper.emplace_back(magnitude.dump_exact());
  json::array initial_row_upper;
  initial_row_upper.reserve(model.initial_row_upper.size());
  for (const auto& magnitude : model.initial_row_upper)
    initial_row_upper.emplace_back(magnitude.dump_exact());
  json::array prescriptions;
  prescriptions.reserve(model.prescriptions.size());
  for (const auto& prescription : model.prescriptions)
    prescriptions.push_back(json::object{
        {"factor_exact", prescription.factor_exact},
        {"sign", prescription.sign},
        {"multiplicity", prescription.multiplicity},
        {"leading_coefficient_sign",
         prescription.leading_coefficient_sign}});
  return json::object{
      {"schema", "diffexp2-regular-taylor-tail-model-v1"},
      {"epsilon", json::object{{"min", model.epsilon.min_power},
                                {"max", model.epsilon.complete_max}}},
      {"dimension", model.dimension},
      {"taylor_complete_max", model.taylor_complete_max},
      {"q_coefficients", std::move(q_coefficients)},
      {"n_coefficients", std::move(n_coefficients)},
      {"n_row_sum_upper_exact", std::move(n_row_sum_upper)},
      {"initial_row_upper_exact", std::move(initial_row_upper)},
      {"chart", json::object{
           {"center_exact", model.chart.center_exact},
           {"scale_exact", model.chart.scale_exact},
           {"radius_exact_ball", checkpoint_ball_record(model.chart.radius)},
           {"infinite_radius", model.chart.infinite_radius}}},
      {"prescriptions", std::move(prescriptions)},
      {"operator_identity", model.operator_identity},
      {"local_checkpoint_identity", model.local_checkpoint_identity},
      {"provenance", model.provenance}};
}

RegularTaylorTailModel parse_checkpoint_regular_tail_model(
    const json::value& raw) {
  const auto& object = as_object(raw, "checkpoint regular tail model");
  require_exact_keys(
      object,
      {"schema", "epsilon", "dimension", "taylor_complete_max",
       "q_coefficients", "n_coefficients", "n_row_sum_upper_exact",
       "initial_row_upper_exact", "chart", "prescriptions",
       "operator_identity", "local_checkpoint_identity", "provenance"},
      "checkpoint regular tail model");
  if (required_string(object, "schema") !=
      "diffexp2-regular-taylor-tail-model-v1")
    throw std::invalid_argument(
        "unsupported checkpoint regular tail-model schema");
  RegularTaylorTailModel model;
  const auto& epsilon = as_object(
      object.at("epsilon"), "checkpoint tail epsilon window");
  require_exact_keys(epsilon, {"min", "max"},
                     "checkpoint tail epsilon window");
  model.epsilon = {
      as_i32(epsilon.at("min"), "checkpoint tail epsilon minimum"),
      as_i32(epsilon.at("max"), "checkpoint tail epsilon maximum")};
  (void)model.epsilon.width();
  model.dimension = as_u32(
      object.at("dimension"), "checkpoint tail dimension");
  if (model.dimension == 0)
    throw std::invalid_argument("checkpoint tail dimension is zero");
  model.taylor_complete_max = as_u32(
      object.at("taylor_complete_max"),
      "checkpoint tail Taylor complete maximum");
  for (const auto& coefficient : as_array(
           object.at("q_coefficients"),
           "checkpoint tail q coefficients"))
    model.q_coefficients.push_back(parse_checkpoint_ball(
        coefficient, "checkpoint tail q coefficient"));
  for (const auto& raw_matrix : as_array(
           object.at("n_coefficients"),
           "checkpoint tail N coefficients")) {
    std::vector<ComplexBall> matrix;
    for (const auto& coefficient : as_array(
             raw_matrix, "checkpoint tail N matrix"))
      matrix.push_back(parse_checkpoint_ball(
          coefficient, "checkpoint tail N coefficient"));
    model.n_coefficients.push_back(std::move(matrix));
  }
  for (const auto& magnitude : as_array(
           object.at("n_row_sum_upper_exact"),
           "checkpoint tail N norms"))
    model.n_row_sum_upper.push_back(parse_checkpoint_magnitude(
        magnitude, "checkpoint tail N norm"));
  for (const auto& magnitude : as_array(
           object.at("initial_row_upper_exact"),
           "checkpoint tail initial magnitudes"))
    model.initial_row_upper.push_back(parse_checkpoint_magnitude(
        magnitude, "checkpoint tail initial magnitude"));
  const auto& chart = as_object(
      object.at("chart"), "checkpoint tail chart");
  require_exact_keys(
      chart,
      {"center_exact", "scale_exact", "radius_exact_ball",
       "infinite_radius"},
      "checkpoint tail chart");
  model.chart.center_exact = required_string(chart, "center_exact");
  model.chart.scale_exact = required_string(chart, "scale_exact");
  model.chart.radius = parse_checkpoint_ball(
      chart.at("radius_exact_ball"), "checkpoint tail chart radius");
  if (!chart.at("infinite_radius").is_bool())
    throw std::invalid_argument(
        "checkpoint tail infinite-radius flag must be Boolean");
  model.chart.infinite_radius = chart.at("infinite_radius").as_bool();
  for (const auto& raw_prescription : as_array(
           object.at("prescriptions"), "checkpoint tail prescriptions")) {
    const auto& prescription = as_object(
        raw_prescription, "checkpoint tail prescription");
    require_exact_keys(
        prescription,
        {"factor_exact", "sign", "multiplicity",
         "leading_coefficient_sign"},
        "checkpoint tail prescription");
    const auto sign = as_i32(
        prescription.at("sign"), "checkpoint tail prescription sign");
    const auto multiplicity = as_u32(
        prescription.at("multiplicity"),
        "checkpoint tail prescription multiplicity");
    const auto leading = as_i32(
        prescription.at("leading_coefficient_sign"),
        "checkpoint tail leading-coefficient sign");
    if ((sign != -1 && sign != 1) || multiplicity == 0 ||
        (leading != -1 && leading != 1))
      throw std::invalid_argument(
          "checkpoint tail prescription is malformed");
    model.prescriptions.push_back(Prescription{
        required_string(prescription, "factor_exact"), sign,
        multiplicity, leading});
  }
  model.operator_identity = required_string(object, "operator_identity");
  model.local_checkpoint_identity = required_string(
      object, "local_checkpoint_identity");
  model.provenance = required_string(object, "provenance");
  return model;
}

template <typename Scalar>
json::array checkpoint_pseudo_hits_record(
    const std::vector<PseudoHit<Scalar>>& hits) {
  json::array output;
  output.reserve(hits.size());
  for (const auto& hit : hits) {
    json::array columns;
    columns.reserve(hit.columns.size());
    for (const auto column : hit.columns) columns.emplace_back(column);
    json::array frames;
    frames.reserve(hit.gamma_frames.size());
    for (const auto& frame : hit.gamma_frames) {
      json::array coefficients;
      coefficients.reserve(frame.size());
      for (const auto& coefficient : frame)
        coefficients.push_back(checkpoint_scalar_record<Scalar>(coefficient));
      frames.push_back(std::move(coefficients));
    }
    json::array validity;
    validity.reserve(hit.gamma_validity.size());
    for (const auto value : hit.gamma_validity)
      validity.push_back(encode_validity(value));
    output.push_back(json::object{
        {"n", hit.n}, {"columns", std::move(columns)},
        {"delta_b", checkpoint_scalar_record<Scalar>(hit.delta_b)},
        {"gamma_frames", std::move(frames)},
        {"gamma_validity", std::move(validity)}});
  }
  return output;
}

template <typename Scalar>
struct NativeLocalRun {
  LocalSolution<Scalar> solution;
  std::vector<PseudoHit<Scalar>> pseudo_hits;
  NativeLocalDiagnostics diagnostics;
  RegularTaylorTailModelResult tail_model = unavailable_tail_model(
      "tail model was not prepared for this native local run");
};

std::uint64_t scoped_handle_id(const std::string& handle,
                               std::string_view prefix,
                               const char* label);

struct SCCColumnProvenance {
  std::string scc_handle;
  std::string scc_exact_identity;
  std::uint32_t seed_block = 0;
  std::uint32_t basis_index = 0;
  std::string exact_column_identity;

  json::object encode() const {
    return json::object{{"scc", scc_handle},
                        {"scc_exact_identity", scc_exact_identity},
                        {"seed_block", seed_block},
                        {"basis_index", basis_index},
                        {"exact_column_identity", exact_column_identity}};
  }
};

struct RationalShadowColumnWitness {
  std::shared_ptr<const LocalSolution<Rational>> solution;
  std::string rational_shadow_identity;
  std::string source_column_identity;
};

EpsilonLatticeSaturationResult<Rational>
rational_shadow_formal_saturation(
    const std::vector<std::shared_ptr<const RationalShadowColumnWitness>>&
        witnesses,
    EpsilonWindow window, const std::string& context);

struct OwnerBoundResidualBinding {
  std::string kind;
  std::string capability;
  std::string operator_identity;
  std::string source_identity;
  std::string local_checkpoint_identity;
  json::object analytic_metadata;
  std::string owner_signature_identity;
  std::string physical_payload_identity;
  std::string provenance_identity;

  json::object encode() const {
    return json::object{
        {"kind", kind},
        {"capability", capability},
        {"operator_identity", operator_identity},
        {"source_identity", source_identity},
        {"local_checkpoint_identity", local_checkpoint_identity},
        {"analytic_metadata", analytic_metadata},
        {"owner_signature_identity", owner_signature_identity},
        {"physical_payload_identity", physical_payload_identity},
        {"provenance_identity", provenance_identity},
        {"equation", "q(t,eps) theta(f) = C(t,eps) f"},
        {"basis", "physical-original-master"},
        {"source_kind", "homogeneous-zero"},
        {"operator_payload_owner",
         "retained-immutable-physical-equation-owner"},
        {"supported_scopes",
         json::array{"stored_truncation", "full_local_solution-inconclusive"}}};
  }
};

template <typename Scalar>
OwnerBoundResidualBinding make_owner_bound_residual_binding(
    const LocalSolution<Scalar>& solution,
    const PreparedPhysicalClearedODE<Scalar>& equation,
    const PhysicalEquationOwnerBase& owner,
    const std::string& expected_operator_identity) {
  static_assert(std::is_same_v<Scalar, Rational> ||
                    std::is_same_v<Scalar, ComplexBall>,
                "owner-bound physical residual bindings support Rational or Acb "
                "locals only");
  physical_ode_detail::validate_ode(equation);
  validate_local_solution(solution, false);
  if (equation.dimension != solution.dimension ||
      owner.equation_operator_identity() != expected_operator_identity ||
      owner.physical_payload_identity() != equation.payload_identity ||
      owner.owner_signature_identity() != equation.owner_signature_identity)
    throw std::invalid_argument(
        "physical residual owner, equation payload, and local provenance disagree");
  auto analytic_metadata =
      checkpoint_local_analytic_metadata_record(solution);
  const json::object source_provenance{
      {"schema", "diffexp2-retained-homogeneous-zero-source-v1"},
      {"operator_identity", expected_operator_identity},
      {"owner_signature_identity", equation.owner_signature_identity},
      {"physical_payload_identity", equation.payload_identity},
      {"local_checkpoint_identity", solution.checkpoint_identity},
      {"analytic_metadata", analytic_metadata}};
  const auto source_identity = json::serialize(
      canonical_json_value(source_provenance));
  json::object provenance{
      {"schema", "diffexp2-owner-bound-physical-residual-binding-v2"},
      {"kind", "physical-homogeneous-cleared-ode"},
      {"capability", kOwnerBoundPhysicalResidualCapability},
      {"operator_identity", expected_operator_identity},
      {"source_identity", source_identity},
      {"local_checkpoint_identity", solution.checkpoint_identity},
      {"analytic_metadata", analytic_metadata},
      {"owner_kind", owner.equation_owner_kind()},
      {"owner_handle", owner.equation_owner_handle()},
      {"owner_signature_identity", equation.owner_signature_identity},
      {"physical_payload_identity", equation.payload_identity},
      {"equation", "q(t,eps) theta(f) = C(t,eps) f"},
      {"basis", "physical-original-master"},
      {"source_kind", "homogeneous-zero"}};
  return OwnerBoundResidualBinding{
      "physical-homogeneous-cleared-ode",
      kOwnerBoundPhysicalResidualCapability,
      expected_operator_identity,
      source_identity,
      solution.checkpoint_identity,
      std::move(analytic_metadata),
      equation.owner_signature_identity,
      equation.payload_identity,
      json::serialize(canonical_json_value(provenance))};
}

class StoredLocalBase {
 public:
  StoredLocalBase(std::string handle, std::string source_chart,
                  std::string source_operator_identity,
                  double create_parse_ms, double create_kernel_ms,
                  std::optional<SCCColumnProvenance> column_provenance =
                      std::nullopt)
      : handle_(std::move(handle)), source_chart_(std::move(source_chart)),
        source_operator_identity_(std::move(source_operator_identity)),
        create_parse_ms_(create_parse_ms),
        create_kernel_ms_(create_kernel_ms),
        column_provenance_(std::move(column_provenance)) {}
  virtual ~StoredLocalBase() = default;

  virtual json::object evaluate(const json::object& request,
                                int output_digits) = 0;
  virtual LocalEvaluation evaluate_retained_point(
      const RealEvaluationPoint& point,
      std::optional<std::int32_t> exact_rim) = 0;
  virtual std::optional<CertifiedLocalEvaluation>
  evaluate_retained_point_with_certified_tail(
      const RealEvaluationPoint& point,
      std::optional<std::int32_t> exact_rim,
      const std::string& witness_radius_exact) = 0;
  virtual json::object certify_residual(const json::object& request,
                                        int output_digits) = 0;
  virtual EndpointLimitResult endpoint_limit(
      const EndpointLimitOptions& options) = 0;
  virtual StoredLineIntegral integrate_planned_line(
      const ExactAffineChart& planning_chart,
      const ChartGeometry& chart, const ComplexBall& scale_numeric,
      const std::string& scale_exact, std::int32_t scale_sign,
      const std::vector<Prescription>& prescriptions,
      const RealEvaluationPoint& local_begin,
      const RealEvaluationPoint& local_end,
      bool reverse_local_orientation, bool certified_zero_length,
      const EpsilonWindow& delivered_epsilon,
      std::optional<std::int32_t> exact_rim,
      bool certify_tail,
      const std::optional<
          StoredLineIntegrationOptions::BoundedDivergentCancellation>&
          divergent_cancellation = std::nullopt) = 0;
  virtual json::object endpoint_metadata() const = 0;
  virtual json::object exact_analytic_metadata() const = 0;
  virtual void require_exact_plan_binding(
      const ChartGeometry& chart,
      const std::vector<Prescription>& prescriptions,
      const std::string& label) const = 0;
  virtual const std::string& checkpoint_identity() const = 0;
  virtual const char* scalar_domain() const = 0;
  virtual json::object summary() const = 0;
  virtual json::object stats_json() const = 0;
  virtual StoredLocalStats stats() const = 0;
  virtual json::object checkpoint_record() const = 0;
  virtual const std::optional<json::object>& retained_derivation() const = 0;
  virtual std::shared_ptr<void> retained_derivation_owner() const = 0;
  virtual json::object seal_plan_match_lineage() = 0;
  virtual bool has_sealed_plan_match_lineage() const = 0;
  virtual std::shared_ptr<PhysicalEquationOwnerBase>
      retained_equation_owner() const = 0;
  virtual std::shared_ptr<const RationalShadowColumnWitness>
      rational_shadow_witness() const = 0;

  const std::string& handle() const { return handle_; }
  const std::string& source_chart() const { return source_chart_; }
  const std::string& source_operator_identity() const {
    return source_operator_identity_;
  }
  const std::optional<SCCColumnProvenance>& column_provenance() const {
    return column_provenance_;
  }

 protected:
  std::string handle_;
  std::string source_chart_;
  std::string source_operator_identity_;
  double create_parse_ms_ = 0.0;
  double create_kernel_ms_ = 0.0;
  std::optional<SCCColumnProvenance> column_provenance_;
};

template <typename Scalar>
class StoredLocal final : public StoredLocalBase {
 public:
  StoredLocal(std::string handle, std::string source_chart,
              std::string source_operator_identity,
              LocalSolution<Scalar>&& solution, slong precision_bits,
              std::vector<PseudoHit<Scalar>>&& pseudo_hits,
              NativeLocalDiagnostics diagnostics,
              std::optional<SCCColumnProvenance> column_provenance =
                  std::nullopt,
              std::optional<json::object> retained_derivation =
                  std::nullopt,
              std::shared_ptr<void> retained_owner = nullptr,
              RegularTaylorTailModelResult tail_model =
                  unavailable_tail_model(
                      "tail model is unavailable for this retained local"),
              std::optional<TailModelCheckpointMarker>
                  tail_checkpoint_marker = std::nullopt,
              bool serialize_tail_checkpoint_fields = true,
              bool serialize_derivation_checkpoint_fields = true,
              std::shared_ptr<PhysicalEquationOwnerBase> equation_owner =
                  nullptr,
              std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
                  physical_equation = nullptr,
              std::string residual_unavailable_reason = {},
              std::shared_ptr<const RationalShadowColumnWitness>
                  rational_shadow_witness = nullptr,
              bool sealed_plan_match_lineage = false)
      : StoredLocalBase(std::move(handle), std::move(source_chart),
                        std::move(source_operator_identity),
                        diagnostics.parse_ms, diagnostics.kernel_ms,
                        std::move(column_provenance)),
        solution_(std::move(solution)), precision_bits_(precision_bits),
        pseudo_hits_(std::move(pseudo_hits)), top_valid_(diagnostics.top_valid),
        retained_derivation_(std::move(retained_derivation)),
        retained_owner_(std::move(retained_owner)),
        tail_model_(std::move(tail_model)),
        tail_checkpoint_marker_(std::move(tail_checkpoint_marker)),
        serialize_tail_checkpoint_fields_(serialize_tail_checkpoint_fields),
        serialize_derivation_checkpoint_fields_(
            serialize_derivation_checkpoint_fields),
        equation_owner_(std::move(equation_owner)),
        physical_equation_(std::move(physical_equation)),
        rational_shadow_witness_(std::move(rational_shadow_witness)),
        sealed_plan_match_lineage_(sealed_plan_match_lineage) {
    validate_local_solution(solution_, false);
    if (retained_derivation_.has_value()) {
      const auto* raw_identity =
          retained_derivation_->if_contains("provenance_identity");
      if (raw_identity == nullptr || !raw_identity->is_string() ||
          raw_identity->as_string().empty())
        throw std::invalid_argument(
            "retained local derivation has no provenance identity");
      const auto& identity = raw_identity->as_string();
      retained_provenance_identity_bytes_ = identity.size();
      retained_provenance_fingerprint_ = public_provenance_fingerprint(
          std::string_view(identity.data(), identity.size()));
    }
    if constexpr (std::is_same_v<Scalar, Rational> ||
                  std::is_same_v<Scalar, ComplexBall>) {
      const bool homogeneous_match_derivation =
          retained_derivation_.has_value() &&
          (required_string(*retained_derivation_, "schema") ==
               "diffexp2-retained-plan-match-local-materialization-v1" ||
           required_string(*retained_derivation_, "schema") ==
               "diffexp2-retained-plan-match-local-materialization-v2" ||
           is_retained_plan_value_handoff_schema(
               required_string(*retained_derivation_, "schema")));
      if (sealed_plan_match_lineage_ &&
          (!homogeneous_match_derivation || retained_owner_ != nullptr))
        throw std::invalid_argument(
            "sealed plan-match lineage requires one owner-free plan-match derivation");
      if (homogeneous_match_derivation && column_provenance_.has_value())
        throw std::invalid_argument(
            "a plan-match combination cannot retain canonical SCC-column provenance");
      if (equation_owner_ != nullptr && physical_equation_ != nullptr) {
        if (retained_derivation_.has_value() &&
            !homogeneous_match_derivation)
          throw std::invalid_argument(
              "only a homogeneous plan-match derivation may retain a physical equation owner");
        const auto owner_kind =
            std::string(equation_owner_->equation_owner_kind());
        if (owner_kind == "prepared-chart") {
          if (column_provenance_.has_value())
            throw std::invalid_argument(
                "SCC column local cannot acquire a primitive chart equation owner");
        } else if (owner_kind == "composite-scc") {
          const bool completed_column =
              column_provenance_.has_value() &&
              column_provenance_->scc_handle ==
                  equation_owner_->equation_owner_handle() &&
              column_provenance_->scc_exact_identity ==
                  equation_owner_->equation_operator_identity();
          if (!completed_column && !homogeneous_match_derivation)
            throw std::invalid_argument(
                "SCC-owned local is neither a completed full column nor a homogeneous plan-match combination");
        } else {
          throw std::invalid_argument(
              "local received an unsupported physical equation owner kind");
        }
        if (equation_owner_->equation_owner_handle() != source_chart_ ||
            equation_owner_->equation_operator_identity() !=
                source_operator_identity_ ||
            std::string(equation_owner_->equation_scalar_domain()) !=
                scalar_domain() ||
            equation_owner_->physical_payload_record() !=
                physical_equation_->exact_payload_record)
          throw std::invalid_argument(
              "physical equation owner disagrees with its local source or q/C provenance");
        residual_binding_ = make_owner_bound_residual_binding(
            solution_, *physical_equation_, *equation_owner_,
            source_operator_identity_);
        residual_binding_detail_ =
            "owner-bound physical homogeneous q/C payload is available";
      } else if (equation_owner_ != nullptr || physical_equation_ != nullptr) {
        throw std::invalid_argument(
            "primitive residual ownership requires both owner and physical payload");
      } else if (!residual_unavailable_reason.empty()) {
        residual_binding_detail_ = std::move(residual_unavailable_reason);
      } else if (retained_derivation_.has_value()) {
        residual_binding_detail_ =
            "owner-bound residual is unsupported for derived locals until their equation/source provenance is retained";
      } else if (column_provenance_.has_value()) {
        residual_binding_detail_ =
            "owner-bound residual is unsupported for composite SCC locals until one full physical q/C certificate is retained";
      } else {
        residual_binding_detail_ =
            "owner-bound residual is unsupported: the prepared chart has no physical q/C payload";
      }
    } else {
      residual_binding_detail_ =
          "owner-bound residual rejects unresolved symbolic coefficients";
    }
  }

  json::object evaluate(const json::object& request,
                        int output_digits) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.evaluate rejects unresolved symbolic coefficients; solve a "
          "numerically specialized chart before native evaluation");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto point = parse_local_evaluation_point(request);
      auto options = parse_local_evaluation_options(request, true);
      const auto tail_witness = parse_certified_tail_witness(request);
      if (tail_witness.has_value())
        options.compute_tail_estimate = false;
      const auto started = std::chrono::steady_clock::now();
      LocalEvaluation result;
      std::optional<RegularTaylorPointTailCertificate> tail_certificate;
      TailMajorantStatus requested_tail_status = tail_model_.status;
      std::string requested_tail_detail = tail_model_.detail;
      if (tail_witness.has_value() && tail_model_.model.has_value()) {
        auto certified = evaluate_local_solution_with_certified_tail(
            solution_, *tail_model_.model, point, *tail_witness, options);
        result = std::move(certified.evaluation);
        requested_tail_status = certified.tail.status;
        requested_tail_detail = certified.tail.detail;
        tail_certificate = std::move(certified.tail);
      } else {
        result = evaluate_local_solution(solution_, point, options);
      }
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      evaluations_.fetch_add(1);
      if (tail_witness.has_value()) {
        tail_certificate_requests_.fetch_add(1);
        switch (requested_tail_status) {
          case TailMajorantStatus::Certified:
            tail_certificate_certified_.fetch_add(1);
            break;
          case TailMajorantStatus::Inconclusive:
            tail_certificate_inconclusive_.fetch_add(1);
            break;
          case TailMajorantStatus::Unsupported:
            tail_certificate_unsupported_.fetch_add(1);
            break;
        }
      }
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        evaluate_ms_ += elapsed;
      }
      json::object response{
          {"point_exact", point.exact_coordinate},
          {"imaginary_sign", result.imaginary_sign.has_value()
               ? json::value(*result.imaginary_sign) : json::value(nullptr)},
          {"arithmetic_enclosed", result.arithmetic_enclosed},
          {"elapsed_ms", elapsed},
          {"value", encode_epsilon_vector(result.value, output_digits)},
          {"theta", encode_epsilon_vector(result.theta_value, output_digits)}};
      if (tail_witness.has_value()) {
        if (tail_certificate.has_value()) {
          response["tail_certificate"] = encode_point_tail_certificate(
              *tail_certificate, *tail_witness);
        } else {
          response["tail_certificate"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"requested", true},
              {"status", tail_majorant_status_name(requested_tail_status)},
              {"model_status", tail_majorant_status_name(tail_model_.status)},
              {"witness_radius_exact", *tail_witness},
              {"detail", requested_tail_detail}};
        }
      }
      return response;
    }
  }

  LocalEvaluation evaluate_retained_point(
      const RealEvaluationPoint& point,
      std::optional<std::int32_t> exact_rim) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "retained-point evaluation rejects unresolved symbolic coefficients");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      EvaluationOptions options;
      options.imaginary_sign = exact_rim;
      options.compute_tail_estimate = false;
      const auto started = std::chrono::steady_clock::now();
      auto result = evaluate_local_solution(solution_, point, options);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      evaluations_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        evaluate_ms_ += elapsed;
      }
      return result;
    }
  }

  std::optional<CertifiedLocalEvaluation>
  evaluate_retained_point_with_certified_tail(
      const RealEvaluationPoint& point,
      std::optional<std::int32_t> exact_rim,
      const std::string& witness_radius_exact) override {
    if constexpr (!std::is_same_v<Scalar, ComplexBall>) {
      return std::nullopt;
    } else {
      if (!tail_model_.model.has_value() ||
          tail_model_.status != TailMajorantStatus::Certified)
        return std::nullopt;
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      EvaluationOptions options;
      options.imaginary_sign = exact_rim;
      options.compute_tail_estimate = false;
      const auto started = std::chrono::steady_clock::now();
      auto result = evaluate_local_solution_with_certified_tail(
          solution_, *tail_model_.model, point, witness_radius_exact,
          options);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      evaluations_.fetch_add(1);
      tail_certificate_requests_.fetch_add(1);
      switch (result.tail.status) {
        case TailMajorantStatus::Certified:
          tail_certificate_certified_.fetch_add(1);
          break;
        case TailMajorantStatus::Inconclusive:
          tail_certificate_inconclusive_.fetch_add(1);
          break;
        case TailMajorantStatus::Unsupported:
          tail_certificate_unsupported_.fetch_add(1);
          break;
      }
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        evaluate_ms_ += elapsed;
      }
      return result;
    }
  }

  json::object certify_residual(const json::object& request,
                                int output_digits) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.certify_residual rejects unresolved symbolic coefficients; "
          "solve a numerically specialized chart first");
    } else {
      if (request.if_contains("output_digits") != nullptr)
        require_exact_keys(
            request,
            {"schema", "op", "session", "local", "point", "options",
             "relative_tolerance", "scope", "include_residual",
             "operator_identity", "source_identity",
             "checkpoint_identity", "analytic_metadata",
             "owner_signature_identity", "physical_payload_identity",
             "provenance_identity", "output_digits"},
            "owner-bound residual request");
      else
        require_exact_keys(
            request,
            {"schema", "op", "session", "local", "point", "options",
             "relative_tolerance", "scope", "include_residual",
             "operator_identity", "source_identity",
             "checkpoint_identity", "analytic_metadata",
             "owner_signature_identity", "physical_payload_identity",
             "provenance_identity"},
            "owner-bound residual request");
      if (!residual_binding_.has_value() || physical_equation_ == nullptr ||
          equation_owner_ == nullptr)
        throw std::domain_error(residual_binding_detail_);
      const auto& binding = *residual_binding_;
      if (required_string(request, "operator_identity") !=
              binding.operator_identity ||
          required_string(request, "source_identity") !=
              binding.source_identity ||
          required_string(request, "checkpoint_identity") !=
              binding.local_checkpoint_identity ||
          required_string(request, "owner_signature_identity") !=
              binding.owner_signature_identity ||
          required_string(request, "physical_payload_identity") !=
              binding.physical_payload_identity ||
          required_string(request, "provenance_identity") !=
              binding.provenance_identity ||
          request.at("analytic_metadata") != binding.analytic_metadata)
        throw std::invalid_argument(
            "owner-bound residual expectations differ from the retained "
            "operator/source/checkpoint/provenance/analytic metadata");

      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto point = parse_local_evaluation_point(request);
      const auto options = parse_local_evaluation_options(request, false);
      if (options.compute_tail_estimate)
        throw std::invalid_argument(
            "stored-truncation residual certification rejects advisory tail estimates");
      if (!solution_.error.empty())
        throw std::domain_error(
            "owner-bound residual certification does not ignore a retained local error envelope");
      if (point.sign < 0 && options.imaginary_sign.has_value()) {
        const auto retained_sign = derive_chart_imaginary_sign(solution_);
        if (retained_sign.has_value() &&
            *retained_sign != *options.imaginary_sign)
          throw std::invalid_argument(
              "explicit residual rim disagrees with the retained analytic prescription");
      }
      const auto tolerance_text = required_string(
          request, "relative_tolerance");
      const auto tolerance = Magnitude::decimal(tolerance_text);
      const auto scope_text = required_string(request, "scope");
      const auto scope = scope_text == "stored_truncation"
          ? ResidualScope::StoredTruncation
          : scope_text == "full_local_solution"
              ? ResidualScope::FullLocalSolution
              : throw std::invalid_argument(
                    "residual scope must be stored_truncation or full_local_solution");
      if (!request.at("include_residual").is_bool())
        throw std::invalid_argument(
            "owner-bound residual include_residual must be Boolean");
      const auto include_residual = request.at("include_residual").as_bool();

      const auto started = std::chrono::steady_clock::now();
      const auto evaluation = evaluate_local_solution(solution_, point, options);
      auto certificate = certify_physical_cleared_ode_residual(
          *physical_equation_, evaluation, point, tolerance, scope);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      residual_certifications_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        residual_certify_ms_ += elapsed;
      }

      const char* verdict = certificate.verdict == ResidualVerdict::Pass
          ? "pass"
          : certificate.verdict == ResidualVerdict::Fail
              ? "fail" : "inconclusive";
      json::object result{
          {"point_exact", point.exact_coordinate},
          {"imaginary_sign", evaluation.imaginary_sign.has_value()
               ? json::value(*evaluation.imaginary_sign)
               : json::value(nullptr)},
          {"arithmetic_enclosed", evaluation.arithmetic_enclosed},
          {"scope", scope_text}, {"verdict", verdict},
          {"detail", certificate.detail},
          {"relative_tolerance", tolerance_text},
          {"binding_kind", binding.kind},
          {"capability", binding.capability},
          {"checkpoint_identity", binding.local_checkpoint_identity},
          {"operator_identity", binding.operator_identity},
          {"source_identity", binding.source_identity},
          {"provenance_identity", binding.provenance_identity},
          {"owner_signature_identity", binding.owner_signature_identity},
          {"physical_payload_identity", binding.physical_payload_identity},
          {"analytic_metadata", binding.analytic_metadata},
          {"equation", "q(t,eps) theta(f) = C(t,eps) f"},
          {"basis", "physical-original-master"},
          {"source_kind", "homogeneous-zero"},
          {"epsilon_min", certificate.residual.epsilon.min_power},
          {"epsilon_max", certificate.residual.epsilon.complete_max},
          {"dimension", certificate.residual.dimension},
          {"residual_upper_approx",
           encode_magnitude_diagnostics(certificate.residual_upper)},
          {"scale_lower_approx",
           encode_magnitude_diagnostics(certificate.scale_lower)},
          {"relative_upper_approx",
           encode_magnitude_diagnostics(certificate.relative_upper)},
          {"bound_encoding", "approximate-double-diagnostics"},
          {"json_coefficients", include_residual
               ? certificate.residual.coefficients.size() : 0},
          {"elapsed_ms", elapsed}};
      if (include_residual)
        result["residual"] = encode_epsilon_vector(
            certificate.residual, output_digits);
      return result;
    }
  }

  EndpointLimitResult endpoint_limit(
      const EndpointLimitOptions& options) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "local.endpoint_limit rejects unresolved symbolic coefficients; "
          "solve an explicitly specialized chart before endpoint evaluation");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto started = std::chrono::steady_clock::now();
      auto result = endpoint_sector_limit(solution_, options);
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      endpoint_limits_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        endpoint_limit_ms_ += elapsed;
      }
      return result;
    }
  }

  StoredLineIntegral integrate_planned_line(
      const ExactAffineChart& planning_chart,
      const ChartGeometry& chart, const ComplexBall& scale_numeric,
      const std::string& scale_exact, std::int32_t scale_sign,
      const std::vector<Prescription>& prescriptions,
      const RealEvaluationPoint& local_begin,
      const RealEvaluationPoint& local_end,
      bool reverse_local_orientation, bool certified_zero_length,
      const EpsilonWindow& delivered_epsilon,
      std::optional<std::int32_t> exact_rim,
      bool certify_tail,
      const std::optional<
          StoredLineIntegrationOptions::BoundedDivergentCancellation>&
          divergent_cancellation) override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "native planned line integration rejects unresolved symbolic "
          "coefficients; specialize the retained local first");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      if (solution_.chart.center_exact != chart.center_exact ||
          solution_.chart.scale_exact != chart.scale_exact ||
          solution_.chart.infinite_radius != chart.infinite_radius ||
          (!solution_.chart.infinite_radius &&
           !acb_equal(solution_.chart.radius.raw(), chart.radius.raw())))
        throw std::invalid_argument(
            "retained local chart geometry differs from its native tile plan");
      const auto same_prescription = [](const Prescription& left,
                                        const Prescription& right) {
        return left.factor_exact == right.factor_exact &&
               left.sign == right.sign &&
               left.multiplicity == right.multiplicity &&
               left.leading_coefficient_sign ==
                   right.leading_coefficient_sign;
      };
      if (solution_.prescriptions.size() != prescriptions.size() ||
          !std::equal(solution_.prescriptions.begin(),
                      solution_.prescriptions.end(), prescriptions.begin(),
                      same_prescription))
        throw std::invalid_argument(
            "retained local prescriptions differ from its native tile plan");

      StoredLineIntegrationOptions options;
      options.delivered_epsilon = delivered_epsilon;
      options.imaginary_sign = exact_rim;
      options.certified_chart_scale_sign = scale_sign;
      options.divergent_cancellation = divergent_cancellation;
      const bool local_zero =
          local_begin.exact_coordinate == local_end.exact_coordinate;
      if (local_zero != certified_zero_length)
        throw std::invalid_argument(
            "native planned line zero-length authorization contradicts its exact local endpoints");
      if (certified_zero_length) {
        std::optional<std::int32_t> chart_sign;
        try {
          chart_sign = derive_chart_imaginary_sign(solution_, scale_sign);
        } catch (const std::domain_error& error) {
          throw NativeIntegrationError(
              NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
              std::string("invalid prepared chart branch prescription: ") +
                  error.what());
        }
        if (chart_sign.has_value() && exact_rim.has_value() &&
            *chart_sign != *exact_rim)
          throw NativeIntegrationError(
              NativeIntegrationErrorCode::MissingBranchPrescription, "E3",
              "explicit line branch sign conflicts with the prepared chart");
        auto zero = certified_zero_physical_line(
            delivered_epsilon, solution_.dimension,
            exact_rim.has_value() ? exact_rim : chart_sign,
            local_begin.sign == 0, certify_tail);
        line_integrations_.fetch_add(1);
        return zero;
      }
      const auto& primitive_begin =
          reverse_local_orientation ? local_end : local_begin;
      const auto& primitive_end =
          reverse_local_orientation ? local_begin : local_end;
      const auto started = std::chrono::steady_clock::now();
      const auto& begin_point = primitive_begin;
      const auto& end_point = primitive_end;
      StoredLineIntegral result;
      if (certify_tail && !begin_point.certified_algebraic &&
          !end_point.certified_algebraic &&
          (tail_model_.model.has_value() ||
                           (rational_row_tail_model_.has_value() &&
                            rational_row_tail_model_->model.has_value()))) {
        const Rational begin_exact(begin_point.exact_coordinate);
        const Rational end_exact(end_point.exact_coordinate);
        const auto begin_modulus = begin_exact.sign() < 0
            ? -begin_exact : begin_exact;
        const auto end_modulus = end_exact.sign() < 0
            ? -end_exact : end_exact;
        const auto outer = begin_modulus < end_modulus
            ? end_modulus : begin_modulus;
        const auto& chart_radius = planning_chart.radius;
        if (!(outer < chart_radius))
          throw std::invalid_argument(
              "planned line endpoint is not strictly inside its exact chart radius");
        const auto witness =
            (outer + chart_radius) / Rational(2);
        auto certified = tail_model_.model.has_value()
            ? integrate_regular_local_line_with_certified_tail(
                  solution_, *tail_model_.model, begin_point, end_point,
                  options, witness.str())
            : integrate_rational_row_local_line_with_certified_tail(
                  solution_, *rational_row_tail_model_->model,
                  begin_point, end_point, options, witness.str());
        result = std::move(certified.integral);
      } else {
        result = integrate_stored_local_line(
            solution_, begin_point, end_point, options);
        if (certify_tail) {
          const auto requested_status = rational_row_tail_model_.has_value()
              ? rational_row_tail_model_->status : tail_model_.status;
          const auto& requested_detail = rational_row_tail_model_.has_value()
              ? rational_row_tail_model_->detail : tail_model_.detail;
          result.diagnostics.tail_certificate_requested = true;
          result.diagnostics.tail_certificate_status =
              tail_majorant_status_name(requested_status);
          result.value.error.provenance =
              std::string(tail_majorant_status_name(requested_status)) +
              ": " + requested_detail +
              "; returned value remains stored Taylor truncation only";
          result.diagnostics.detail = result.value.error.provenance;
        }
      }
      // integrate_stored_local_line is deliberately a local-coordinate
      // primitive.  A retained physical tile has dx = scale dt, so apply the
      // exact affine Jacobian before publishing the physical line result.
      const auto jacobian = reverse_local_orientation
          ? -scale_numeric : scale_numeric;
      const auto oriented_jacobian_exact = reverse_local_orientation
          ? "-(" + scale_exact + ")" : scale_exact;
      for (auto& coefficient : result.value.coefficients)
        coefficient *= jacobian;
      if (!result.value.error.empty()) {
        const auto jacobian_upper = Magnitude::upper_abs(jacobian);
        for (auto& bound : result.value.error.absolute)
          bound = bound * jacobian_upper;
        result.value.error.provenance +=
            "; physical_jacobian_exact=" + oriented_jacobian_exact;
      }
      const auto elapsed = std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - started).count();
      line_integrations_.fetch_add(1);
      {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        line_integration_ms_ += elapsed;
      }
      return result;
    }
  }

  json::object endpoint_metadata() const override { return metadata_json(); }

  json::object exact_analytic_metadata() const override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "exact checkpoint analytic metadata does not serialize symbolic coefficients");
    } else {
      return checkpoint_local_analytic_metadata_record(solution_);
    }
  }

  void require_exact_plan_binding(
      const ChartGeometry& chart,
      const std::vector<Prescription>& prescriptions,
      const std::string& label) const override {
    AcbPrecisionLease lease(precision_bits_);
    ComplexBall::set_precision(precision_bits_);
    if (solution_.chart.center_exact != chart.center_exact ||
        solution_.chart.scale_exact != chart.scale_exact ||
        solution_.chart.infinite_radius != chart.infinite_radius ||
        (!solution_.chart.infinite_radius &&
         !acb_equal(solution_.chart.radius.raw(), chart.radius.raw())))
      throw std::invalid_argument(
          label +
          " retained local geometry differs from its exact tile-plan chart");
    const auto same_prescription = [](const Prescription& left,
                                      const Prescription& right) {
      return left.factor_exact == right.factor_exact &&
             left.sign == right.sign &&
             left.multiplicity == right.multiplicity &&
             left.leading_coefficient_sign ==
                 right.leading_coefficient_sign;
    };
    if (solution_.prescriptions.size() != prescriptions.size() ||
        !std::equal(solution_.prescriptions.begin(),
                    solution_.prescriptions.end(), prescriptions.begin(),
                    same_prescription))
      throw std::invalid_argument(
          label +
          " retained local prescriptions differ from its exact tile-plan chart");
  }

  const std::string& checkpoint_identity() const override {
    return solution_.checkpoint_identity;
  }

  const char* scalar_domain() const override {
    if constexpr (std::is_same_v<Scalar, Rational>) return "rational";
    if constexpr (std::is_same_v<Scalar, ComplexBall>) return "acb";
    return "symbolic";
  }

  json::object summary() const override {
    json::object result{
        {"local", handle_}, {"chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"dimension", solution_.dimension},
        {"epsilon_min", solution_.epsilon.min_power},
        {"epsilon_max", solution_.epsilon.complete_max},
        {"taylor_complete_max", solution_.taylor_complete_max},
        {"sectors", solution_.sectors.size()},
        {"coefficient_count", coefficient_count()},
        {"pseudo_hit_count", pseudo_hits_.size()},
        {"top_valid", encode_validity(top_valid_)},
        {"checkpoint_identity", solution_.checkpoint_identity},
        {"tail_majorant", encode_tail_model_status(tail_model_)},
        {"residual_binding", residual_binding_.has_value()
             ? json::value(json::object{
                   {"status", "available"},
                   {"binding", residual_binding_->encode()}})
             : json::value(json::object{
                   {"status", "unsupported"},
                   {"kind", "none"},
                   {"capability", kOwnerBoundPhysicalResidualCapability},
                   {"reason", residual_binding_detail_}})},
        {"metadata", metadata_json()},
        {"create_parse_ms", create_parse_ms_},
        {"create_kernel_ms", create_kernel_ms_}};
    if (column_provenance_.has_value())
      result["column_provenance"] = column_provenance_->encode();
    if (retained_derivation_.has_value()) {
      const auto& derivation = *retained_derivation_;
      if (!retained_provenance_fingerprint_.has_value())
        throw std::logic_error(
            "retained local lost its public provenance fingerprint");
      result["retained_derivation"] = json::object{
          {"schema", derivation.at("schema")},
          {"capability", derivation.at("capability")},
          {"provenance", json::object{
               {"algorithm", "fnv1a64-v1"},
               {"fingerprint", *retained_provenance_fingerprint_},
               {"identity_bytes", retained_provenance_identity_bytes_}}},
          {"checkpoint_identity", solution_.checkpoint_identity},
          {"ownership", retained_owner_ != nullptr
               ? "strong"
               : sealed_plan_match_lineage_
                     ? "sealed-plan-match-lineage"
                     : "unowned"}};
      result["strong_derivation_ownership"] =
          retained_owner_ != nullptr;
    }
    if (rational_row_tail_model_.has_value()) {
      json::object projected_tail{
          {"capability",
           "retained-rational-row-analytic-line-tail-v1"},
          {"status", tail_majorant_status_name(
                         rational_row_tail_model_->status)},
          {"attached", rational_row_tail_model_->model.has_value()},
          {"detail", rational_row_tail_model_->detail},
          // The derived model is scratch state.  Completed certified line
          // envelopes are checkpointed exactly; a fresh transport
          // contraction after restore re-derives this model from its restored
          // source and caller-supplied prepared row.
          {"checkpoint_serialized", false}};
      if (rational_row_tail_model_->model.has_value()) {
        projected_tail["source_checkpoint_identity"] =
            rational_row_tail_model_->model->source_checkpoint_identity;
        projected_tail["local_checkpoint_identity"] =
            rational_row_tail_model_->model->local_checkpoint_identity;
        projected_tail["row_exact_identity"] =
            rational_row_tail_model_->model->row_exact_identity;
      }
      result["rational_row_line_tail_majorant"] =
          std::move(projected_tail);
    }
    return result;
  }

  json::object stats_json() const override {
    auto out = summary();
    const auto current = stats();
    out["evaluations"] = current.evaluations;
    out["residual_certifications"] = current.residual_certifications;
    out["endpoint_limits"] = current.endpoint_limits;
    out["line_integrations"] = current.line_integrations;
    out["evaluate_ms"] = current.evaluate_ms;
    out["residual_certify_ms"] = current.residual_certify_ms;
    out["endpoint_limit_ms"] = current.endpoint_limit_ms;
    out["line_integration_ms"] = current.line_integration_ms;
    out["tail_certificate_requests"] = tail_certificate_requests_.load();
    out["tail_certificate_certified"] =
        tail_certificate_certified_.load();
    out["tail_certificate_inconclusive"] =
        tail_certificate_inconclusive_.load();
    out["tail_certificate_unsupported"] =
        tail_certificate_unsupported_.load();
    return out;
  }

  StoredLocalStats stats() const override {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return {evaluations_.load(), residual_certifications_.load(),
            endpoint_limits_.load(), line_integrations_.load(), evaluate_ms_,
            residual_certify_ms_, endpoint_limit_ms_, line_integration_ms_,
            create_parse_ms_, create_kernel_ms_, coefficient_count(),
            tail_certificate_requests_.load(),
            tail_certificate_certified_.load(),
            tail_certificate_inconclusive_.load(),
            tail_certificate_unsupported_.load()};
  }

  json::object checkpoint_record() const override {
    if constexpr (std::is_same_v<Scalar, SymbolicRational>) {
      throw std::domain_error(
          "native checkpoint does not serialize symbolic-coefficient local state");
    } else {
      AcbPrecisionLease lease(precision_bits_);
      ComplexBall::set_precision(precision_bits_);
      const auto current = stats();
      if (retained_derivation_.has_value()) {
        const auto schema = required_string(
            *retained_derivation_, "schema");
        if (schema !=
                "diffexp2-retained-plan-match-local-materialization-v1" &&
            schema !=
                "diffexp2-retained-plan-match-local-materialization-v2" &&
            !is_retained_plan_value_handoff_schema(schema) &&
            schema !=
                "diffexp2-retained-rational-row-local-application-v1")
          throw std::domain_error(
              "native checkpoint does not serialize this retained local derivation kind");
        if (retained_owner_ == nullptr && !sealed_plan_match_lineage_)
          throw std::logic_error(
              "derived local lost its strong derivation owner before checkpointing");
      } else if (retained_owner_ != nullptr) {
        throw std::logic_error(
            "primitive local unexpectedly retains a derivation owner before checkpointing");
      }
      json::value equation_owner_restore = nullptr;
      if (equation_owner_ != nullptr) {
        if (physical_equation_ == nullptr ||
            equation_owner_->physical_payload_identity() !=
                physical_equation_->payload_identity ||
            equation_owner_->owner_signature_identity() !=
                physical_equation_->owner_signature_identity ||
            equation_owner_->physical_payload_record() !=
                physical_equation_->exact_payload_record)
          throw std::logic_error(
              "primitive local physical equation owner changed before checkpointing");
        equation_owner_restore = json::object{
            {"owner_kind", equation_owner_->equation_owner_kind()},
            {"owner_handle", equation_owner_->equation_owner_handle()},
            {"operator_identity",
             equation_owner_->equation_operator_identity()},
            {"owner_signature_identity",
             equation_owner_->owner_signature_identity()},
            {"physical_payload_identity",
             equation_owner_->physical_payload_identity()}};
      } else if (physical_equation_ != nullptr) {
        throw std::logic_error(
            "primitive local lost its physical equation owner before checkpointing");
      }
      json::value owner_lineage = nullptr;
      if (retained_derivation_.has_value()) {
        const auto& derivation = *retained_derivation_;
        const auto derivation_schema = required_string(derivation, "schema");
        if (derivation_schema ==
            "diffexp2-retained-plan-match-local-materialization-v1") {
          owner_lineage = json::object{
              {"match", derivation.at("source_match")},
              {"match_checkpoint_identity",
               derivation.at("source_match_checkpoint_identity")},
              {"match_provenance_identity",
               derivation.at("source_match_provenance_identity")},
              {"planned_hop_provenance_identity",
               derivation.at("planned_hop_provenance_identity")},
              {"derivation_provenance_identity",
               derivation.at("provenance_identity")}};
        } else if (derivation_schema ==
            "diffexp2-retained-plan-match-local-materialization-v2") {
          owner_lineage = json::object{
              {"match", derivation.at("source_match")},
              {"match_checkpoint_identity",
               derivation.at("source_match_checkpoint_identity")},
              {"tile_plan", derivation.at("tile_plan")},
              {"tile_plan_checkpoint_identity",
               derivation.at("tile_plan_checkpoint_identity")},
              {"arm", derivation.at("arm")},
              {"match_index", derivation.at("match")}};
        } else if (is_retained_plan_value_handoff_schema(
                       derivation_schema)) {
          owner_lineage = json::object{
              {"tile_plan", derivation.at("tile_plan")},
              {"tile_plan_checkpoint_identity",
               derivation.at("tile_plan_checkpoint_identity")},
              {"arm", derivation.at("arm")},
              {"match_index", derivation.at("match")},
              {"incoming_checkpoint_identity",
               as_object(derivation.at("incoming"),
                         "checkpoint value-handoff incoming")
                   .at("checkpoint_identity")},
              {"handoff_provenance_identity",
               derivation.at("provenance_identity")}};
        } else {
          owner_lineage = rational_row_owner_lineage();
        }
      }
      json::object runtime{
          {"evaluations", current.evaluations},
          {"residual_certifications", current.residual_certifications},
          {"endpoint_limits", current.endpoint_limits},
          {"line_integrations", current.line_integrations},
          {"evaluate_ms", current.evaluate_ms},
          {"residual_certify_ms", current.residual_certify_ms},
          {"endpoint_limit_ms", current.endpoint_limit_ms},
          {"line_integration_ms", current.line_integration_ms},
          {"coefficient_count", current.coefficient_count}};
      if (serialize_tail_checkpoint_fields_) {
        runtime["tail_certificate_requests"] =
            current.tail_certificate_requests;
        runtime["tail_certificate_certified"] =
            current.tail_certificate_certified;
        runtime["tail_certificate_inconclusive"] =
            current.tail_certificate_inconclusive;
        runtime["tail_certificate_unsupported"] =
            current.tail_certificate_unsupported;
      }
      json::object record{
        {"schema", sealed_plan_match_lineage_
             ? "diffexp2-retained-local-v5"
             : "diffexp2-retained-local-v4"},
        {"handle", handle_},
        {"source_chart", source_chart_},
        {"source_operator_identity", source_operator_identity_},
        {"scalar_domain", scalar_domain()},
        {"precision_bits", precision_bits_},
        {"solution", checkpoint_local_solution_record(solution_)},
        {"pseudo_hits", checkpoint_pseudo_hits_record(pseudo_hits_)},
        {"diagnostics",
         json::object{{"top_valid", encode_validity(top_valid_)},
                      {"create_parse_ms", create_parse_ms_},
                      {"create_kernel_ms", create_kernel_ms_}}},
        {"runtime_stats", std::move(runtime)},
        {"column_provenance", column_provenance_.has_value()
             ? json::value(column_provenance_->encode())
             : json::value(nullptr)},
        {"equation_owner_restore", std::move(equation_owner_restore)}};
      if (sealed_plan_match_lineage_)
        record["derivation_owner_restore"] =
            retained_derivation_.has_value() &&
                    is_retained_plan_value_handoff_schema(
                        required_string(*retained_derivation_, "schema"))
                ? "sealed-plan-value-handoff-lineage"
                : "sealed-plan-match-lineage";
      if (serialize_derivation_checkpoint_fields_) {
        record["retained_derivation"] = retained_derivation_.has_value()
            ? json::value(*retained_derivation_) : json::value(nullptr);
        record["retained_owner_lineage"] = std::move(owner_lineage);
      }
      if (serialize_tail_checkpoint_fields_) {
        if (tail_model_.model.has_value()) {
          if (tail_model_.status != TailMajorantStatus::Certified ||
              tail_checkpoint_marker_.has_value())
            throw std::logic_error(
                "attached regular tail model has inconsistent checkpoint status");
          tail_majorant_detail::validate_restored_regular_taylor_tail_model(
              *tail_model_.model, solution_, source_operator_identity_);
          record["tail_model_restore"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"serialized", true},
              {"status", "certified"},
              {"attached_before_save", true},
              {"model", checkpoint_regular_tail_model_record(
                   *tail_model_.model)}};
        } else {
          record["tail_model_restore"] = json::object{
              {"capability", kRegularTailMajorantCapability},
              {"serialized", false},
              {"status", tail_checkpoint_marker_.has_value()
                   ? tail_checkpoint_marker_->saved_status
                   : tail_majorant_status_name(tail_model_.status)},
              {"attached_before_save", tail_checkpoint_marker_.has_value()
                   ? tail_checkpoint_marker_->attached_before_save
                   : false}};
        }
      }
      if (residual_binding_.has_value()) {
        if (equation_owner_ == nullptr || physical_equation_ == nullptr)
          throw std::logic_error(
              "owner-bound residual lost its retained physical q/C owner before checkpointing");
        const auto recomputed = make_owner_bound_residual_binding(
            solution_, *physical_equation_, *equation_owner_,
            source_operator_identity_);
        if (recomputed.encode() != residual_binding_->encode())
          throw std::logic_error(
              "owner-bound residual binding changed before checkpointing");
        record["residual_operator_restore"] = json::object{
            {"capability", kOwnerBoundPhysicalResidualCapability},
            {"kind", residual_binding_->kind},
            {"serialized", true},
            {"status", "available"},
            {"operator_payload_owner",
             "retained-immutable-physical-equation-owner"},
            {"binding", residual_binding_->encode()},
            {"reason", nullptr}};
      } else {
        record["residual_operator_restore"] = json::object{
            {"capability", kOwnerBoundPhysicalResidualCapability},
            {"kind", "none"},
            {"serialized", false},
            {"status", "unsupported"},
            {"operator_payload_owner", nullptr},
            {"binding", nullptr},
            {"reason", residual_binding_detail_}};
      }
      return record;
    }
  }

  void restore_runtime_stats(const StoredLocalStats& state) {
    if (state.coefficient_count != coefficient_count())
      throw std::invalid_argument(
          "checkpoint local coefficient count does not match its tensor");
    evaluations_.store(state.evaluations);
    residual_certifications_.store(state.residual_certifications);
    endpoint_limits_.store(state.endpoint_limits);
    line_integrations_.store(state.line_integrations);
    tail_certificate_requests_.store(state.tail_certificate_requests);
    tail_certificate_certified_.store(state.tail_certificate_certified);
    tail_certificate_inconclusive_.store(
        state.tail_certificate_inconclusive);
    tail_certificate_unsupported_.store(state.tail_certificate_unsupported);
    std::lock_guard<std::mutex> lock(stats_mutex_);
    evaluate_ms_ = state.evaluate_ms;
    residual_certify_ms_ = state.residual_certify_ms;
    endpoint_limit_ms_ = state.endpoint_limit_ms;
    line_integration_ms_ = state.line_integration_ms;
  }

  const LocalSolution<Scalar>& solution() const { return solution_; }
  const RegularTaylorTailModelResult& tail_model() const {
    return tail_model_;
  }
  void attach_rational_row_line_tail_model(
      RationalRowLineTailModelResult<Scalar> model) {
    if (!retained_derivation_.has_value() ||
        required_string(*retained_derivation_, "schema") !=
            "diffexp2-retained-rational-row-local-application-v1" ||
        retained_owner_ == nullptr || rational_row_tail_model_.has_value())
      throw std::invalid_argument(
          "rational-row line-tail model requires one newly derived retained scalar local");
    if (model.model.has_value() &&
        model.model->local_checkpoint_identity !=
            solution_.checkpoint_identity)
      throw std::invalid_argument(
          "rational-row line-tail model checkpoint binding differs from its scalar local");
    rational_row_tail_model_ = std::move(model);
  }
  const std::optional<OwnerBoundResidualBinding>& residual_binding() const {
    return residual_binding_;
  }
  void restore_residual_binding_detail(std::string detail) {
    if (residual_binding_.has_value() || detail.empty())
      throw std::invalid_argument(
          "checkpoint unsupported residual reason is incompatible with an available binding");
    residual_binding_detail_ = std::move(detail);
  }
  std::int32_t top_valid() const { return top_valid_; }
  const std::optional<json::object>& retained_derivation() const override {
    return retained_derivation_;
  }
  std::shared_ptr<void> retained_derivation_owner() const override {
    return retained_owner_;
  }
  json::object seal_plan_match_lineage() override {
    if (sealed_plan_match_lineage_)
      throw std::logic_error(
          "plan-match local lineage is already sealed");
    if (!retained_derivation_.has_value() || retained_owner_ == nullptr)
      throw std::invalid_argument(
          "only a strongly owned plan-match local lineage can be sealed");
    const auto& derivation = *retained_derivation_;
    const auto derivation_schema = required_string(derivation, "schema");
    const bool value_handoff = is_retained_plan_value_handoff_schema(
        derivation_schema);
    if (!value_handoff && derivation_schema !=
            "diffexp2-retained-plan-match-local-materialization-v1" &&
        derivation_schema !=
            "diffexp2-retained-plan-match-local-materialization-v2")
      throw std::invalid_argument(
          "only a plan-match materialization lineage can be sealed");
    const auto& output = as_object(
        derivation.at("output"), "sealed plan-match output");
    if (required_string(output, "checkpoint_identity") !=
            solution_.checkpoint_identity ||
        required_string(output, "chart") != source_chart_ ||
        required_string(output, "source_operator_identity") !=
            source_operator_identity_)
      throw std::logic_error(
          "plan-match derivation output changed before sealing");
    if (equation_owner_ == nullptr || physical_equation_ == nullptr ||
        equation_owner_->equation_owner_handle() != source_chart_ ||
        equation_owner_->equation_operator_identity() !=
            source_operator_identity_ ||
        equation_owner_->physical_payload_identity() !=
            physical_equation_->payload_identity ||
        equation_owner_->owner_signature_identity() !=
            physical_equation_->owner_signature_identity ||
        equation_owner_->physical_payload_record() !=
            physical_equation_->exact_payload_record)
      throw std::logic_error(
          "plan-match local lost its exact physical equation owner before sealing");
    json::object sealed;
    if (value_handoff) {
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
          "sealed plan-value handoff derivation");
      if (required_string(derivation, "capability") !=
              (derivation_schema ==
                       "diffexp2-retained-plan-value-handoff-v2"
                   ? "retained-native-regular-value-handoff-v2"
                   : "retained-native-regular-value-handoff-v1") ||
          required_string(derivation, "scope") !=
              "single-regular-to-regular-transport-hop" ||
          required_string(derivation, "coefficient_transport") !=
              "native-retained-only" ||
          !derivation.at("whole_arm_complete").is_bool() ||
          derivation.at("whole_arm_complete").as_bool() ||
          required_string(derivation,
              "equation_owner_signature_identity") !=
              equation_owner_->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner_->physical_payload_identity())
        throw std::logic_error(
            "plan-value handoff changed its scope or equation-owner binding before sealing");
      auto identity_input = derivation;
      const auto identity = required_string(
          identity_input, "provenance_identity");
      identity_input.erase("provenance_identity");
      if (json::serialize(canonical_json_value(identity_input)) != identity)
        throw std::logic_error(
            "plan-value handoff derivation identity changed before sealing");
      const auto& incoming = as_object(
          derivation.at("incoming"), "sealed plan-value incoming");
      sealed = json::object{
          {"schema", "diffexp2-sealed-plan-value-handoff-lineage-v1"},
          {"local", handle_},
          {"local_checkpoint_identity", solution_.checkpoint_identity},
          {"source_operator_identity", source_operator_identity_},
          {"tile_plan", derivation.at("tile_plan")},
          {"tile_plan_checkpoint_identity",
           derivation.at("tile_plan_checkpoint_identity")},
          {"arm", derivation.at("arm")},
          {"match_index", derivation.at("match")},
          {"incoming_checkpoint_identity",
           incoming.at("checkpoint_identity")},
          {"handoff_provenance_identity", identity},
          {"equation_owner_signature_identity",
           equation_owner_->owner_signature_identity()},
          {"equation_payload_identity",
           equation_owner_->physical_payload_identity()}};
    } else if (derivation_schema ==
        "diffexp2-retained-plan-match-local-materialization-v2") {
      if (required_string(derivation,
              "equation_owner_signature_identity") !=
              equation_owner_->owner_signature_identity() ||
          required_string(derivation, "equation_payload_identity") !=
              equation_owner_->physical_payload_identity())
        throw std::logic_error(
            "compact plan-match derivation changed its equation-owner binding before sealing");
      sealed = json::object{
        {"schema", "diffexp2-sealed-plan-match-local-lineage-v2"},
        {"local", handle_},
        {"local_checkpoint_identity", solution_.checkpoint_identity},
        {"source_operator_identity", source_operator_identity_},
        {"match", derivation.at("source_match")},
        {"match_checkpoint_identity",
         derivation.at("source_match_checkpoint_identity")},
        {"tile_plan", derivation.at("tile_plan")},
        {"tile_plan_checkpoint_identity",
         derivation.at("tile_plan_checkpoint_identity")},
        {"arm", derivation.at("arm")},
        {"match_index", derivation.at("match")},
        {"incoming_checkpoint_identity",
         as_object(derivation.at("incoming"),
                   "compact derivation incoming").at("checkpoint_identity")},
        {"equation_owner_signature_identity",
         equation_owner_->owner_signature_identity()},
        {"equation_payload_identity",
         equation_owner_->physical_payload_identity()}};
    } else sealed = json::object{
        {"schema", "diffexp2-sealed-plan-match-local-lineage-v1"},
        {"local", handle_},
        {"local_checkpoint_identity", solution_.checkpoint_identity},
        {"source_operator_identity", source_operator_identity_},
        {"match", derivation.at("source_match")},
        {"match_checkpoint_identity",
         derivation.at("source_match_checkpoint_identity")},
        {"match_provenance_identity",
         derivation.at("source_match_provenance_identity")},
        {"planned_hop_provenance_identity",
         derivation.at("planned_hop_provenance_identity")},
        {"derivation_provenance_identity",
         derivation.at("provenance_identity")},
        {"equation_owner_signature_identity",
         equation_owner_->owner_signature_identity()},
        {"equation_payload_identity",
         equation_owner_->physical_payload_identity()}};
    retained_owner_.reset();
    sealed_plan_match_lineage_ = true;
    return sealed;
  }
  bool has_sealed_plan_match_lineage() const override {
    return sealed_plan_match_lineage_;
  }
  std::shared_ptr<PhysicalEquationOwnerBase> retained_equation_owner()
      const override {
    return equation_owner_;
  }
  std::shared_ptr<const RationalShadowColumnWitness>
  rational_shadow_witness() const override {
    return rational_shadow_witness_;
  }
  const std::vector<PseudoHit<Scalar>>& pseudo_hits() const {
    return pseudo_hits_;
  }

 private:
  json::object metadata_json() const {
    json::array prescriptions;
    for (const auto& prescription : solution_.prescriptions) {
      prescriptions.push_back(json::object{
          {"factor_exact", prescription.factor_exact},
          {"sign", prescription.sign},
          {"multiplicity", prescription.multiplicity},
          {"leading_coefficient_sign",
           prescription.leading_coefficient_sign}});
    }
    json::array sectors;
    sectors.reserve(solution_.sectors.size());
    for (const auto& sector : solution_.sectors)
      sectors.push_back(json::object{
          {"a", encode_exact_descriptor(sector.a)},
          {"b", encode_exact_descriptor(sector.b)},
          {"log_power", sector.log_power}});
    const auto& sector = solution_.sectors.front();
    return json::object{
        {"chart", json::object{
            {"center_exact", solution_.chart.center_exact},
            {"scale_exact", solution_.chart.scale_exact},
            {"infinite_radius", solution_.chart.infinite_radius},
            {"radius_ball", encode_scalar(solution_.chart.radius, 30)}}},
        {"tag", json::object{{"a", encode_exact_descriptor(sector.a)},
                             {"b", encode_exact_descriptor(sector.b)}}},
        {"sectors", std::move(sectors)},
        {"prescriptions", std::move(prescriptions)}};
  }

  std::size_t coefficient_count() const {
    std::size_t count = 0;
    for (const auto& sector : solution_.sectors)
      count += sector.coefficients.size();
    return count;
  }

  json::object rational_row_owner_lineage() const {
    if (!retained_derivation_.has_value() || retained_owner_ == nullptr)
      throw std::logic_error(
          "rational-row local lost its derivation or source owner");
    const auto& derivation = *retained_derivation_;
    require_exact_keys(
        derivation,
        {"schema", "capability", "source", "row", "output",
         "analytic_prescriptions", "coefficient_transport",
         "provenance_identity"},
        "retained rational-row local derivation");
    if (required_string(derivation, "schema") !=
            "diffexp2-retained-rational-row-local-application-v1" ||
        required_string(derivation, "capability") !=
            "retained-native-rational-row-local-application-v1" ||
        required_string(derivation, "analytic_prescriptions") !=
            "preserved-exactly" ||
        required_string(derivation, "coefficient_transport") !=
            "native-retained-only")
      throw std::invalid_argument(
          "retained rational-row derivation changes its certified scope");

    auto erased_source =
        std::static_pointer_cast<StoredLocalBase>(retained_owner_);
    auto source = std::dynamic_pointer_cast<StoredLocal<Scalar>>(
        erased_source);
    if (!source || std::string(source->scalar_domain()) != scalar_domain())
      throw std::invalid_argument(
          "retained rational-row derivation source domain is inconsistent");
    const auto& source_record = as_object(
        derivation.at("source"), "retained rational-row source");
    require_exact_keys(
        source_record,
        {"local", "chart", "source_operator_identity",
         "checkpoint_identity", "dimension", "epsilon",
         "taylor_complete_max"},
        "retained rational-row source");
    const auto& source_epsilon = as_object(
        source_record.at("epsilon"), "retained rational-row source epsilon");
    require_exact_keys(source_epsilon, {"min", "max"},
                       "retained rational-row source epsilon");
    if (required_string(source_record, "local") != source->handle() ||
        required_string(source_record, "chart") != source->source_chart() ||
        required_string(source_record, "source_operator_identity") !=
            source->source_operator_identity() ||
        required_string(source_record, "checkpoint_identity") !=
            source->checkpoint_identity() ||
        as_u32(source_record.at("dimension"),
               "rational-row source dimension") !=
            source->solution().dimension ||
        as_i32(source_epsilon.at("min"),
               "rational-row source epsilon minimum") !=
            source->solution().epsilon.min_power ||
        as_i32(source_epsilon.at("max"),
               "rational-row source epsilon maximum") !=
            source->solution().epsilon.complete_max ||
        as_u32(source_record.at("taylor_complete_max"),
               "rational-row source Taylor maximum") !=
            source->solution().taylor_complete_max)
      throw std::invalid_argument(
          "retained rational-row source provenance disagrees with its strong owner");

    const auto& row = as_object(
        derivation.at("row"), "retained rational-row identity");
    require_exact_keys(
        row, {"exact_identity", "columns", "active_entries",
              "structurally_zero"},
        "retained rational-row identity");
    const auto row_identity = required_string(row, "exact_identity");
    if (as_u32(row.at("columns"), "retained rational-row columns") !=
        source->solution().dimension)
      throw std::invalid_argument(
          "retained rational-row column count differs from its source");
    if (!row.at("structurally_zero").is_bool())
      throw std::invalid_argument(
          "retained rational-row zero fact must be Boolean");
    const auto& active = as_array(
        row.at("active_entries"), "retained rational-row active entries");
    if (row.at("structurally_zero").as_bool() != active.empty())
      throw std::invalid_argument(
          "retained rational-row zero fact disagrees with its active entries");
    std::optional<std::uint32_t> previous_column;
    for (const auto& raw_entry : active) {
      const auto& entry = as_object(
          raw_entry, "retained rational-row active entry");
      require_exact_keys(
          entry, {"column", "epsilon_shift", "center_pole_order",
                  "exact_identity"},
          "retained rational-row active entry");
      const auto column = as_u32(
          entry.at("column"), "retained rational-row active column");
      if (column >= source->solution().dimension ||
          (previous_column.has_value() && *previous_column >= column))
        throw std::invalid_argument(
            "retained rational-row active columns are not canonical");
      previous_column = column;
      (void)as_i32(entry.at("epsilon_shift"),
                   "retained rational-row epsilon shift");
      (void)as_u32(entry.at("center_pole_order"),
                   "retained rational-row pole order");
      (void)required_string(entry, "exact_identity");
    }

    const auto& output = as_object(
        derivation.at("output"), "retained rational-row output");
    require_exact_keys(
        output, {"checkpoint_identity", "dimension", "epsilon",
                 "taylor_complete_max"},
        "retained rational-row output");
    const auto& output_epsilon = as_object(
        output.at("epsilon"), "retained rational-row output epsilon");
    require_exact_keys(output_epsilon, {"min", "max"},
                       "retained rational-row output epsilon");
    if (required_string(output, "checkpoint_identity") !=
            solution_.checkpoint_identity ||
        as_u32(output.at("dimension"),
               "rational-row output dimension") != solution_.dimension ||
        as_i32(output_epsilon.at("min"),
               "rational-row output epsilon minimum") !=
            solution_.epsilon.min_power ||
        as_i32(output_epsilon.at("max"),
               "rational-row output epsilon maximum") !=
            solution_.epsilon.complete_max ||
        as_u32(output.at("taylor_complete_max"),
               "rational-row output Taylor maximum") !=
            solution_.taylor_complete_max)
      throw std::invalid_argument(
          "retained rational-row output provenance disagrees with its tensor");
    if (source_chart_ != source->source_chart() ||
        !local_algebra_detail::same_chart(
            solution_.chart, source->solution().chart) ||
        !local_algebra_detail::same_prescriptions(
            solution_.prescriptions, source->solution().prescriptions))
      throw std::invalid_argument(
          "retained rational-row output left its source analytic chart");

    auto identity_input = derivation;
    const auto derivation_identity = required_string(
        derivation, "provenance_identity");
    identity_input.erase("provenance_identity");
    if (json::serialize(canonical_json_value(identity_input)) !=
        derivation_identity)
      throw std::invalid_argument(
          "retained rational-row derivation identity is inconsistent");
    const json::object operator_provenance{
        {"schema", "diffexp2-rational-row-derived-operator-v1"},
        {"source_operator_identity", source->source_operator_identity()},
        {"row_exact_identity", row_identity},
        {"provenance_identity", derivation_identity}};
    if (json::serialize(canonical_json_value(operator_provenance)) !=
        source_operator_identity_)
      throw std::invalid_argument(
          "retained rational-row derived operator identity is inconsistent");

    return json::object{
        {"source_local", source->handle()},
        {"source_chart", source->source_chart()},
        {"source_operator_identity", source->source_operator_identity()},
        {"source_checkpoint_identity", source->checkpoint_identity()},
        {"row_exact_identity", row_identity},
        {"derivation_provenance_identity", derivation_identity},
        {"derived_operator_identity", source_operator_identity_}};
  }

  LocalSolution<Scalar> solution_;
  slong precision_bits_ = 256;
  std::vector<PseudoHit<Scalar>> pseudo_hits_;
  std::int32_t top_valid_ = kCompleteInfinity;
  std::optional<json::object> retained_derivation_;
  std::shared_ptr<void> retained_owner_;
  std::optional<std::string> retained_provenance_fingerprint_;
  std::size_t retained_provenance_identity_bytes_ = 0;
  RegularTaylorTailModelResult tail_model_ = unavailable_tail_model(
      "tail model is unavailable for this retained local");
  std::optional<RationalRowLineTailModelResult<Scalar>>
      rational_row_tail_model_;
  std::shared_ptr<PhysicalEquationOwnerBase> equation_owner_;
  std::shared_ptr<const PreparedPhysicalClearedODE<Scalar>>
      physical_equation_;
  std::shared_ptr<const RationalShadowColumnWitness>
      rational_shadow_witness_;
  bool sealed_plan_match_lineage_ = false;
  std::optional<OwnerBoundResidualBinding> residual_binding_;
  std::string residual_binding_detail_ =
      "owner-bound residual payload was not prepared";
  std::optional<TailModelCheckpointMarker> tail_checkpoint_marker_;
  bool serialize_tail_checkpoint_fields_ = true;
  bool serialize_derivation_checkpoint_fields_ = true;
  std::atomic<std::uint64_t> evaluations_{0};
  std::atomic<std::uint64_t> residual_certifications_{0};
  std::atomic<std::uint64_t> endpoint_limits_{0};
  std::atomic<std::uint64_t> line_integrations_{0};
  std::atomic<std::uint64_t> tail_certificate_requests_{0};
  std::atomic<std::uint64_t> tail_certificate_certified_{0};
  std::atomic<std::uint64_t> tail_certificate_inconclusive_{0};
  std::atomic<std::uint64_t> tail_certificate_unsupported_{0};
  mutable std::mutex stats_mutex_;
  double evaluate_ms_ = 0.0;
  double residual_certify_ms_ = 0.0;
  double endpoint_limit_ms_ = 0.0;
  double line_integration_ms_ = 0.0;
};
