// Shared decoder for every native rational multiplier payload.  Producers
// may supply an already-expanded finite rectangle, a compact list of exact
// numerator/denominator polynomials, or both.  Keeping this policy in one
// place prevents public rows, SCC couplings, gauge frames, and spectral
// frames from independently reintroducing eager O(epsilon*taylor) JSON.
template <typename Scalar>
PreparedRationalTaylorMultiplier<Scalar>
parse_prepared_rational_taylor_multiplier(
    const json::object& raw, std::size_t epsilon_width,
    std::size_t taylor_width, bool exact_rectangle, const char* label) {
  const bool has_kernels = raw.if_contains("kernels") != nullptr;
  const bool has_analytic =
      raw.if_contains("analytic_coefficients") != nullptr;
  if (has_kernels && has_analytic)
    require_exact_keys(raw,
        {"epsilon_shift", "center_pole_order", "kernels",
         "exact_identity", "proven_zero", "analytic_coefficients"}, label);
  else if (has_kernels)
    require_exact_keys(raw,
        {"epsilon_shift", "center_pole_order", "kernels",
         "exact_identity", "proven_zero"}, label);
  else if (has_analytic)
    require_exact_keys(raw,
        {"epsilon_shift", "center_pole_order", "exact_identity",
         "proven_zero", "analytic_coefficients"}, label);
  else
    throw std::invalid_argument(
        std::string(label) +
        " needs Taylor kernels or a compact analytic rational source");
  if (!raw.at("proven_zero").is_bool())
    throw std::invalid_argument(
        std::string(label) + " structural-zero fact must be boolean");

  PreparedRationalTaylorMultiplier<Scalar> multiplier;
  multiplier.epsilon_shift = as_i32(
      raw.at("epsilon_shift"), "rational multiplier epsilon shift");
  multiplier.center_pole_order = as_u32(
      raw.at("center_pole_order"), "rational multiplier center-pole order");
  multiplier.exact_identity = required_string(raw, "exact_identity");
  multiplier.proven_zero = raw.at("proven_zero").as_bool();
  if (multiplier.exact_identity.empty())
    throw std::invalid_argument(
        std::string(label) + " exact identity must be nonempty");

  // A structural zero carries no numerical rectangle.  SCC coupling
  // records retain those entries because their zero proof is checked
  // against the parent matrix, whereas sparse public rows omit them.  Make
  // that one representation explicit here instead of forcing each caller
  // to manufacture epsilon_width*taylor_width zero coefficients.
  if (multiplier.proven_zero) {
    if (has_kernels &&
        !as_array(raw.at("kernels"),
                  "zero rational multiplier epsilon kernels").empty())
      throw std::invalid_argument(
          std::string(label) +
          " structural zero must not carry Taylor kernels");
    if (has_analytic &&
        !as_array(raw.at("analytic_coefficients"),
                  "zero rational multiplier analytic coefficients").empty())
      throw std::invalid_argument(
          std::string(label) +
          " structural zero must not carry analytic coefficients");
    return multiplier;
  }

  const auto covers = [exact_rectangle](std::size_t actual,
                                         std::size_t required) {
    return exact_rectangle ? actual == required : actual >= required;
  };
  std::vector<PreparedRationalAnalyticCoefficient<Scalar>> analytic;
  if (has_analytic) {
    const auto& raw_coefficients = as_array(
        raw.at("analytic_coefficients"),
        "rational multiplier analytic coefficients");
    if (!covers(raw_coefficients.size(), epsilon_width))
      throw std::invalid_argument(
          std::string(label) +
          " analytic source does not cover its epsilon rectangle");
    analytic.reserve(raw_coefficients.size());
    for (const auto& raw_coefficient : raw_coefficients) {
      const auto& coefficient = as_object(
          raw_coefficient, "rational multiplier analytic coefficient");
      require_exact_keys(coefficient, {"numerator", "denominator"},
                         "rational multiplier analytic coefficient");
      PreparedRationalAnalyticCoefficient<Scalar> parsed;
      for (const auto& value : as_array(
               coefficient.at("numerator"),
               "rational multiplier analytic numerator"))
        parsed.numerator.push_back(parse_scalar<Scalar>(value));
      for (const auto& value : as_array(
               coefficient.at("denominator"),
               "rational multiplier analytic denominator"))
        parsed.denominator.push_back(parse_scalar<Scalar>(value));
      if (parsed.numerator.empty() || parsed.denominator.empty())
        throw std::invalid_argument(
            std::string(label) +
            " analytic numerator/denominator cannot be empty");
      analytic.push_back(std::move(parsed));
    }
  }

  if (has_kernels) {
    const auto& raw_kernels = as_array(
        raw.at("kernels"), "rational multiplier epsilon kernels");
    if (!covers(raw_kernels.size(), epsilon_width))
      throw std::invalid_argument(
          std::string(label) + " kernels do not cover its epsilon rectangle");
    if (has_analytic && analytic.size() != raw_kernels.size())
      throw std::invalid_argument(
          std::string(label) +
          " analytic coefficient count differs from its kernels");
    multiplier.kernels.reserve(raw_kernels.size());
    for (const auto& raw_kernel : raw_kernels) {
      const auto& coefficients = as_array(
          raw_kernel, "rational multiplier Taylor kernel");
      if (!covers(coefficients.size(), taylor_width))
        throw std::invalid_argument(
            std::string(label) +
            " kernels do not cover its Taylor rectangle");
      std::vector<Scalar> kernel;
      kernel.reserve(coefficients.size());
      for (const auto& coefficient : coefficients)
        kernel.push_back(parse_scalar<Scalar>(coefficient));
      multiplier.kernels.push_back(std::move(kernel));
    }
  } else {
    multiplier.kernels.reserve(analytic.size());
    for (const auto& rational : analytic)
      multiplier.kernels.push_back(
          local_algebra_detail::expand_rational_taylor(
              rational, taylor_width));
  }
  if (has_analytic)
    multiplier.analytic_coefficients = std::move(analytic);
  return multiplier;
}
