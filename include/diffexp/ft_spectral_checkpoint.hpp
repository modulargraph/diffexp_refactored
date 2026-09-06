#pragma once
#include "diffexp/adjoint_checkpoint.hpp"
#include "diffexp/ft_spectral.hpp"

namespace diffexp::ft_spectral_checkpoint {
namespace json = boost::json;
namespace fs = std::filesystem;
struct Storage {
  fs::path directory;
  std::size_t max_bytes = 64 * 1024 * 1024;
  std::size_t *reuse_counter = nullptr;
};
// Bump this revision for changes to sampling, gauges, selection, or tail
// estimates.
inline constexpr const char *algorithm = "DiffExp.FTSpectralCompletedArm/v1";
inline std::string identity(const ExactEpsilonMatrix &matrix,
                            const LaurentRows &initial,
                            const ExactEpsilonMatrix &forcing,
                            const std::vector<Exact> &vertices,
                            const ft_spectral::Options &spectral,
                            const AdjointOptions &native = {}) {
  const json::object payload{
      {"algorithm", algorithm},
      {"adjoint_identity", adjoint_checkpoint::identity(
                               matrix, initial, forcing, vertices, native)},
      {"diagonal_gauge", spectral.diagonal_gauge},
      {"endpoint_clustering", spectral.endpoint_clustering},
      {"conservative", spectral.conservative},
      {"accuracy_goal", spectral.accuracy_goal},
      {"max_nodes", spectral.max_nodes},
      {"max_block_size", spectral.max_block_size},
      {"max_block_nodes", spectral.max_block_nodes},
      {"max_cells", spectral.max_cells},
      {"seconds_budget_ieee_bits",
       std::to_string(std::bit_cast<std::uint64_t>(spectral.seconds_budget))}};
  return artifacts::detail::sha256(artifacts::detail::canonical(payload));
}
namespace detail {
inline unsigned leg_count(const std::vector<Exact> &path) {
  unsigned count = 0;
  for (std::size_t i = 0; i + 1 < path.size(); ++i)
    if (path[i] != path[i + 1])
      ++count;
  return count;
}
inline int expected_low(const LaurentRows &initial,
                        const ExactEpsilonMatrix &forcing,
                        std::size_t epsilon) {
  int low = std::min(0, initial.low);
  for (const auto &row : forcing)
    for (const auto &q : row)
      if (!q.is_zero()) {
        auto value = *exact_epsilon_valuation(q, epsilon);
        if (value < -1000)
          throw std::invalid_argument("FT spectral checkpoint forcing window");
        if (value < low)
          low = static_cast<int>(value);
      }
  return low;
}
inline void validate(const LaurentRows &result, const LaurentRows &initial,
                     int low) {
  if (result.columns() != initial.columns() ||
      result.coefficients.size() != initial.coefficients.size() ||
      result.low != low || result.high != initial.high)
    throw std::invalid_argument(
        "FT spectral checkpoint result shape/window mismatch");
  for (const auto &row : result.coefficients)
    for (const auto &series : row)
      for (const auto &value : series)
        if (!value.is_finite())
          throw std::invalid_argument(
              "FT spectral checkpoint nonfinite result");
}
inline json::array integers(const std::vector<unsigned> &v) {
  json::array a;
  for (auto i : v)
    a.emplace_back(i);
  return a;
}
inline json::object encode(const LaurentRows &rows, const std::string &key,
                           const ft_spectral::Diagnostics &diagnostics) {
  return {
      {"schema", algorithm},
      {"identity", key},
      {"certificate", "uncertified retained numerical checkpoint"},
      {"rows", numerical_rows_io::exact_rows(rows)},
      {"diagnostics",
       json::object{{"legs", diagnostics.legs},
                    {"nodes", integers(diagnostics.nodes)},
                    {"block_sizes", integers(diagnostics.block_sizes)},
                    {"normalized_diagonals", diagnostics.normalized_diagonals},
                    {"clustered_legs", diagnostics.clustered_legs},
                    {"absolute_stability_components",
                     diagnostics.absolute_stability_components}}}};
}
inline unsigned count(const json::value &value) {
  auto n = artifacts::detail::integer(value);
  if (n < 0 ||
      static_cast<std::uint64_t>(n) > std::numeric_limits<unsigned>::max())
    throw std::invalid_argument("FT spectral checkpoint diagnostic count");
  return static_cast<unsigned>(n);
}
inline LaurentRows decode(const json::value &envelope, const std::string &key,
                          const LaurentRows &initial, int low, unsigned legs,
                          const ft_spectral::Options &options,
                          ft_spectral::Diagnostics &diagnostics) {
  const auto &root = envelope.as_object();
  artifacts::detail::keys(root, {"payload", "sha256"});
  const auto &payload = root.at("payload");
  if (artifacts::detail::string(root.at("sha256")) !=
      artifacts::detail::sha256(artifacts::detail::canonical(payload)))
    throw std::invalid_argument("FT spectral checkpoint checksum mismatch");
  const auto &p = payload.as_object();
  artifacts::detail::keys(
      p, {"schema", "identity", "certificate", "rows", "diagnostics"});
  if (p.at("schema") != algorithm ||
      artifacts::detail::string(p.at("identity")) != key ||
      p.at("certificate") != "uncertified retained numerical checkpoint")
    throw std::invalid_argument(
        "FT spectral checkpoint identity/schema mismatch");
  auto rows = numerical_rows_io::read_rows(p.at("rows"));
  validate(rows, initial, low);
  const auto &d = p.at("diagnostics").as_object();
  artifacts::detail::keys(d, {"legs", "nodes", "block_sizes",
                              "normalized_diagonals", "clustered_legs",
                              "absolute_stability_components"});
  ft_spectral::Diagnostics loaded;
  loaded.legs = count(d.at("legs"));
  if (loaded.legs != legs)
    throw std::invalid_argument(
        "FT spectral checkpoint is not a completed arm");
  const auto &nodes = d.at("nodes").as_array();
  const auto &blocks = d.at("block_sizes").as_array();
  if (nodes.size() > 9ULL * legs ||
      blocks.size() > initial.columns() * static_cast<std::size_t>(legs))
    throw std::invalid_argument("FT spectral checkpoint diagnostic shape");
  for (const auto &value : nodes) {
    auto n = count(value);
    if (n < 8 || n > options.max_nodes)
      throw std::invalid_argument("FT spectral checkpoint node count");
    loaded.nodes.push_back(n);
  }
  for (const auto &value : blocks) {
    auto n = count(value);
    if (!n || n > options.max_block_size)
      throw std::invalid_argument("FT spectral checkpoint block size");
    loaded.block_sizes.push_back(n);
  }
  loaded.normalized_diagonals = count(d.at("normalized_diagonals"));
  loaded.clustered_legs = count(d.at("clustered_legs"));
  loaded.absolute_stability_components =
      count(d.at("absolute_stability_components"));
  if (loaded.normalized_diagonals > initial.columns() ||
      loaded.clustered_legs > legs)
    throw std::invalid_argument("FT spectral checkpoint gauge count");
  // No factorization or numerical work occurs on reuse; retained structure is
  // reported.
  diagnostics = std::move(loaded);
  return rows;
}
} // namespace detail
inline std::optional<LaurentRows> try_transport(
    const ExactEpsilonMatrix &matrix, const LaurentRows &initial,
    const ExactEpsilonMatrix &forcing, const std::vector<Exact> &vertices,
    const ft_spectral::Options &spectral, ft_spectral::Diagnostics &diagnostics,
    const AdjointOptions &native = {}, const Storage &storage = {}) {
  if (storage.directory.empty())
    return ft_spectral::try_transport(matrix, initial, forcing, vertices,
                                      spectral, diagnostics);
  adjoint_checkpoint::detail::limits(
      {storage.directory, storage.max_bytes, nullptr});
  if (native.continuation)
    throw std::invalid_argument(
        "FT spectral completed-arm cache cannot use a partial continuation");
  const auto columns = initial.columns();
  if (matrix.size() != columns ||
      forcing.size() != initial.coefficients.size() || vertices.empty())
    throw std::invalid_argument("FT spectral checkpoint input shape");
  for (const auto *m : {&matrix, &forcing})
    for (const auto &row : *m)
      if (row.size() != columns)
        throw std::invalid_argument("FT spectral checkpoint matrix shape");
  auto [xi, ei] = path_epsilon_variables(vertices.front());
  (void)xi;
  const auto low = detail::expected_low(initial, forcing, ei);
  const auto legs = detail::leg_count(vertices);
  const auto key =
      identity(matrix, initial, forcing, vertices, spectral, native);
  const auto file = storage.directory / (key + ".json");
  if (fs::exists(file)) {
    auto rows = detail::decode(
        adjoint_checkpoint::detail::read(file, storage.max_bytes), key, initial,
        low, legs, spectral, diagnostics);
    if (storage.reuse_counter)
      ++*storage.reuse_counter;
    return rows;
  }
  auto rows = ft_spectral::try_transport(matrix, initial, forcing, vertices,
                                         spectral, diagnostics);
  if (!rows)
    return std::nullopt;
  // A mixed local/spectral fallback cannot be published as a spectral arm.
  if (diagnostics.legs != legs)
    throw std::logic_error(
        "FT spectral checkpoint requires every arm leg to be spectral");
  detail::validate(*rows, initial, low);
  const auto payload = detail::encode(*rows, key, diagnostics);
  const auto bytes = artifacts::detail::canonical(json::object{
      {"payload", payload},
      {"sha256",
       artifacts::detail::sha256(artifacts::detail::canonical(payload))}});
  if (bytes.size() > storage.max_bytes)
    throw std::length_error("FT spectral checkpoint exceeds file budget");
  adjoint_checkpoint::detail::publish(file, bytes);
  return rows;
}
} // namespace diffexp::ft_spectral_checkpoint
