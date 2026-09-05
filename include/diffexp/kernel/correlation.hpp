#pragma once

#include <boost/json.hpp>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp::kernel::correlation {

inline constexpr const char* kSchema = "DiffExp3.Correlation/v1";

enum class Model { Exact, SharedLinearSources, IndependentHull };

inline const char* model_name(Model model) {
  switch (model) {
    case Model::Exact: return "Exact";
    case Model::SharedLinearSources: return "SharedLinearSources";
    case Model::IndependentHull: return "IndependentHull";
  }
  throw std::logic_error("unknown correlation model");
}

inline Model parse_model(const std::string& name) {
  if (name == "Exact") return Model::Exact;
  if (name == "SharedLinearSources") return Model::SharedLinearSources;
  if (name == "IndependentHull") return Model::IndependentHull;
  throw std::invalid_argument("unsupported correlation model: " + name);
}

struct Descriptor {
  std::string schema = kSchema;
  Model model = Model::IndependentHull;
  std::string identity;
  std::vector<std::string> source_ids;
  std::vector<std::string> exact_linear_map_ids;
  bool lossy_hull_explicit = true;

  [[nodiscard]] bool valid() const {
    if (schema != kSchema || identity.empty()) return false;
    const auto unique_nonempty = [](const auto& values) {
      std::set<std::string> unique;
      return std::all_of(values.begin(), values.end(), [&](const auto& value) {
        return !value.empty() && unique.insert(value).second;
      });
    };
    if (!unique_nonempty(source_ids) ||
        !unique_nonempty(exact_linear_map_ids))
      return false;
    if (model == Model::IndependentHull)
      return lossy_hull_explicit && exact_linear_map_ids.empty();
    if (lossy_hull_explicit) return false;
    if (model == Model::Exact) return source_ids.empty();
    return !source_ids.empty();
  }

  friend bool operator==(const Descriptor&, const Descriptor&) = default;
};

inline Descriptor exact(std::string identity) {
  Descriptor output;
  output.model = Model::Exact;
  output.identity = std::move(identity);
  output.lossy_hull_explicit = false;
  if (!output.valid())
    throw std::invalid_argument("invalid exact correlation identity");
  return output;
}

inline Descriptor shared_linear_sources(
    std::string identity, std::vector<std::string> source_ids) {
  Descriptor output;
  output.model = Model::SharedLinearSources;
  output.identity = std::move(identity);
  output.source_ids = std::move(source_ids);
  output.lossy_hull_explicit = false;
  if (!output.valid())
    throw std::invalid_argument("invalid shared-source correlation record");
  return output;
}

inline Descriptor exact_pushforward(const Descriptor& input,
                                    std::string output_identity,
                                    std::string exact_map_identity) {
  if (!input.valid() || input.model == Model::IndependentHull ||
      output_identity.empty() || exact_map_identity.empty())
    throw std::invalid_argument(
        "exact correlation pushforward requires a correlated input and exact map");
  Descriptor output = input;
  output.identity = std::move(output_identity);
  output.exact_linear_map_ids.push_back(std::move(exact_map_identity));
  if (!output.valid())
    throw std::invalid_argument("correlation pushforward is not canonical");
  return output;
}

inline Descriptor independent_hull(const Descriptor& input,
                                   std::string output_identity) {
  if (!input.valid() || output_identity.empty())
    throw std::invalid_argument(
        "independent hull conversion requires an input and output identity");
  Descriptor output;
  output.model = Model::IndependentHull;
  output.identity = std::move(output_identity);
  output.source_ids = input.source_ids;
  output.lossy_hull_explicit = true;
  if (!output.valid())
    throw std::invalid_argument("independent hull record is not canonical");
  return output;
}

// Exact/shared correlation is invariant information.  An independent hull
// can satisfy only an explicitly independent demand; it can never satisfy a
// shared-source demand even when its numerical radii are smaller.
inline bool dominates(const Descriptor& guarantee,
                      const Descriptor& demand) {
  if (!guarantee.valid() || !demand.valid()) return false;
  if (demand.model == Model::IndependentHull)
    return guarantee.model == Model::IndependentHull &&
           guarantee.identity == demand.identity;
  return guarantee.model == demand.model &&
         guarantee.identity == demand.identity &&
         guarantee.source_ids == demand.source_ids &&
         guarantee.exact_linear_map_ids == demand.exact_linear_map_ids;
}

inline boost::json::object to_json(const Descriptor& descriptor) {
  if (!descriptor.valid())
    throw std::invalid_argument("cannot encode invalid correlation descriptor");
  boost::json::array sources;
  for (const auto& source : descriptor.source_ids)
    sources.push_back(boost::json::value(source));
  boost::json::array maps;
  for (const auto& map : descriptor.exact_linear_map_ids)
    maps.push_back(boost::json::value(map));
  return boost::json::object{
      {"schema", descriptor.schema},
      {"model", model_name(descriptor.model)},
      {"identity", descriptor.identity},
      {"source_ids", std::move(sources)},
      {"exact_linear_map_ids", std::move(maps)},
      {"lossy_hull_explicit", descriptor.lossy_hull_explicit}};
}

inline Descriptor from_json(const boost::json::value& value) {
  if (!value.is_object())
    throw std::invalid_argument("correlation descriptor must be an object");
  const auto& object = value.as_object();
  const std::set<std::string> expected{
      "schema", "model", "identity", "source_ids",
      "exact_linear_map_ids", "lossy_hull_explicit"};
  std::set<std::string> actual;
  for (const auto& entry : object) actual.emplace(entry.key());
  if (actual != expected)
    throw std::invalid_argument("correlation descriptor has unexpected keys");
  if (!object.at("schema").is_string() ||
      !object.at("model").is_string() ||
      !object.at("identity").is_string() ||
      !object.at("source_ids").is_array() ||
      !object.at("exact_linear_map_ids").is_array() ||
      !object.at("lossy_hull_explicit").is_bool())
    throw std::invalid_argument("correlation descriptor has invalid field types");
  Descriptor output;
  output.schema = std::string(object.at("schema").as_string());
  output.model = parse_model(std::string(object.at("model").as_string()));
  output.identity = std::string(object.at("identity").as_string());
  for (const auto& source : object.at("source_ids").as_array()) {
    if (!source.is_string())
      throw std::invalid_argument("correlation source identity is not a string");
    output.source_ids.emplace_back(source.as_string());
  }
  for (const auto& map : object.at("exact_linear_map_ids").as_array()) {
    if (!map.is_string())
      throw std::invalid_argument("correlation map identity is not a string");
    output.exact_linear_map_ids.emplace_back(map.as_string());
  }
  output.lossy_hull_explicit =
      object.at("lossy_hull_explicit").as_bool();
  if (!output.valid())
    throw std::invalid_argument("correlation descriptor is inconsistent");
  return output;
}

struct ResidualBinding {
  std::string identity;
  Descriptor correlation;
  std::string equation_operator_identity;
  std::string input_frame_identity;
  std::string residual_frame_identity;
  std::string theorem_identity;

  [[nodiscard]] bool valid() const {
    return !identity.empty() && correlation.valid() &&
           !equation_operator_identity.empty() &&
           !input_frame_identity.empty() && !residual_frame_identity.empty() &&
           !theorem_identity.empty();
  }
};

}  // namespace diffexp::kernel::correlation
