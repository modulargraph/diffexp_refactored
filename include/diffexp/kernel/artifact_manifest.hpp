#pragma once

#include "diffexp/kernel/scalar.hpp"

#include <boost/json/value.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp::kernel::artifact {

inline constexpr const char* kSemanticKeySchema =
    "DiffExp3.SemanticKey/v1";
inline constexpr const char* kArtifactManifestSchema =
    "DiffExp3.ArtifactManifest/v1";

struct LaurentSupport {
  std::int32_t min_power = 0;
  std::int32_t complete_max = -1;

  [[nodiscard]] bool valid() const noexcept {
    return min_power <= complete_max;
  }

  [[nodiscard]] bool dominates(const LaurentSupport& demand) const noexcept {
    return valid() && demand.valid() && min_power <= demand.min_power &&
           complete_max >= demand.complete_max;
  }
};

// Equality fields name mathematical semantics.  Resource fields form a
// product order.  They are deliberately separate so a larger numerical
// object can never compensate for a different branch, frame, correlation
// model, theorem, geometry, scalar domain, or master ordering.
struct ResourceContract {
  std::optional<LaurentSupport> laurent;
  std::map<std::string, LaurentSupport> component_laurent;
  std::optional<std::uint32_t> digits_lower;
  std::optional<Rational> absolute_radius_upper;
  std::optional<Rational> relative_radius_upper;
  std::optional<std::uint32_t> taylor_complete_max;
  std::optional<Rational> certified_radius_lower;
  std::set<std::string> observables;
  std::set<std::string> capabilities;

  std::string frame;
  std::string correlation;
  std::string branch;
  std::string geometry;
  std::string theorem;
  std::string proof_scope;
  std::string scalar_domain;
  std::string master_ordering;
};

inline bool optional_equality_axis_dominates(const std::string& guarantee,
                                             const std::string& demand) {
  return demand.empty() || guarantee == demand;
}

template <class T>
inline bool optional_lower_bound_dominates(const std::optional<T>& guarantee,
                                           const std::optional<T>& demand) {
  return !demand.has_value() ||
         (guarantee.has_value() && *guarantee >= *demand);
}

template <class T>
inline bool optional_upper_bound_dominates(const std::optional<T>& guarantee,
                                           const std::optional<T>& demand) {
  return !demand.has_value() ||
         (guarantee.has_value() && *guarantee <= *demand);
}

inline bool set_dominates(const std::set<std::string>& guarantee,
                          const std::set<std::string>& demand) {
  return std::includes(guarantee.begin(), guarantee.end(), demand.begin(),
                       demand.end());
}

inline bool dominates(const ResourceContract& guarantee,
                      const ResourceContract& demand) {
  if (demand.laurent.has_value() &&
      (!guarantee.laurent.has_value() ||
       !guarantee.laurent->dominates(*demand.laurent)))
    return false;
  for (const auto& [component, needed] : demand.component_laurent) {
    const auto found = guarantee.component_laurent.find(component);
    if (found == guarantee.component_laurent.end() ||
        !found->second.dominates(needed))
      return false;
  }
  return optional_lower_bound_dominates(guarantee.digits_lower,
                                        demand.digits_lower) &&
         optional_upper_bound_dominates(guarantee.absolute_radius_upper,
                                        demand.absolute_radius_upper) &&
         optional_upper_bound_dominates(guarantee.relative_radius_upper,
                                        demand.relative_radius_upper) &&
         optional_lower_bound_dominates(guarantee.taylor_complete_max,
                                        demand.taylor_complete_max) &&
         optional_lower_bound_dominates(guarantee.certified_radius_lower,
                                        demand.certified_radius_lower) &&
         set_dominates(guarantee.observables, demand.observables) &&
         set_dominates(guarantee.capabilities, demand.capabilities) &&
         optional_equality_axis_dominates(guarantee.frame, demand.frame) &&
         optional_equality_axis_dominates(guarantee.correlation,
                                          demand.correlation) &&
         optional_equality_axis_dominates(guarantee.branch, demand.branch) &&
         optional_equality_axis_dominates(guarantee.geometry,
                                          demand.geometry) &&
         optional_equality_axis_dominates(guarantee.theorem,
                                          demand.theorem) &&
         optional_equality_axis_dominates(guarantee.proof_scope,
                                          demand.proof_scope) &&
         optional_equality_axis_dominates(guarantee.scalar_domain,
                                          demand.scalar_domain) &&
         optional_equality_axis_dominates(guarantee.master_ordering,
                                          demand.master_ordering);
}

struct SemanticKey {
  std::string schema = kSemanticKeySchema;
  std::string domain;
  boost::json::object record;
  std::string identity;

  [[nodiscard]] bool valid() const {
    return schema == kSemanticKeySchema && !domain.empty() &&
           !identity.empty();
  }

  friend bool operator==(const SemanticKey& left, const SemanticKey& right) {
    return left.schema == right.schema && left.domain == right.domain &&
           left.record == right.record && left.identity == right.identity;
  }
};

struct Manifest {
  std::string schema = kArtifactManifestSchema;
  SemanticKey semantic_key;
  ResourceContract guarantee;
  std::vector<std::string> parent_content_ids;
  std::string payload_identity;
  std::string content_id;

  [[nodiscard]] bool valid() const {
    if (schema != kArtifactManifestSchema || !semantic_key.valid() ||
        payload_identity.empty() || content_id.empty())
      return false;
    std::set<std::string> unique;
    for (const auto& parent : parent_content_ids)
      if (parent.empty() || !unique.insert(parent).second) return false;
    return true;
  }
};

inline bool satisfies(const Manifest& artifact, const SemanticKey& key,
                      const ResourceContract& demand) {
  return artifact.valid() && artifact.semantic_key == key &&
         dominates(artifact.guarantee, demand);
}

inline std::int64_t laurent_excess(const ResourceContract& guarantee,
                                   const ResourceContract& demand) {
  std::int64_t score = 0;
  if (guarantee.laurent && demand.laurent)
    score += static_cast<std::int64_t>(demand.laurent->min_power) -
             guarantee.laurent->min_power +
             static_cast<std::int64_t>(guarantee.laurent->complete_max) -
             demand.laurent->complete_max;
  if (guarantee.digits_lower && demand.digits_lower)
    score += *guarantee.digits_lower - *demand.digits_lower;
  if (guarantee.taylor_complete_max && demand.taylor_complete_max)
    score += *guarantee.taylor_complete_max - *demand.taylor_complete_max;
  score += static_cast<std::int64_t>(guarantee.observables.size() -
                                     demand.observables.size());
  score += static_cast<std::int64_t>(guarantee.capabilities.size() -
                                     demand.capabilities.size());
  return score;
}

inline const Manifest* select_least_dominator(
    const std::vector<Manifest>& artifacts, const SemanticKey& key,
    const ResourceContract& demand) {
  const Manifest* selected = nullptr;
  std::int64_t selected_score = 0;
  for (const auto& artifact : artifacts) {
    if (!satisfies(artifact, key, demand)) continue;
    const auto score = laurent_excess(artifact.guarantee, demand);
    if (selected == nullptr || score < selected_score ||
        (score == selected_score && artifact.content_id < selected->content_id)) {
      selected = &artifact;
      selected_score = score;
    }
  }
  return selected;
}

inline std::set<std::string> descendant_invalidation_set(
    const std::vector<Manifest>& artifacts,
    const std::set<std::string>& changed) {
  std::set<std::string> invalid = changed;
  bool grew = true;
  while (grew) {
    grew = false;
    for (const auto& artifact : artifacts) {
      if (!artifact.valid() || invalid.contains(artifact.content_id)) continue;
      if (std::any_of(artifact.parent_content_ids.begin(),
                      artifact.parent_content_ids.end(),
                      [&](const auto& parent) { return invalid.contains(parent); })) {
        invalid.insert(artifact.content_id);
        grew = true;
      }
    }
  }
  return invalid;
}

}  // namespace diffexp::kernel::artifact
