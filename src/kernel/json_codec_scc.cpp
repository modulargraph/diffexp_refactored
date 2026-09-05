#include "json_codec_scc.hpp"

#include "json_codec_support.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace diffexp::kernel::json_codec_detail {
namespace json = boost::json;

namespace {

std::pair<std::uint32_t, std::uint32_t> parse_edge(
    const json::value& value, std::uint32_t bound, const char* label) {
  const auto& edge = as_array(value, label);
  if (edge.size() != 2)
    throw std::invalid_argument(std::string(label) +
                                " must have two vertices");
  const auto source = as_u32(edge[0], "edge source");
  const auto target = as_u32(edge[1], "edge target");
  if (source >= bound || target >= bound)
    throw std::invalid_argument(std::string(label) +
                                " vertex is outside range");
  return {source, target};
}

}  // namespace

SCCCertificate validate_scc_certificate(const json::value& raw,
                                        std::uint32_t dimension) {
  const auto& object = as_object(raw, "SCC certificate");
  const auto& raw_components = as_array(
      object.at("components"), "SCC components");
  if (raw_components.empty())
    throw std::invalid_argument("SCC certificate has no components");

  std::vector<std::vector<std::uint32_t>> components;
  std::vector<std::uint32_t> component_of(
      dimension, std::numeric_limits<std::uint32_t>::max());
  for (std::uint32_t block = 0; block < raw_components.size(); ++block) {
    std::vector<std::uint32_t> component;
    for (const auto& raw_vertex :
         as_array(raw_components[block], "SCC component")) {
      const auto vertex = as_u32(raw_vertex, "SCC vertex");
      if (vertex >= dimension ||
          component_of[vertex] !=
              std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(
            "SCC components must partition the recurrence dimension");
      component_of[vertex] = block;
      component.push_back(vertex);
    }
    if (component.empty())
      throw std::invalid_argument(
          "SCC certificate contains an empty component");
    std::sort(component.begin(), component.end());
    components.push_back(std::move(component));
  }
  if (std::any_of(component_of.begin(), component_of.end(), [](auto value) {
        return value == std::numeric_limits<std::uint32_t>::max();
      }))
    throw std::invalid_argument(
        "SCC components do not cover the recurrence dimension");

  std::vector<std::vector<std::uint32_t>> adjacency(dimension);
  std::vector<std::vector<std::uint32_t>> reverse(dimension);
  std::set<std::pair<std::uint32_t, std::uint32_t>> structural_edges;
  for (const auto& raw_edge : as_array(
           object.at("structural_edges"), "exact structural edges")) {
    const auto edge = parse_edge(
        raw_edge, dimension, "exact structural edge");
    if (structural_edges.insert(edge).second) {
      adjacency[edge.first].push_back(edge.second);
      reverse[edge.second].push_back(edge.first);
    }
  }

  std::vector<std::uint8_t> seen(dimension, 0);
  std::vector<std::uint32_t> order;
  std::function<void(std::uint32_t)> visit = [&](std::uint32_t vertex) {
    seen[vertex] = 1;
    for (const auto next : adjacency[vertex])
      if (!seen[next]) visit(next);
    order.push_back(vertex);
  };
  for (std::uint32_t vertex = 0; vertex < dimension; ++vertex)
    if (!seen[vertex]) visit(vertex);
  std::fill(seen.begin(), seen.end(), 0);
  std::vector<std::vector<std::uint32_t>> derived_components;
  std::function<void(std::uint32_t, std::vector<std::uint32_t>&)> collect =
      [&](std::uint32_t vertex, std::vector<std::uint32_t>& component) {
        seen[vertex] = 1;
        component.push_back(vertex);
        for (const auto next : reverse[vertex])
          if (!seen[next]) collect(next, component);
      };
  for (auto iterator = order.rbegin(); iterator != order.rend(); ++iterator) {
    if (seen[*iterator]) continue;
    std::vector<std::uint32_t> component;
    collect(*iterator, component);
    std::sort(component.begin(), component.end());
    derived_components.push_back(std::move(component));
  }
  const auto canonical_components = [](auto value) {
    std::sort(value.begin(), value.end());
    return value;
  };
  if (canonical_components(components) !=
      canonical_components(derived_components))
    throw std::invalid_argument(
        "SCC partition does not match the exact structural graph");

  const auto block_count = static_cast<std::uint32_t>(components.size());
  std::set<std::pair<std::uint32_t, std::uint32_t>> derived_condensation;
  for (const auto [source, target] : structural_edges) {
    const auto from = component_of[source];
    const auto to = component_of[target];
    if (from != to) derived_condensation.insert({from, to});
  }
  std::set<std::pair<std::uint32_t, std::uint32_t>> supplied_condensation;
  for (const auto& raw_edge : as_array(
           object.at("condensation_edges"), "condensation edges"))
    supplied_condensation.insert(
        parse_edge(raw_edge, block_count, "condensation edge"));
  if (derived_condensation != supplied_condensation)
    throw std::invalid_argument(
        "SCC condensation edges do not match the exact structural graph");

  std::vector<std::uint32_t> topological_order;
  std::vector<std::uint32_t> topological_position(
      block_count, std::numeric_limits<std::uint32_t>::max());
  for (const auto& raw_block : as_array(
           object.at("topological_order"), "SCC topological order")) {
    const auto block = as_u32(raw_block, "SCC topological block");
    if (block >= block_count ||
        topological_position[block] !=
            std::numeric_limits<std::uint32_t>::max())
      throw std::invalid_argument("invalid SCC topological order");
    topological_position[block] =
        static_cast<std::uint32_t>(topological_order.size());
    topological_order.push_back(block);
  }
  if (topological_order.size() != block_count)
    throw std::invalid_argument("SCC topological order is incomplete");
  for (const auto [source, target] : derived_condensation)
    if (topological_position[source] >= topological_position[target])
      throw std::invalid_argument(
          "SCC topological order contradicts a condensation edge");

  std::vector<std::uint32_t> depth(block_count, 0);
  for (const auto block : topological_order) {
    for (const auto [source, target] : derived_condensation)
      if (source == block)
        depth[target] = std::max(depth[target], depth[source] + 1);
  }
  const auto coupling_depth = depth.empty()
      ? 0
      : *std::max_element(depth.begin(), depth.end());
  if (as_u32(object.at("coupling_depth"), "SCC coupling depth") !=
      coupling_depth)
    throw std::invalid_argument(
        "SCC coupling depth does not match the exact condensation graph");

  SCCCertificate certificate;
  certificate.component_count = block_count;
  certificate.coupling_depth = coupling_depth;
  certificate.components = std::move(components);
  certificate.component_of = std::move(component_of);
  certificate.structural_edges.assign(
      structural_edges.begin(), structural_edges.end());
  certificate.condensation_edges.assign(
      derived_condensation.begin(), derived_condensation.end());
  certificate.topological_order = std::move(topological_order);
  certificate.exact_record = json::serialize(raw);
  return certificate;
}

json::array encode_indices(const std::vector<std::uint32_t>& values) {
  json::array encoded;
  encoded.reserve(values.size());
  for (const auto value : values) encoded.push_back(value);
  return encoded;
}

}  // namespace diffexp::kernel::json_codec_detail
