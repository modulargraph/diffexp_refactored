#pragma once

#include <boost/json.hpp>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace diffexp::kernel::json_codec_detail {

struct SCCCertificate {
  std::uint32_t component_count = 0;
  std::uint32_t coupling_depth = 0;
  std::vector<std::vector<std::uint32_t>> components;
  std::vector<std::uint32_t> component_of;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> structural_edges;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> condensation_edges;
  std::vector<std::uint32_t> topological_order;
  std::string exact_record;
};

SCCCertificate validate_scc_certificate(const boost::json::value& raw,
                                        std::uint32_t dimension);
boost::json::array encode_indices(
    const std::vector<std::uint32_t>& values);

}  // namespace diffexp::kernel::json_codec_detail
