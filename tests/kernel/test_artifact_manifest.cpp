#include "diffexp/kernel/artifact_manifest.hpp"

#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using diffexp::kernel::artifact::LaurentSupport;
using diffexp::kernel::artifact::Manifest;
using diffexp::kernel::artifact::ResourceContract;
using diffexp::kernel::artifact::SemanticKey;

namespace {

SemanticKey key(std::string branch = "+i0") {
  return SemanticKey{diffexp::kernel::artifact::kSemanticKeySchema,
                     "receiving-basis",
                     boost::json::object{{"family", "h48"},
                                         {"level", 3},
                                         {"branch", branch}},
                     "semantic-h48-" + branch};
}

ResourceContract contract(std::int32_t complete_max,
                          std::uint32_t digits = 90,
                          std::string branch = "+i0") {
  ResourceContract result;
  result.laurent = LaurentSupport{-2, complete_max};
  result.digits_lower = digits;
  result.taylor_complete_max = 80;
  result.observables = {"lower", "upper"};
  result.capabilities = {"match"};
  result.frame = "relative";
  result.correlation = "joint";
  result.branch = std::move(branch);
  result.theorem = "ordinary-recurrence-v1";
  result.master_ordering = "h48-order";
  return result;
}

Manifest artifact(const SemanticKey& semantic, ResourceContract guarantee,
                  std::string content,
                  std::vector<std::string> parents = {}) {
  return Manifest{diffexp::kernel::artifact::kArtifactManifestSchema,
                  semantic, std::move(guarantee), std::move(parents),
                  "payload-" + content, std::move(content)};
}

}  // namespace

int main() {
  const auto semantic = key();
  const auto demand = contract(21, 80);
  const auto narrow = artifact(semantic, contract(22, 90), "narrow");
  const auto wide = artifact(semantic, contract(30, 140), "wide");
  assert(narrow.valid());
  assert(diffexp::kernel::artifact::dominates(narrow.guarantee, demand));

  auto wrong_branch = contract(40, 200, "-i0");
  assert(!diffexp::kernel::artifact::dominates(wrong_branch, demand));
  auto independent = contract(40, 200);
  independent.correlation = "independent-hull";
  assert(!diffexp::kernel::artifact::dominates(independent, demand));

  const std::vector<Manifest> candidates{wide, narrow};
  const auto* selected = diffexp::kernel::artifact::select_least_dominator(
      candidates, semantic, demand);
  assert(selected != nullptr && selected->content_id == "narrow");

  const auto child = artifact(key(), contract(22), "child", {"narrow"});
  const auto grandchild =
      artifact(key(), contract(22), "grandchild", {"child"});
  const auto unrelated = artifact(key(), contract(22), "unrelated");
  const auto invalid = diffexp::kernel::artifact::descendant_invalidation_set(
      {narrow, child, grandchild, unrelated}, {"narrow"});
  assert(invalid == std::set<std::string>({"narrow", "child", "grandchild"}));

  std::cout << "artifact manifest contracts passed\n";
  return 0;
}
