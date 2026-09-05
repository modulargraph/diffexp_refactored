#include "diffexp/kernel/checkpoint.hpp"
#include "diffexp/kernel/correlation.hpp"

#include <boost/json.hpp>

#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>

int main() {
  namespace correlation = diffexp::kernel::correlation;
  auto reduced = correlation::shared_linear_sources(
      "b4-reduced-tail", {"common-seed", "exact-fuchsian-tail"});
  auto physical = correlation::exact_pushforward(
      reduced, "b4-physical-tail", "nonidentity-shear");
  physical = correlation::exact_pushforward(
      physical, "b4-residual-frame", "spectral-normal-frame");
  assert(physical.source_ids == reduced.source_ids);
  assert(physical.exact_linear_map_ids.size() == 2);

  auto restored = correlation::from_json(correlation::to_json(physical));
  assert(restored == physical);
  assert(correlation::dominates(restored, physical));

  auto hull = correlation::independent_hull(physical, "b4-lossy-hull");
  assert(hull.model == correlation::Model::IndependentHull);
  assert(!correlation::dominates(hull, physical));

  correlation::ResidualBinding residual{
      "b4-residual-binding", physical, "physical-qc-operator",
      "reduced-fuchsian-frame", "normal-residual-frame",
      "full-ball-forward-residual-v1"};
  assert(residual.valid());

  const auto path = std::filesystem::temp_directory_path() /
      ("diffexp3-correlation-" + std::to_string(::getpid()) + ".de2cp");
  diffexp::kernel::checkpoint::write_atomic(
      path.string(), boost::json::object{{"schema", 1}},
      correlation::to_json(physical));
  const auto checkpoint = diffexp::kernel::checkpoint::read_json(path.string());
  const auto checkpoint_correlation =
      correlation::from_json(checkpoint.payload);
  assert(checkpoint_correlation == physical);
  std::filesystem::remove(path);

  std::cout << "correlation artifact contracts passed\n";
  return 0;
}
