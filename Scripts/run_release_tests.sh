#!/usr/bin/env bash
# Sequential release verification. One Wolfram process is run at a time.

set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"

wolfram_tests=(
  Tests/test_tolerances.m
  Tests/test_config.m
  Tests/test_eps_series.m
  Tests/test_indicial.m
  Tests/test_sectorseries.m
  Tests/test_solve.m
  Tests/test_transport.m
  Tests/test_path_planner_algebraic.m
  Tests/test_native_algebraic_chart_bridge.m
  Tests/test_native_algebraic_scale_bridge.m
  Tests/test_native_regular_basis_batch.m
  Tests/test_native_regular_prepare_owner.m
  Tests/test_integrate.m
  Tests/test_api.m
  Tests/test_public_api.m
  Tests/test_paclet_metadata.m
  Tests/test_cpp_backend.m
  Tests/test_cpp_arm_batch.m
  Tests/test_feynmantrick_algebra.m
  Tests/test_feynmantrick_family_spec.m
  Tests/test_fast_tadpole_boundary.m
  Tests/test_ft_example_specs.m
  Tests/test_feynmantrick_failure_semantics.m
  Tests/test_fire7_modular_contract.m
  Tests/test_fire7_pipeline_plumbing.m
  Tests/test_fire7_rational_reconstruction.m
  Tests/test_fire7_lock_cache_contract.m
  Tests/test_fire7_classical_exact_output.m
  Tests/test_fire_inmemory_reduction_cache.m
  Tests/test_fire_level_reduction_batch.m
  Tests/test_ft_pipeline_facade.m
  Tests/test_ft_pipeline_family_request.m
  Tests/test_ft_pipeline_process.m
  Tests/test_ft_runner_family_request_contract.m
  Tests/test_ft_ladder_checkpointing.m
  Tests/test_double_box_planar_oracle.m
  Tests/test_pentagon_massive_oracle.m
  Tests/test_banana4_bessel_oracle.m
)

for test_file in "${wolfram_tests[@]}"; do
  echo "=== $test_file"
  if [[ "$test_file" == Tests/test_cpp_backend.m ||
        "$test_file" == Tests/test_cpp_arm_batch.m ]]; then
    DE2_REQUIRE_CPP=1 wolframscript -file "$test_file"
  else
    wolframscript -file "$test_file"
  fi
done

echo "=== direct examples"
wolframscript -file Examples/Direct/MinimalTransport.wl
wolframscript -file Examples/Direct/SingularEndpointAndSegments.wl

echo "=== Feynman-trick shell syntax"
for script in Examples/FeynmanTrick/*.sh; do
  sh -n "$script"
done

if [[ "${DE2_RUN_FIRE_TESTS:-0}" == 1 ]]; then
  echo "=== optional FIRE7 integration"
  wolframscript -file Tests/test_feynmantrick_fire.m
  wolframscript -file Tests/test_feynmantrick_iteration.m
  wolframscript -file Tests/test_fire7_real_bubble_parity.m
  wolframscript -file Tests/test_fire7_real_sunrise_parity.m
fi

if [[ "${DE2_RUN_BANANA_TIMING:-0}" == 1 ]]; then
  echo "=== optional full banana timing regression"
  Scripts/run_banana_timing_regression.sh
fi

if [[ "${DE2_RUN_BANANA4_BOUNDARY:-0}" == 1 ]]; then
  echo "=== optional banana4 production chain through level-2 boundary"
  Scripts/run_banana4_boundary_regression.sh
fi

if [[ "${DE2_RUN_BANANA4_FULL:-0}" == 1 ]]; then
  echo "=== optional full banana4 production chain and Bessel oracle"
  BANANA4_BOUNDARY_TARGET_LEVEL=final \
    Scripts/run_banana4_boundary_regression.sh
fi

if [[ "${DE2_RUN_PENTAGON_TIMING:-0}" == 1 ]]; then
  echo "=== optional full pentagon timing regression"
  Scripts/run_pentagon_timing_regression.sh
fi

if [[ "${DE2_RUN_DOUBLE_BOX_TIMING:-0}" == 1 ]]; then
  echo "=== optional full double-box timing regression"
  Scripts/run_double_box_timing_regression.sh
fi

echo "Release test suite passed."
