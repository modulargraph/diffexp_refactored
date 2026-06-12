#!/bin/zsh
# Sequential test battery (ONE Wolfram kernel license - never parallelize).
export WolframKernel="${WolframKernel:-/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel}"
cd "$(dirname "$0")/.."

tests=(
  "Tests/test_tolerances.m"
  "Tests/test_config.m"
  "Tests/test_eps_series.m"
  "Tests/test_indicial.m"
  "Tests/test_sectorseries.m"
  "Tests/test_solve.m"
  "Tests/test_transport.m"
  "Tests/test_integrate.m"
  "Tests/test_api.m"
  "Tests/test_pins.m"
  Tests/test_package_loading.m
  Tests/test_symbols_namespacing.m
  Tests/test_recurrence_no_fallback.m
  Tests/test_local_series_edge_cases.m
  Tests/test_resonant_2f1.m
  Tests/test_singular_recurrence.m
  Tests/test_singular_endpoint.m
  Tests/test_regularized_integration.m
  Tests/test_regularized_integration_edge_cases.m
  Tests/test_interior_singular_integration.m
  Tests/test_integration_log_depth.m
  Tests/test_multisector_fit.m
  Tests/test_feynmantrick_algebra.m
  Tests/test_feynmantrick_fire.m
  Tests/test_feynmantrick_iteration.m
  Tests/test_feynmantrick_integration_edge_cases.m
  Tests/test_feynmantrick_pipeline.m
  Tests/test_feynmantrick_endtoend.m
)

overall=0
for t in "${tests[@]}"; do
  echo "=== $t"
  out=$(wolframscript -file "$t" 2>&1)
  code=$?
  summary=$(echo "$out" | grep -E "tests passed|PASSED|FAILED|Results:" | tail -3)
  echo "$summary"
  if [ $code -ne 0 ] || echo "$out" | grep -Eq "Some tests FAILED|SOME FAILED|TEST FAILED"; then
    echo "FAILED: $t (exit $code)"
    echo "$out" | tail -30
    overall=1
  fi
done
echo "=== battery exit: $overall"
exit $overall
