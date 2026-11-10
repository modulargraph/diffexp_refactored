execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" --version
  OUTPUT_VARIABLE version RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT version MATCHES "DiffExp 2.1.1")
  message(FATAL_ERROR "Version command failed: ${version}")
endif()
set(request "${DIFFEXP_SOURCE_DIR}/examples/logarithm.json")
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" series "${request}"
  OUTPUT_VARIABLE from_file ERROR_VARIABLE errors RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT errors STREQUAL "")
  message(FATAL_ERROR "File request failed: ${errors}")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" series - INPUT_FILE "${request}"
  OUTPUT_VARIABLE from_stdin ERROR_VARIABLE errors RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT from_file STREQUAL from_stdin)
  message(FATAL_ERROR "File/stdin responses differ: ${errors}")
endif()
string(JSON coefficient GET "${from_stdin}" coefficients 3 0 1)
if(NOT coefficient STREQUAL "1/3")
  message(FATAL_ERROR "Exact logarithm coefficient changed: ${coefficient}")
endif()
string(JSON schema GET "${from_stdin}" schema)
if(NOT schema STREQUAL "DiffExp.TaylorSeries/v1")
  message(FATAL_ERROR "Unexpected series schema: ${schema}")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" series -
  INPUT_FILE "${DIFFEXP_SOURCE_DIR}/LICENSE"
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
if(status EQUAL 0 OR NOT output STREQUAL "" OR errors STREQUAL "")
  message(FATAL_ERROR "Invalid JSON was not rejected cleanly")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft not-a-family --json
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON error_status GET "${output}" status)
if(status EQUAL 0 OR NOT error_status STREQUAL "error")
  message(FATAL_ERROR "FT failure is not a JSON error: ${output}")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" backend-info
  OUTPUT_VARIABLE output RESULT_VARIABLE status)
string(JSON standalone GET "${output}" standalone_cpp)
if(NOT status EQUAL 0 OR NOT standalone)
  message(FATAL_ERROR "Backend is not standalone C++")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" transport "${DIFFEXP_SOURCE_DIR}/examples/transport/ordinary.json"
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON schema GET "${output}" schema)
if(NOT status EQUAL 0 OR NOT schema STREQUAL "DiffExp.TransportResult/v1")
  message(FATAL_ERROR "Generic transport failed: ${errors} ${output}")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/bubble.json" --epsilon-order 0 --json
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON completed GET "${output}" status)
string(JSON timing GET "${output}" timings total_seconds)
if(NOT status EQUAL 0 OR NOT completed STREQUAL "completed" OR timing LESS 0)
  message(FATAL_ERROR "Configured FT dispatch or timing failed: ${errors} ${output}")
endif()

# This mode must preserve the answer without creating numerical checkpoints.
string(JSON cached_coefficients GET "${output}" coefficients)
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef transient_id)
set(transient_cache "${CMAKE_CURRENT_BINARY_DIR}/cli-transient-${transient_id}")
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/bubble.json"
  --epsilon-order 0 --cache "${transient_cache}" --no-numerical-cache --json
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON transient_coefficients GET "${output}" coefficients)
string(JSON persistent GET "${output}" persistent_numerical_cache)
if(NOT status EQUAL 0 OR persistent OR NOT cached_coefficients STREQUAL transient_coefficients)
  message(FATAL_ERROR "Transient FT mode changed the result: ${errors} ${output}")
endif()
if(EXISTS "${transient_cache}/affine-series" OR EXISTS "${transient_cache}/ordinary-transport")
  message(FATAL_ERROR "Transient FT mode created numerical checkpoint directories")
endif()
file(REMOVE_RECURSE "${transient_cache}")

# Complete native IBP recursion must actually call the bundled solver on a cold cache.
set(ibp_cache "${CMAKE_CURRENT_BINARY_DIR}/cli-ibp-${transient_id}")
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/sunrise.json"
  --ibp-provider ibp-solver --cache "${ibp_cache}" --no-numerical-cache --json
  OUTPUT_VARIABLE native ERROR_VARIABLE errors RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Native IBP recursion failed: ${errors} ${native}")
endif()
string(JSON provider GET "${native}" ibp_provider)
string(JSON probes GET "${native}" ibp_statistics fresh_probes)
if(NOT provider STREQUAL "ibp-solver" OR probes LESS 1)
  message(FATAL_ERROR "Native IBP provider was not executed")
endif()
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/sunrise.json"
  --cache "${ibp_cache}" --no-numerical-cache --json
  OUTPUT_VARIABLE resumed ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON reused GET "${resumed}" systems_reused)
string(JSON probes GET "${resumed}" ibp_statistics fresh_probes)
if(NOT status EQUAL 0 OR reused LESS 1 OR NOT probes EQUAL 0)
  message(FATAL_ERROR "Native IBP closure cache was not reused")
endif()
# A reused closure ignores the IBP provider and seed options. That must be
# visible, on stderr and in the report, and absent on the cold run.
string(JSON cold_ignored GET "${native}" ibp_cache provider_options_ignored)
string(JSON warm_ignored GET "${resumed}" ibp_cache provider_options_ignored)
string(JSON warm_reused GET "${resumed}" ibp_cache systems_reused)
if(cold_ignored OR NOT warm_ignored OR NOT warm_reused EQUAL reused OR
    NOT errors MATCHES "reused from the verified cache" OR NOT errors MATCHES "did not apply")
  message(FATAL_ERROR "IBP cache reuse did not report the ignored provider options: ${errors} ${resumed}")
endif()
# The phases preceding transport are reported separately.
string(JSON series_seconds GET "${native}" timings endpoint_series_seconds)
string(JSON operator_seconds GET "${native}" timings endpoint_operator_seconds)
string(JSON selection_seconds GET "${native}" timings basis_selection_seconds)
string(JSON skipped_scans GET "${native}" basis_selection scans_skipped)
if(series_seconds LESS 0 OR operator_seconds LESS 0 OR selection_seconds LESS 0 OR skipped_scans LESS 0)
  message(FATAL_ERROR "FT phase timings missing: ${native}")
endif()
file(REMOVE_RECURSE "${ibp_cache}")

# Numerical options are rejected before the exact preparation starts: no IBP
# output, and not even the cache directory, may precede the error.
foreach(invalid "--precision-bits;32" "--leaf-digits;0" "--endpoint-order;0" "--ordinary-order;4")
  execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/sunrise.json"
    ${invalid} --cache "${ibp_cache}-invalid" --no-numerical-cache --json
    OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
  string(JSON error_status GET "${output}" status)
  if(status EQUAL 0 OR NOT error_status STREQUAL "error" OR errors MATCHES "IBP solver" OR EXISTS "${ibp_cache}-invalid")
    message(FATAL_ERROR "Invalid FT option ${invalid} was not rejected before preparation: ${errors} ${output}")
  endif()
endforeach()
# The command line accepts every ordinary solver the configuration accepts.
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/bubble.json"
  --epsilon-order 0 --ft-transport ultraspherical --no-numerical-cache --cache "${transient_cache}" --json
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
string(JSON completed GET "${output}" status)
if(NOT status EQUAL 0 OR NOT completed STREQUAL "completed")
  message(FATAL_ERROR "--ft-transport ultraspherical rejected: ${errors} ${output}")
endif()
# An executable the selected provider never runs is reported, not ignored silently.
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${DIFFEXP_SOURCE_DIR}/examples/feynman/bubble.json"
  --epsilon-order 0 --fire /nonexistent/FIRE --ibp-provider ibp-solver --no-numerical-cache --cache "${transient_cache}" --json
  OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT errors MATCHES "--fire is ignored")
  message(FATAL_ERROR "Ignored --fire was not reported: ${errors} ${output}")
endif()
file(REMOVE_RECURSE "${transient_cache}")
execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" --help OUTPUT_VARIABLE output RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT output MATCHES "--basis-scan-seconds" OR NOT output MATCHES "--no-numerical-cache keeps only the exact IBP level closures" OR
    NOT output MATCHES "ultraspherical")
  message(FATAL_ERROR "FT help omits documented options: ${output}")
endif()
