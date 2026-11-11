# --boundary: a level whose endpoint qualifies matches there by default, the
# result is that of the interior anchors, and a family without such a level is
# left exactly as it was.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef run_id)
set(work "${CMAKE_CURRENT_BINARY_DIR}/cli-endpoint-boundary-${run_id}")
file(MAKE_DIRECTORY "${work}")
set(examples "${DIFFEXP_SOURCE_DIR}/examples/feynman")

function(run_ft family tag)
  execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${examples}/${family}.json" --cache "${work}/cache-${tag}" --no-numerical-cache --json ${ARGN}
    OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
  set(ft_output "${output}" PARENT_SCOPE)
  set(ft_errors "${errors}" PARENT_SCOPE)
  set(ft_status "${status}" PARENT_SCOPE)
endfunction()
# The first `digits` characters of the real part of the finite coefficient.
function(finite_part variable digits)
  if(NOT ft_status EQUAL 0)
    message(FATAL_ERROR "run failed: ${ft_errors} ${ft_output}")
  endif()
  string(JSON low GET "${ft_output}" epsilon_low)
  math(EXPR index "0-(${low})")
  string(JSON real GET "${ft_output}" coefficients 0 ${index} real)
  string(SUBSTRING "${real}" 0 ${digits} leading)
  set(${variable} "${leading}" PARENT_SCOPE)
endfunction()
function(boundaries variable)
  string(JSON count LENGTH "${ft_output}" ladder)
  set(found "")
  math(EXPR last "${count}-1")
  foreach(level RANGE ${last})
    string(JSON site GET "${ft_output}" ladder ${level} boundary)
    string(APPEND found "${site};")
  endforeach()
  set(${variable} "${found}" PARENT_SCOPE)
endfunction()

# The massive kite: the first two levels merge massive propagators on one
# loop momentum and match at x = 0; the third merges different loop momenta.
run_ft(kite auto)
finite_part(endpoint_value 30)
boundaries(endpoint_sites)
string(JSON endpoint_masters GET "${ft_output}" ladder 2 masters)
if(NOT endpoint_sites STREQUAL "lower;lower;anchor;" OR NOT ft_errors MATCHES "Level 1: boundary at the endpoint x = 0")
  message(FATAL_ERROR "kite: the default did not take the endpoints: ${endpoint_sites} ${ft_errors}")
endif()
run_ft(kite anchor --boundary anchor)
finite_part(anchor_value 30)
boundaries(anchor_sites)
string(JSON anchor_masters GET "${ft_output}" ladder 2 masters)
if(NOT anchor_sites STREQUAL "anchor;anchor;anchor;")
  message(FATAL_ERROR "kite: --boundary anchor took an endpoint: ${anchor_sites}")
endif()
# "[0." and 27 digits.
if(NOT endpoint_value STREQUAL anchor_value OR NOT endpoint_value MATCHES "^\\[0\\.223983919107442530381929718$")
  message(FATAL_ERROR "kite: endpoint ${endpoint_value} and anchor ${anchor_value} runs differ")
endif()
# The child of the second level is a pinched family at its original
# kinematics: its level has fewer masters than at an interior anchor.
if(NOT endpoint_masters LESS anchor_masters)
  message(FATAL_ERROR "kite: the pinched level has ${endpoint_masters} masters, the anchored one ${anchor_masters}")
endif()

# The massless box has no qualifying level: nothing changes, and an endpoint
# that was demanded is an error.
run_ft(box auto)
boundaries(box_sites)
if(NOT ft_status EQUAL 0 OR NOT box_sites STREQUAL "anchor;anchor;")
  message(FATAL_ERROR "box: a massless family took an endpoint: ${box_sites} ${ft_errors}")
endif()
run_ft(box endpoint --boundary endpoint)
if(ft_status EQUAL 0 OR NOT ft_errors MATCHES "no level of this family qualifies")
  message(FATAL_ERROR "box: --boundary endpoint was not refused: ${ft_errors}")
endif()
run_ft(box unknown --boundary sideways)
if(ft_status EQUAL 0 OR NOT ft_errors MATCHES "--boundary must be anchor, endpoint or auto")
  message(FATAL_ERROR "an unknown --boundary value was accepted: ${ft_errors}")
endif()
file(REMOVE_RECURSE "${work}")
