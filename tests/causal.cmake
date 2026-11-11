# The Feynman i0 is derived, not supplied: families in physical regions run
# from their definition alone and reproduce closed forms; a contour that
# cannot be derived is refused with a message that names the point.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef run_id)
set(work "${CMAKE_CURRENT_BINARY_DIR}/cli-causal-${run_id}")
file(MAKE_DIRECTORY "${work}")
set(examples "${DIFFEXP_SOURCE_DIR}/examples/feynman")

function(run_ft family_file)
  execute_process(COMMAND "${DIFFEXP_EXECUTABLE}" ft "${family_file}" --cache "${work}/cache" --no-numerical-cache --json ${ARGN}
    OUTPUT_VARIABLE output ERROR_VARIABLE errors RESULT_VARIABLE status)
  set(ft_output "${output}" PARENT_SCOPE)
  set(ft_errors "${errors}" PARENT_SCOPE)
  set(ft_status "${status}" PARENT_SCOPE)
endfunction()
# expect(<label> <epsilon power> <real prefix regex> <imaginary prefix regex>):
# the leading digits of one Laurent coefficient of the first component.
function(expect label power real imaginary)
  if(NOT ft_status EQUAL 0)
    message(FATAL_ERROR "${label}: run failed: ${ft_errors} ${ft_output}")
  endif()
  string(JSON low GET "${ft_output}" epsilon_low)
  math(EXPR index "${power}-(${low})")
  string(JSON found_real GET "${ft_output}" coefficients 0 ${index} real)
  string(JSON found_imaginary GET "${ft_output}" coefficients 0 ${index} imaginary)
  if(NOT found_real MATCHES "${real}" OR NOT found_imaginary MATCHES "${imaginary}")
    message(FATAL_ERROR "${label}: epsilon^${power} = ${found_real} + ${found_imaginary} i, expected ${real} + ${imaginary} i")
  endif()
endfunction()
function(expect_mode label mode)
  string(JSON found GET "${ft_output}" causal mode)
  if(NOT found STREQUAL "${mode}")
    message(FATAL_ERROR "${label}: causal mode ${found}, expected ${mode}")
  endif()
endfunction()

# Massless bubble at p^2=1, d=2-2eps: Gamma(1+eps) (-1-i0)^(-1-eps) B(-eps,-eps).
run_ft("${examples}/bubble_timelike.json" --epsilon-order 0)
expect("massless timelike bubble" -1 "^\\[?2\\.0000000000000000" "^(0|\\[\\+/-)")
expect("massless timelike bubble" 0 "^\\[-1\\.15443132980306572121" "^\\[6\\.28318530717958647692")
expect_mode("massless timelike bubble" derived)

# Equal-mass bubble at p^2=9 m^2, d=2-2eps: with beta=sqrt(5)/3 the finite
# term is 2/(9 beta) (log((1-beta)/(1+beta)) + i pi); the order-epsilon term
# is from an independent contour quadrature. The leaf passes two zeros of F.
run_ft("${examples}/bubble_above_threshold.json")
expect("bubble above threshold" 0 "^\\[-0\\.57387858795200538518591097" "^\\[0\\.93664196413876351908751661")
expect("bubble above threshold" 1 "^\\[-1\\.34953334547654819456090200" "^\\[-2\\.04811150156670114602461229")
expect_mode("bubble above threshold" derived)

# One-loop triangle whose first level needs opposite sides at its two cuts;
# the value does not depend on the cell of the level's anchor. Reference: an
# independent two-fold Feynman-parameter integral.
file(READ "${examples}/triangle_physical.json" triangle)
foreach(anchor "1/10" "1/2" "19/20")
  string(JSON moved SET "${triangle}" preparation "{\"anchors\":[\"${anchor}\",\"1/2\"]}")
  string(REPLACE "/" "o" tag "${anchor}")
  file(WRITE "${work}/triangle_${tag}.json" "${moved}")
  run_ft("${work}/triangle_${tag}.json")
  expect("triangle, anchor ${anchor}" 0 "^\\[0\\.1960034378504855817166" "^\\[0\\.2657392476831952020141")
  expect_mode("triangle, anchor ${anchor}" derived)
  string(JSON sides LENGTH "${ft_output}" causal oriented_points)
  if(sides LESS 2)
    message(FATAL_ERROR "triangle, anchor ${anchor}: the two cuts of level 1 were not oriented: ${ft_output}")
  endif()
endforeach()
string(JSON first GET "${ft_output}" causal oriented_points 0 side)
string(JSON last_index LENGTH "${ft_output}" causal oriented_points)
math(EXPR last_index "${last_index}-1")
string(JSON last GET "${ft_output}" causal oriented_points ${last_index} side)
if(NOT first EQUAL 1 OR NOT last EQUAL -1)
  message(FATAL_ERROR "triangle: level 1 must pass its first cut above and its last below: ${ft_output}")
endif()
# A contour height that would put the first two cuts under one detour: they
# ask for opposite sides and keep separate ones.
string(JSON tall SET "${triangle}" preparation "{\"anchors\":[\"1/10\",\"1/2\"]}")
string(JSON tall SET "${tall}" numerical "{\"contour_height\":\"1/2\"}")
file(WRITE "${work}/triangle_tall.json" "${tall}")
run_ft("${work}/triangle_tall.json")
expect("triangle, contour height 1/2" 0 "^\\[0\\.1960034378504855817166" "^\\[0\\.2657392476831952020141")
# The same through epsilon sampling, whose lifted path crosses the axis
# between the two cuts.
run_ft("${examples}/triangle_physical.json" --epsilon-strategy sample --sample-jobs 2)
expect("triangle, epsilon sampling" 0 "^\\[0\\.196003437850" "^\\[0\\.265739247683")
# An explicit object without level_signs still derives; with them it is the
# stated contour, here one side for both cuts, which is not the Feynman one.
string(JSON derived_rim SET "${triangle}" causal "{\"f_rim\":-1}")
file(WRITE "${work}/triangle_rim.json" "${derived_rim}")
run_ft("${work}/triangle_rim.json")
expect("triangle, causal object without level_signs" 0 "^\\[0\\.1960034378504855817166" "^\\[0\\.2657392476831952020141")
expect_mode("triangle, causal object without level_signs" derived)
string(JSON one_side SET "${triangle}" causal "{\"f_rim\":-1,\"level_signs\":[1,1],\"provenance\":\"test: one side for the whole level\"}")
string(JSON one_side SET "${one_side}" preparation "{\"anchors\":[\"1/2\",\"1/2\"]}")
file(WRITE "${work}/triangle_one_side.json" "${one_side}")
run_ft("${work}/triangle_one_side.json")
expect_mode("triangle, explicit level_signs" supplied)
string(JSON supplied_imaginary GET "${ft_output}" coefficients 0 1 imaginary)
if(supplied_imaginary MATCHES "^\\[0\\.2657392476831952020141")
  message(FATAL_ERROR "explicit level_signs were not applied as written: ${ft_output}")
endif()
# Sides stated point by point complete the same explicit contour: the last
# two cuts below the axis.
string(JSON stated SET "${one_side}" causal
  "{\"f_rim\":-1,\"level_signs\":[1,1],\"provenance\":\"test\",\"points\":[{\"level\":1,\"x\":\"2/3\",\"side\":-1},{\"level\":1,\"x\":\"0.872677996\",\"side\":-1}]}")
file(WRITE "${work}/triangle_points.json" "${stated}")
run_ft("${work}/triangle_points.json")
expect("triangle, explicit sides point by point" 0 "^\\[0\\.1960034378504855817166" "^\\[0\\.2657392476831952020141")
expect_mode("triangle, explicit sides point by point" supplied)
string(JSON misplaced SET "${triangle}" causal "{\"points\":[{\"level\":1,\"x\":\"0.3\",\"side\":-1}]}")
file(WRITE "${work}/triangle_misplaced.json" "${misplaced}")
run_ft("${work}/triangle_misplaced.json")
if(ft_status EQUAL 0 OR NOT ft_output MATCHES "level 1 has no real singular point at x = 3/10")
  message(FATAL_ERROR "a stated side that names no singular point was accepted: ${ft_output} ${ft_errors}")
endif()

# Massless box at s=1, t=-1/3, d=4-2eps:
# r_Gamma/(s t) {2/eps^2 [(-s)^-eps + (-t)^-eps] - log^2(s/t) - pi^2}, s -> s+i0.
# Default merges (the scalar leaf passes a zero of F), and a ladder whose
# second level has a real singular point of its own.
run_ft("${examples}/box_physical.json")
expect("physical box" -2 "^\\[-12\\.00000000000000000000" "^(0|\\[\\+/-)")
expect("physical box" -1 "^\\[0\\.33491424680973617890667" "^\\[-18\\.849555921538759430775")
expect("physical box" 0 "^\\[41\\.284167397546468572091" "^\\[-9\\.8280948169895949709241")
expect_mode("physical box" derived)
file(READ "${examples}/box_physical.json" box)
string(JSON ladder SET "${box}" preparation "{\"merges\":[[0,1],[1,2],[0,1]]}")
file(WRITE "${work}/box_ladder.json" "${ladder}")
run_ft("${work}/box_ladder.json")
expect("physical box, second ladder" -1 "^\\[0\\.33491424680973617890667" "^\\[-18\\.849555921538759430775")
expect("physical box, second ladder" 0 "^\\[41\\.284167397546468572091" "^\\[-9\\.8280948169895949709241")
string(JSON level GET "${ft_output}" causal oriented_points 0 level)
string(JSON side GET "${ft_output}" causal oriented_points 0 side)
if(NOT level EQUAL 2 OR NOT side EQUAL 1)
  message(FATAL_ERROR "physical box: level 2 must pass x=2/5 above the axis: ${ft_output}")
endif()

# A Euclidean family keeps the positivity fast path.
run_ft("${examples}/bubble.json" --epsilon-order 0)
expect_mode("Euclidean bubble" euclidean)

# At threshold the two zeros of F coincide: the contour is pinched.
file(READ "${examples}/bubble_above_threshold.json" bubble)
string(JSON threshold SET "${bubble}" external_gram "[[\"4\"]]")
file(WRITE "${work}/bubble_threshold.json" "${threshold}")
run_ft("${work}/bubble_threshold.json")
if(ft_status EQUAL 0 OR NOT ft_output MATCHES "repeated zero at x = 0\\.5000" OR NOT ft_output MATCHES "pinched")
  message(FATAL_ERROR "bubble at threshold was not refused as a pinched contour: ${ft_output} ${ft_errors}")
endif()
# A two-loop level above threshold is not decided by the implemented
# statements: no side is guessed.
file(READ "${examples}/sunrise.json" sunrise)
string(JSON above SET "${sunrise}" external_gram "[[\"20\"]]")
file(WRITE "${work}/sunrise_above.json" "${above}")
run_ft("${work}/sunrise_above.json")
if(ft_status EQUAL 0 OR NOT ft_output MATCHES "causal side of the real singular point x = 0\\.[0-9]+ of level 1 could not be derived"
    OR NOT ft_output MATCHES "explicit `causal` object that states its side" OR NOT ft_output MATCHES "points")
  message(FATAL_ERROR "undetermined two-loop contour was not refused with the documented message: ${ft_output} ${ft_errors}")
endif()
# Below threshold (but not Euclidean) the same level is free of cuts and real.
string(JSON below SET "${sunrise}" external_gram "[[\"5\"]]")
file(WRITE "${work}/sunrise_below.json" "${below}")
run_ft("${work}/sunrise_below.json")
expect("sunrise below threshold" 0 "^\\[3\\.27180532142958952052" "^(0|\\[\\+/-|\\[-?[0-9.]+e-[0-9][0-9])")
expect_mode("sunrise below threshold" derived)
file(REMOVE_RECURSE "${work}")
