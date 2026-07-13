find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)
find_package(Boost 1.80 REQUIRED COMPONENTS json)

add_library(diffexp2_core STATIC
  src/checkpoint.cpp
  src/json_codec.cpp
)
set_target_properties(diffexp2_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(diffexp2_core PUBLIC cxx_std_20)
target_include_directories(diffexp2_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(diffexp2_core PUBLIC PkgConfig::FLINT Boost::json)
target_compile_options(diffexp2_core PRIVATE
  $<$<CONFIG:Release>:-O3>)

set(WOLFRAM_LIBRARYLINK_INCLUDE
  "${WOLFRAM_INSTALL_DIR}/SystemFiles/IncludeFiles/C")
if(NOT EXISTS "${WOLFRAM_LIBRARYLINK_INCLUDE}/WolframLibrary.h")
  message(FATAL_ERROR
    "WolframLibrary.h not found below ${WOLFRAM_INSTALL_DIR}; set WOLFRAM_INSTALL_DIR")
endif()

add_library(diffexp2_librarylink SHARED src/librarylink.cpp)
target_compile_features(diffexp2_librarylink PRIVATE cxx_std_20)
target_include_directories(diffexp2_librarylink PRIVATE
  ${WOLFRAM_LIBRARYLINK_INCLUDE})
target_link_libraries(diffexp2_librarylink PRIVATE diffexp2_core)
set_target_properties(diffexp2_librarylink PROPERTIES
  PREFIX ""
  OUTPUT_NAME "diffexp2_librarylink")

if(DE2_BUILD_TESTS)
  enable_testing()
  add_executable(test_diffexp2_recurrence tests/test_recurrence.cpp)
  target_compile_features(test_diffexp2_recurrence PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_recurrence PRIVATE diffexp2_core)
  add_test(NAME cpp_recurrence COMMAND test_diffexp2_recurrence)

  add_executable(test_diffexp2_singular_indicial
    tests/test_singular_indicial.cpp)
  target_compile_features(test_diffexp2_singular_indicial PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_singular_indicial PRIVATE diffexp2_core)
  add_test(NAME cpp_singular_indicial COMMAND test_diffexp2_singular_indicial)

  add_executable(test_diffexp2_local_solution tests/test_local_solution.cpp)
  target_compile_features(test_diffexp2_local_solution PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_local_solution PRIVATE diffexp2_core)
  add_test(NAME cpp_local_solution COMMAND test_diffexp2_local_solution)

  add_executable(test_diffexp2_integration tests/test_integration.cpp)
  target_compile_features(test_diffexp2_integration PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_integration PRIVATE diffexp2_core)
  add_test(NAME cpp_integration COMMAND test_diffexp2_integration)

  add_executable(test_diffexp2_scalar_row_endpoint
    tests/test_scalar_row_endpoint.cpp)
  target_compile_features(test_diffexp2_scalar_row_endpoint PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_scalar_row_endpoint PRIVATE diffexp2_core)
  add_test(NAME cpp_scalar_row_endpoint
    COMMAND test_diffexp2_scalar_row_endpoint)

  add_executable(test_diffexp2_line_integration
    tests/test_line_integration.cpp)
  target_compile_features(test_diffexp2_line_integration PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_line_integration PRIVATE diffexp2_core)
  add_test(NAME cpp_line_integration COMMAND test_diffexp2_line_integration)

  add_executable(test_diffexp2_tail_majorant
    tests/test_tail_majorant.cpp)
  target_compile_features(test_diffexp2_tail_majorant PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_tail_majorant PRIVATE diffexp2_core)
  add_test(NAME cpp_tail_majorant COMMAND test_diffexp2_tail_majorant)

  add_executable(test_diffexp2_persistent_tail_majorant
    tests/test_persistent_tail_majorant.cpp)
  target_compile_features(test_diffexp2_persistent_tail_majorant
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_tail_majorant
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_tail_majorant
    COMMAND test_diffexp2_persistent_tail_majorant)

  add_executable(test_diffexp2_persistent_rational_row_tail
    tests/test_persistent_rational_row_tail.cpp)
  target_compile_features(test_diffexp2_persistent_rational_row_tail
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_rational_row_tail
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_rational_row_tail
    COMMAND test_diffexp2_persistent_rational_row_tail)

  add_executable(test_diffexp2_persistent_matched_tail_majorant
    tests/test_persistent_matched_tail_majorant.cpp)
  target_compile_features(test_diffexp2_persistent_matched_tail_majorant
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_matched_tail_majorant
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_matched_tail_majorant
    COMMAND test_diffexp2_persistent_matched_tail_majorant)

  add_executable(test_diffexp2_persistent_local
    tests/test_persistent_local.cpp)
  target_compile_features(test_diffexp2_persistent_local PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_local PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_local COMMAND test_diffexp2_persistent_local)

  add_executable(test_diffexp2_persistent_owner_bound_residual
    tests/test_persistent_owner_bound_residual.cpp)
  target_compile_features(test_diffexp2_persistent_owner_bound_residual
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_owner_bound_residual
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_owner_bound_residual
    COMMAND test_diffexp2_persistent_owner_bound_residual)

  add_executable(test_diffexp2_persistent_scc_owner_bound_residual
    tests/test_persistent_scc_owner_bound_residual.cpp)
  target_compile_features(test_diffexp2_persistent_scc_owner_bound_residual
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_owner_bound_residual
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_owner_bound_residual
    COMMAND test_diffexp2_persistent_scc_owner_bound_residual)

  add_executable(test_diffexp2_session_solve_many
    tests/test_session_solve_many.cpp)
  target_compile_features(test_diffexp2_session_solve_many PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_session_solve_many PRIVATE diffexp2_core)
  add_test(NAME cpp_session_solve_many COMMAND test_diffexp2_session_solve_many)

  add_executable(test_diffexp2_local_algebra
    tests/test_local_algebra.cpp)
  target_compile_features(test_diffexp2_local_algebra PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_local_algebra PRIVATE diffexp2_core)
  add_test(NAME cpp_local_algebra COMMAND test_diffexp2_local_algebra)

  add_executable(test_diffexp2_persistent_scc_prepare
    tests/test_persistent_scc_prepare.cpp)
  target_compile_features(test_diffexp2_persistent_scc_prepare PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_prepare PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_prepare
    COMMAND test_diffexp2_persistent_scc_prepare)

  add_executable(test_diffexp2_persistent_scc_solve_column
    tests/test_persistent_scc_solve_column.cpp)
  target_compile_features(test_diffexp2_persistent_scc_solve_column
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_solve_column
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_solve_column
    COMMAND test_diffexp2_persistent_scc_solve_column)

  add_executable(test_diffexp2_persistent_scc_solve_block_column
    tests/test_persistent_scc_solve_block_column.cpp)
  target_compile_features(test_diffexp2_persistent_scc_solve_block_column
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_solve_block_column
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_solve_block_column
    COMMAND test_diffexp2_persistent_scc_solve_block_column)

  add_executable(test_diffexp2_persistent_scc_spectral_frames
    tests/test_persistent_scc_spectral_frames.cpp)
  target_compile_features(test_diffexp2_persistent_scc_spectral_frames
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_spectral_frames
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_spectral_frames
    COMMAND test_diffexp2_persistent_scc_spectral_frames)

  add_executable(test_diffexp2_persistent_scc_solve_singular_column
    tests/test_persistent_scc_solve_singular_column.cpp)
  target_compile_features(
    test_diffexp2_persistent_scc_solve_singular_column PRIVATE cxx_std_20)
  target_link_libraries(
    test_diffexp2_persistent_scc_solve_singular_column PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_solve_singular_column
    COMMAND test_diffexp2_persistent_scc_solve_singular_column)

  add_executable(test_diffexp2_persistent_scc_pseudo_compensation
    tests/test_persistent_scc_pseudo_compensation.cpp)
  target_compile_features(
    test_diffexp2_persistent_scc_pseudo_compensation PRIVATE cxx_std_20)
  target_link_libraries(
    test_diffexp2_persistent_scc_pseudo_compensation PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_pseudo_compensation
    COMMAND test_diffexp2_persistent_scc_pseudo_compensation)

  add_executable(test_diffexp2_matching tests/test_matching.cpp)
  target_compile_features(test_diffexp2_matching PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_matching PRIVATE diffexp2_core)
  add_test(NAME cpp_matching COMMAND test_diffexp2_matching)

  add_executable(test_diffexp2_persistent_local_match
    tests/test_persistent_local_match.cpp)
  target_compile_features(test_diffexp2_persistent_local_match
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_local_match
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_local_match
    COMMAND test_diffexp2_persistent_local_match)

  add_executable(test_diffexp2_persistent_acb_local_match
    tests/test_persistent_acb_local_match.cpp)
  target_compile_features(test_diffexp2_persistent_acb_local_match
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_acb_local_match
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_acb_local_match
    COMMAND test_diffexp2_persistent_acb_local_match)

  add_executable(test_diffexp2_persistent_checkpoint
    tests/test_persistent_checkpoint.cpp)
  target_compile_features(test_diffexp2_persistent_checkpoint
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_checkpoint
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_checkpoint
    COMMAND test_diffexp2_persistent_checkpoint)

  add_executable(test_diffexp2_persistent_checkpoint_deferred_state
    tests/test_persistent_checkpoint_deferred_state.cpp)
  target_compile_features(
    test_diffexp2_persistent_checkpoint_deferred_state PRIVATE cxx_std_20)
  target_link_libraries(
    test_diffexp2_persistent_checkpoint_deferred_state PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_checkpoint_deferred_state
    COMMAND test_diffexp2_persistent_checkpoint_deferred_state)

  add_executable(test_diffexp2_persistent_endpoint_limit
    tests/test_persistent_endpoint_limit.cpp)
  target_compile_features(test_diffexp2_persistent_endpoint_limit
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_endpoint_limit
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_endpoint_limit
    COMMAND test_diffexp2_persistent_endpoint_limit)

  add_executable(test_diffexp2_persistent_planned_endpoint
    tests/test_persistent_planned_endpoint.cpp)
  target_compile_features(test_diffexp2_persistent_planned_endpoint
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_planned_endpoint
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_planned_endpoint
    COMMAND test_diffexp2_persistent_planned_endpoint)

  add_executable(test_diffexp2_persistent_tile_line
    tests/test_persistent_tile_line.cpp)
  target_compile_features(test_diffexp2_persistent_tile_line
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_tile_line
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_tile_line
    COMMAND test_diffexp2_persistent_tile_line)

  add_executable(test_diffexp2_persistent_single_arm_tile_plan
    tests/test_persistent_single_arm_tile_plan.cpp)
  target_compile_features(test_diffexp2_persistent_single_arm_tile_plan
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_single_arm_tile_plan
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_single_arm_tile_plan
    COMMAND test_diffexp2_persistent_single_arm_tile_plan)

  add_executable(test_diffexp2_persistent_analytic_prescription
    tests/test_persistent_analytic_prescription.cpp)
  target_compile_features(test_diffexp2_persistent_analytic_prescription
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_analytic_prescription
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_analytic_prescription
    COMMAND test_diffexp2_persistent_analytic_prescription)

  add_executable(test_diffexp2_persistent_phase_planned_match
    tests/test_persistent_phase_planned_match.cpp)
  target_compile_features(test_diffexp2_persistent_phase_planned_match
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_phase_planned_match
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_phase_planned_match
    COMMAND test_diffexp2_persistent_phase_planned_match)

  add_executable(test_diffexp2_persistent_plan_match_hop
    tests/test_persistent_plan_match_hop.cpp)
  target_compile_features(test_diffexp2_persistent_plan_match_hop
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_plan_match_hop
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_plan_match_hop
    COMMAND test_diffexp2_persistent_plan_match_hop)

  add_executable(test_diffexp2_persistent_scc_plan_match
    tests/test_persistent_scc_plan_match.cpp)
  target_compile_features(test_diffexp2_persistent_scc_plan_match
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_scc_plan_match
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_plan_match
    COMMAND test_diffexp2_persistent_scc_plan_match)

  add_executable(test_diffexp2_persistent_parallel_arms
    tests/test_persistent_parallel_arms.cpp)
  target_compile_features(test_diffexp2_persistent_parallel_arms
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_parallel_arms
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_parallel_arms
    COMMAND test_diffexp2_persistent_parallel_arms)

  add_executable(test_diffexp2_persistent_transport_run_arms
    tests/test_persistent_transport_run_arms.cpp)
  target_compile_features(test_diffexp2_persistent_transport_run_arms
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_transport_run_arms
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_transport_run_arms
    COMMAND test_diffexp2_persistent_transport_run_arms)

  add_executable(test_diffexp2_persistent_singular_scc_arms
    tests/test_persistent_singular_scc_arms.cpp)
  target_compile_features(test_diffexp2_persistent_singular_scc_arms
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_singular_scc_arms
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_singular_scc_arms
    COMMAND test_diffexp2_persistent_singular_scc_arms)

  add_executable(test_diffexp2_persistent_transport_contract_pair
    tests/test_persistent_transport_contract_pair.cpp)
  target_compile_features(test_diffexp2_persistent_transport_contract_pair
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_transport_contract_pair
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_transport_contract_pair
    COMMAND test_diffexp2_persistent_transport_contract_pair)

  add_executable(test_diffexp2_persistent_transport_endpoint_batch
    tests/test_persistent_transport_endpoint_batch.cpp)
  target_compile_features(test_diffexp2_persistent_transport_endpoint_batch
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_transport_endpoint_batch
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_transport_endpoint_batch
    COMMAND test_diffexp2_persistent_transport_endpoint_batch)

  add_executable(test_diffexp2_persistent_rational_row
    tests/test_persistent_rational_row.cpp)
  target_compile_features(test_diffexp2_persistent_rational_row
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_rational_row
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_rational_row
    COMMAND test_diffexp2_persistent_rational_row)

  add_executable(test_diffexp2_path_planner tests/test_path_planner.cpp)
  target_compile_features(test_diffexp2_path_planner PRIVATE cxx_std_20)
  target_include_directories(test_diffexp2_path_planner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(test_diffexp2_path_planner PRIVATE PkgConfig::FLINT)
  add_test(NAME cpp_path_planner COMMAND test_diffexp2_path_planner)
endif()
