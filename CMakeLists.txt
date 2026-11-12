cmake_minimum_required(VERSION 3.24)
project(DiffExp VERSION 2.2.1 LANGUAGES CXX)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint>=3.4)
find_package(Boost 1.80 REQUIRED COMPONENTS json)
find_package(Threads REQUIRED)

add_library(diffexp_core INTERFACE)
add_library(DiffExp::core ALIAS diffexp_core)
set_target_properties(diffexp_core PROPERTIES EXPORT_NAME core)
target_compile_features(diffexp_core INTERFACE cxx_std_20)
target_include_directories(diffexp_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/ibp-solver/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(diffexp_core INTERFACE PkgConfig::FLINT Boost::json Threads::Threads)
# The executables built here are compiled for the CPU of the build machine
# (the library is header-only, so consumers choose their own flags).
# Floating-point contraction stays off: the radii of the compact recurrence
# are doubles and must not depend on whether the CPU has fused multiply-add.
option(DIFFEXP_NATIVE "Compile this build's executables for the host CPU" ON)
if(DIFFEXP_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-march=native DIFFEXP_HAS_MARCH_NATIVE)
  check_cxx_compiler_flag(-mcpu=native DIFFEXP_HAS_MCPU_NATIVE)
  if(DIFFEXP_HAS_MARCH_NATIVE)
    add_compile_options(-march=native -ffp-contract=off)
  elseif(DIFFEXP_HAS_MCPU_NATIVE)
    add_compile_options(-mcpu=native -ffp-contract=off)
  endif()
endif()
# OpenMP is optional: it expands large Taylor charts on several threads. The
# results are identical with and without it.
option(DIFFEXP_OPENMP "Use OpenMP threads for large Taylor charts" ON)
set(DIFFEXP_HAS_OPENMP OFF)
if(DIFFEXP_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(diffexp_core INTERFACE OpenMP::OpenMP_CXX)
    set(DIFFEXP_HAS_OPENMP ON)
  endif()
endif()

set(diffexp_install_targets diffexp_core diffexp)

add_executable(diffexp src/main.cpp)
target_link_libraries(diffexp PRIVATE DiffExp::core)

if(BUILD_TESTING)
  set(native_tests transport spectral_transport finite_lag_transport frontend data_geometry tadpole_bubble ibp sunrise_ibp
    deepest_beta sunrise_pipeline frobenius certified_tail certified_bubble
    banana_oracle henn_observables certified_deepest_beta merge_pullback fuchsify affine_frobenius
    epsilon_gauge laurent_transport affine_matching direct_examples singular_transport)
  if(UNIX)
    list(APPEND native_tests ibp_solver_provider family_config fire_native level_preparation recursion_graph artifact_store certified_rational_transport
      gaussian_boundary level_cache scalar_functional recursion_leaf fire_session recursion_pipeline fire_batch_cache)
    list(APPEND native_tests box_oracle henn_boundary causal causal_contour pentagon_oracle affine_operator adjoint_transport linear_boundary
      henn_scalar_leaf box_bubble_oracle double_box_oracle box_triangle_oracle fire_checkpoint adjoint_scaling
      endpoint_gauge scalar_accumulator affine_large fire_backward_checkpoint affine_finite_lag conditioned_adjoint
      factored_transport factored_pipeline pipeline_preflight centered_adjoint prepared_stage adjoint_checkpoint modular_reconstruction fire_modular univariate_rational affine_endpoint_domain)
    # Explicit acceptance runner: large FT examples are opt-in, with a caller
    # supplied durable cache and finite preparation budget.
    add_executable(diffexp_check_ft_laptop tools/check_ft_laptop.cpp)
    target_link_libraries(diffexp_check_ft_laptop PRIVATE DiffExp::core)
    add_executable(diffexp_ft_examples tests/test_ft_examples.cpp)
    target_link_libraries(diffexp_ft_examples PRIVATE DiffExp::core)
    add_executable(diffexp_henn_ft tools/check_henn_ft.cpp)
    target_link_libraries(diffexp_henn_ft PRIVATE DiffExp::core)
    add_executable(diffexp_probe_recursion_stage tools/probe_recursion_stage.cpp)
    target_link_libraries(diffexp_probe_recursion_stage PRIVATE DiffExp::core)
    add_executable(diffexp_replay_adjoint_stage tools/replay_adjoint_stage.cpp)
    target_link_libraries(diffexp_replay_adjoint_stage PRIVATE DiffExp::core)
    add_executable(diffexp_compose_recursion_stage tools/compose_recursion_stage.cpp)
    target_link_libraries(diffexp_compose_recursion_stage PRIVATE DiffExp::core)
    add_executable(diffexp_master_basis_report tools/master_basis_report.cpp)
    target_link_libraries(diffexp_master_basis_report PRIVATE DiffExp::core)
  endif()
  list(APPEND native_tests rational_circuit rational_circuit_dot compact_recurrence compact_adjoint uncertainty_support basis_persistence
    component_demands physical_basis_candidates ft_functional_reduction exact_jet_composition basis_screening moving_poles finite_masters epsilon_sampling critical_points endpoint_boundary)
  foreach(name IN LISTS native_tests)
    add_executable(diffexp_test_${name} tests/test_${name}.cpp)
    target_link_libraries(diffexp_test_${name} PRIVATE DiffExp::core)
    add_test(NAME ${name} COMMAND diffexp_test_${name})
    set_tests_properties(${name} PROPERTIES TIMEOUT 120)
  endforeach()
  add_executable(diffexp_original_banana_unequal_routes tests/test_original_banana_unequal_routes.cpp)
  target_link_libraries(diffexp_original_banana_unequal_routes PRIVATE DiffExp::core)
  add_executable(diffexp_original_banana_32 tools/check_original_banana_32.cpp)
  target_link_libraries(diffexp_original_banana_32 PRIVATE DiffExp::core)
  add_test(NAME original_mpl COMMAND diffexp mpl)
  set_tests_properties(original_mpl PROPERTIES TIMEOUT 180 LABELS original)
  set(DIFFEXP_ORIGINAL_EXAMPLES_DIR "" CACHE PATH "Optional original-example data for full reference and route comparisons")
  if(DIFFEXP_ORIGINAL_EXAMPLES_DIR)
    add_test(NAME original_henn_nonplanar COMMAND diffexp henn-nonplanar
      "${DIFFEXP_ORIGINAL_EXAMPLES_DIR}/Data/HennNonplanar")
    set_tests_properties(original_henn_nonplanar PROPERTIES TIMEOUT 180 LABELS original)
    foreach(family 1loop zmz mzz zzz)
      add_test(NAME original_planar.${family} COMMAND diffexp planar
        "${DIFFEXP_ORIGINAL_EXAMPLES_DIR}/Data/PlanarOneMass" ${family})
      set_tests_properties(original_planar.${family} PROPERTIES TIMEOUT 600 LABELS original)
    endforeach()
    add_test(NAME singular_transport.original COMMAND diffexp_test_singular_transport "${DIFFEXP_ORIGINAL_EXAMPLES_DIR}")
    add_test(NAME original_banana_equal.historical32 COMMAND diffexp_original_banana_32 "${DIFFEXP_ORIGINAL_EXAMPLES_DIR}")
    set_tests_properties(original_banana_equal.historical32 PROPERTIES TIMEOUT 120)
    add_test(NAME original_banana_unequal.routes COMMAND diffexp_original_banana_unequal_routes "${DIFFEXP_ORIGINAL_EXAMPLES_DIR}")
    set_tests_properties(singular_transport.original PROPERTIES TIMEOUT 120)
    set_tests_properties(original_banana_unequal.routes PROPERTIES TIMEOUT 600)
    set_tests_properties(singular_transport.original original_banana_equal.historical32
      original_banana_unequal.routes PROPERTIES LABELS original)
  endif()
  set(DIFFEXP_FIRE_EXECUTABLE "" CACHE FILEPATH "Optional FIRE7 executable for cold integration tests")
  set(DIFFEXP_FIRE_PRIME_EXECUTABLE "" CACHE FILEPATH "Optional FIRE7p executable for modular reconstruction integration tests")
  if(UNIX AND DIFFEXP_FIRE_PRIME_EXECUTABLE)
    add_test(NAME fire_modular.external COMMAND diffexp_test_fire_modular "${DIFFEXP_FIRE_PRIME_EXECUTABLE}")
    set_tests_properties(fire_modular.external PROPERTIES TIMEOUT 180)
  endif()
  if(UNIX AND DIFFEXP_FIRE_EXECUTABLE)
    add_test(NAME ibp_solver_provider.external COMMAND diffexp_test_ibp_solver_provider "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(ibp_solver_provider.external PROPERTIES TIMEOUT 120)
    add_test(NAME family_config.external COMMAND diffexp_test_family_config "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(family_config.external PROPERTIES TIMEOUT 180)
    add_test(NAME fire_native.external COMMAND diffexp_test_fire_native "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fire_native.external PROPERTIES TIMEOUT 180)
    add_test(NAME level_preparation.external COMMAND diffexp_test_level_preparation "${DIFFEXP_FIRE_EXECUTABLE}" --banana4)
    set_tests_properties(level_preparation.external PROPERTIES TIMEOUT 180)
    add_test(NAME recursion_graph.external COMMAND diffexp_test_recursion_graph "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(recursion_graph.external PROPERTIES TIMEOUT 180)
    add_test(NAME fuchsify.external COMMAND diffexp_test_fuchsify "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fuchsify.external PROPERTIES TIMEOUT 180)
    add_test(NAME affine_frobenius.external COMMAND diffexp_test_affine_frobenius "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(affine_frobenius.external PROPERTIES TIMEOUT 180)
    add_test(NAME affine_matching.external COMMAND diffexp_test_affine_matching "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(affine_matching.external PROPERTIES TIMEOUT 180)
    add_test(NAME fire_session.external COMMAND diffexp_test_fire_session "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fire_session.external PROPERTIES TIMEOUT 180)
    add_test(NAME fire_batch_cache.external COMMAND diffexp_test_fire_batch_cache "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fire_batch_cache.external PROPERTIES TIMEOUT 180)
    add_test(NAME fire_checkpoint.external COMMAND diffexp_test_fire_checkpoint "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fire_checkpoint.external PROPERTIES TIMEOUT 180)
    add_test(NAME fire_backward_checkpoint.external COMMAND diffexp_test_fire_backward_checkpoint "${DIFFEXP_FIRE_EXECUTABLE}")
    set_tests_properties(fire_backward_checkpoint.external PROPERTIES TIMEOUT 180)
  endif()
  set(DIFFEXP_HENN_BASIS_FILE "" CACHE FILEPATH "Optional published dlogBasisXB.txt for all-108 import tests")
  if(DIFFEXP_HENN_BASIS_FILE)
    add_test(NAME henn_observables.ancillary COMMAND diffexp_test_henn_observables "${DIFFEXP_HENN_BASIS_FILE}")
    add_test(NAME merge_pullback.ancillary COMMAND diffexp_test_merge_pullback "${DIFFEXP_HENN_BASIS_FILE}")
    if(UNIX)
      add_test(NAME henn_boundary.ancillary COMMAND diffexp_test_henn_boundary "${DIFFEXP_HENN_BASIS_FILE}")
    endif()
  endif()
endif()

install(TARGETS ${diffexp_install_targets} EXPORT DiffExpTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY third_party/ibp-solver/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES third_party/ibp-solver/LICENSE third_party/ibp-solver/README.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/diffexp/third_party/ibp-solver)
install(EXPORT DiffExpTargets NAMESPACE DiffExp:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DiffExp)
configure_package_config_file(cmake/DiffExpConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/DiffExpConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DiffExp)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/DiffExpConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/DiffExpConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/DiffExpConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DiffExp)

install(DIRECTORY examples/feynman/ DESTINATION ${CMAKE_INSTALL_DATADIR}/diffexp/examples/feynman)
install(DIRECTORY examples/transport/ DESTINATION ${CMAKE_INSTALL_DATADIR}/diffexp/examples/transport)
install(DIRECTORY Mathematica/ DESTINATION ${CMAKE_INSTALL_DATADIR}/diffexp/Mathematica)
install(FILES LICENSE README.md DESTINATION ${CMAKE_INSTALL_DATADIR}/diffexp)
if(BUILD_TESTING)
  add_test(NAME cli.interface COMMAND ${CMAKE_COMMAND}
    -DDIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>
    -DDIFFEXP_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli.cmake)
  set_tests_properties(cli.interface PROPERTIES TIMEOUT 30)
  add_test(NAME cli.causal COMMAND ${CMAKE_COMMAND}
    -DDIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>
    -DDIFFEXP_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/causal.cmake)
  set_tests_properties(cli.causal PROPERTIES TIMEOUT 120)
  add_test(NAME cli.endpoint_boundary COMMAND ${CMAKE_COMMAND}
    -DDIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>
    -DDIFFEXP_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/endpoint_boundary.cmake)
  set_tests_properties(cli.endpoint_boundary PROPERTIES TIMEOUT 120)
  set(DIFFEXP_WOLFRAMSCRIPT "" CACHE FILEPATH "Optional wolframscript executable for the process-wrapper smoke test")
  set(DIFFEXP_WOLFRAM_KERNEL "" CACHE FILEPATH "Optional WolframKernel override for wolframscript")
  if(DIFFEXP_WOLFRAMSCRIPT)
    set(diffexp_wolfram_arguments)
    if(DIFFEXP_WOLFRAM_KERNEL)
      list(APPEND diffexp_wolfram_arguments -local "${DIFFEXP_WOLFRAM_KERNEL}")
    endif()
    add_test(NAME mathematica.interface COMMAND ${CMAKE_COMMAND} -E env
      "DIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>" ${DIFFEXP_WOLFRAMSCRIPT} ${diffexp_wolfram_arguments}
      -file ${CMAKE_CURRENT_SOURCE_DIR}/tests/mathematica/smoke.wl)
    add_test(NAME mathematica.compatibility COMMAND ${CMAKE_COMMAND} -E env
      "DIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>" ${DIFFEXP_WOLFRAMSCRIPT} ${diffexp_wolfram_arguments}
      -file ${CMAKE_CURRENT_SOURCE_DIR}/tests/mathematica/compatibility.wl)
    add_test(NAME mathematica.ibp_solver COMMAND ${CMAKE_COMMAND} -E env
      "DIFFEXP_EXECUTABLE=$<TARGET_FILE:diffexp>" ${DIFFEXP_WOLFRAMSCRIPT} ${diffexp_wolfram_arguments}
      -file ${CMAKE_CURRENT_SOURCE_DIR}/tests/mathematica/ibp-solver.wl)
    set_tests_properties(mathematica.interface mathematica.compatibility mathematica.ibp_solver PROPERTIES TIMEOUT 120 RESOURCE_LOCK wolfram_kernel)
  endif()
endif()
