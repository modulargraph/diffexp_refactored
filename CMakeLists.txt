cmake_minimum_required(VERSION 3.24)
project(DiffExp VERSION 2.1.1 LANGUAGES CXX)
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
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(diffexp_core INTERFACE PkgConfig::FLINT Boost::json Threads::Threads)

option(DIFFEXP_BUILD_KERNEL_RUNTIME "Build the extracted prepared-input compatibility runtime" OFF)
set(diffexp_install_targets diffexp_core diffexp)
if(DIFFEXP_BUILD_KERNEL_RUNTIME)
add_library(diffexp_kernel STATIC
  src/kernel/checkpoint.cpp
  src/kernel/json_codec_backend_info.cpp
  src/kernel/json_codec.cpp
  src/kernel/json_codec_recurrence.cpp
  src/kernel/json_codec_scc.cpp
  src/kernel/json_codec_support.cpp)
add_library(DiffExp::kernel ALIAS diffexp_kernel)
set_target_properties(diffexp_kernel PROPERTIES EXPORT_NAME kernel
  POSITION_INDEPENDENT_CODE ON)
target_link_libraries(diffexp_kernel PUBLIC DiffExp::core)
list(APPEND diffexp_install_targets diffexp_kernel)
endif()

add_executable(diffexp src/main.cpp)
target_link_libraries(diffexp PRIVATE DiffExp::core)
if(DIFFEXP_BUILD_KERNEL_RUNTIME)
  target_link_libraries(diffexp PRIVATE DiffExp::kernel)
  target_compile_definitions(diffexp PRIVATE DIFFEXP_WITH_KERNEL_RUNTIME=1)
endif()

if(BUILD_TESTING)
  set(native_tests transport frontend data_geometry tadpole_bubble ibp sunrise_ibp
    deepest_beta sunrise_pipeline frobenius certified_tail certified_bubble
    banana_oracle henn_observables certified_deepest_beta merge_pullback fuchsify affine_frobenius
    epsilon_gauge laurent_transport affine_matching direct_examples singular_transport)
  if(UNIX)
    list(APPEND native_tests family_config fire_native level_preparation recursion_graph artifact_store certified_rational_transport
      gaussian_boundary level_cache scalar_functional recursion_leaf fire_session recursion_pipeline fire_batch_cache)
    list(APPEND native_tests box_oracle henn_boundary causal causal_contour pentagon_oracle affine_operator adjoint_transport linear_boundary
      polynomial_transport henn_scalar_leaf box_bubble_oracle double_box_oracle box_triangle_oracle fire_checkpoint adjoint_scaling
      endpoint_gauge scalar_accumulator affine_large fire_backward_checkpoint affine_finite_lag conditioned_adjoint cached_affine
      factored_transport factored_pipeline pipeline_preflight centered_adjoint prepared_stage adjoint_checkpoint modular_reconstruction fire_modular univariate_rational cached_affine_columns affine_endpoint_domain)
    # Explicit acceptance runner: large FT examples are opt-in, with a caller
    # supplied durable cache and finite preparation budget.
    add_executable(diffexp_ft_examples tests/test_ft_examples.cpp)
    target_link_libraries(diffexp_ft_examples PRIVATE DiffExp::core)
    add_executable(diffexp_henn_ft tools/check_henn_ft.cpp)
    target_link_libraries(diffexp_henn_ft PRIVATE DiffExp::core)
    add_executable(diffexp_benchmark_affine_domain tools/benchmark_affine_domain.cpp)
    target_link_libraries(diffexp_benchmark_affine_domain PRIVATE DiffExp::core)
    add_executable(diffexp_verify_affine_cache tools/verify_affine_cache.cpp)
    target_link_libraries(diffexp_verify_affine_cache PRIVATE DiffExp::core)
    add_executable(diffexp_diagnose_endpoint_sectors tools/diagnose_endpoint_sectors.cpp)
    target_link_libraries(diffexp_diagnose_endpoint_sectors PRIVATE DiffExp::core)
    add_executable(diffexp_affine_finite_lag_cached tests/test_affine_finite_lag_cached.cpp)
    target_link_libraries(diffexp_affine_finite_lag_cached PRIVATE DiffExp::core)
    add_executable(diffexp_probe_recursion_stage tools/probe_recursion_stage.cpp)
    target_link_libraries(diffexp_probe_recursion_stage PRIVATE DiffExp::core)
    add_executable(diffexp_replay_adjoint_stage tools/replay_adjoint_stage.cpp)
    target_link_libraries(diffexp_replay_adjoint_stage PRIVATE DiffExp::core)
    add_executable(diffexp_compose_recursion_stage tools/compose_recursion_stage.cpp)
    target_link_libraries(diffexp_compose_recursion_stage PRIVATE DiffExp::core)
  endif()
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
  if(DIFFEXP_BUILD_KERNEL_RUNTIME)
  file(GLOB kernel_tests CONFIGURE_DEPENDS tests/kernel/test_*.cpp)
  foreach(source IN LISTS kernel_tests)
    get_filename_component(name ${source} NAME_WE)
    add_executable(diffexp_${name} ${source})
    target_link_libraries(diffexp_${name} PRIVATE DiffExp::kernel)
    # Several inherited regressions use assert(); keep them active in Release.
    target_compile_options(diffexp_${name} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-UNDEBUG>
      $<$<CXX_COMPILER_ID:MSVC>:/UNDEBUG>)
    add_test(NAME kernel.${name} COMMAND diffexp_${name})
    set_tests_properties(kernel.${name} PROPERTIES TIMEOUT 120)
  endforeach()
  endif()
endif()

install(TARGETS ${diffexp_install_targets} EXPORT DiffExpTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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
    set_tests_properties(mathematica.interface mathematica.compatibility PROPERTIES TIMEOUT 120 RESOURCE_LOCK wolfram_kernel)
  endif()
endif()
