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

  add_executable(test_diffexp2_local_solution tests/test_local_solution.cpp)
  target_compile_features(test_diffexp2_local_solution PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_local_solution PRIVATE diffexp2_core)
  add_test(NAME cpp_local_solution COMMAND test_diffexp2_local_solution)

  add_executable(test_diffexp2_integration tests/test_integration.cpp)
  target_compile_features(test_diffexp2_integration PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_integration PRIVATE diffexp2_core)
  add_test(NAME cpp_integration COMMAND test_diffexp2_integration)

  add_executable(test_diffexp2_line_integration
    tests/test_line_integration.cpp)
  target_compile_features(test_diffexp2_line_integration PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_line_integration PRIVATE diffexp2_core)
  add_test(NAME cpp_line_integration COMMAND test_diffexp2_line_integration)

  add_executable(test_diffexp2_persistent_local
    tests/test_persistent_local.cpp)
  target_compile_features(test_diffexp2_persistent_local PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_local PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_local COMMAND test_diffexp2_persistent_local)

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

  add_executable(test_diffexp2_persistent_scc_solve_singular_column
    tests/test_persistent_scc_solve_singular_column.cpp)
  target_compile_features(
    test_diffexp2_persistent_scc_solve_singular_column PRIVATE cxx_std_20)
  target_link_libraries(
    test_diffexp2_persistent_scc_solve_singular_column PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_scc_solve_singular_column
    COMMAND test_diffexp2_persistent_scc_solve_singular_column)

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

  add_executable(test_diffexp2_persistent_checkpoint
    tests/test_persistent_checkpoint.cpp)
  target_compile_features(test_diffexp2_persistent_checkpoint
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_checkpoint
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_checkpoint
    COMMAND test_diffexp2_persistent_checkpoint)

  add_executable(test_diffexp2_persistent_endpoint_limit
    tests/test_persistent_endpoint_limit.cpp)
  target_compile_features(test_diffexp2_persistent_endpoint_limit
    PRIVATE cxx_std_20)
  target_link_libraries(test_diffexp2_persistent_endpoint_limit
    PRIVATE diffexp2_core)
  add_test(NAME cpp_persistent_endpoint_limit
    COMMAND test_diffexp2_persistent_endpoint_limit)
endif()
