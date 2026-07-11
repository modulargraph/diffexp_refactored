find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)
find_package(Boost 1.80 REQUIRED COMPONENTS json)

add_library(diffexp2_core STATIC
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
endif()
