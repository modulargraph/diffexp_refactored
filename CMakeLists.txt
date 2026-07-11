cmake_minimum_required(VERSION 3.24)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
project(DiffExp2Cpp LANGUAGES CXX)

option(DE2_BUILD_TESTS "Build the C++ recurrence tests" ON)
if(DE2_BUILD_TESTS)
  enable_testing()
endif()
set(WOLFRAM_INSTALL_DIR "/Applications/Wolfram.app/Contents" CACHE PATH
    "Wolfram application Contents directory")

add_subdirectory(cpp)
