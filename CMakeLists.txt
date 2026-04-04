cmake_minimum_required(VERSION 3.16)

# Read version from VERSION file (single source of truth)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" RESEQ_VERSION_RAW)
string(STRIP "${RESEQ_VERSION_RAW}" RESEQ_VERSION)
string(REPLACE "." ";" RESEQ_VERSION_LIST ${RESEQ_VERSION})
list(GET RESEQ_VERSION_LIST 0 RESEQ_VERSION_MAJOR)
list(GET RESEQ_VERSION_LIST 1 RESEQ_VERSION_MINOR)
list(GET RESEQ_VERSION_LIST 2 RESEQ_VERSION_PATCH)

project(reseq LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Require C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Export compile commands for clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Centralized dependency resolution
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(ReSeqDependencies)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Version (read from VERSION file above)
# Optional: git describe for dev builds
find_package(Git QUIET)
set(RESEQ_GIT_DESCRIBE "v${RESEQ_VERSION}")
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE RESEQ_GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE GIT_RESULT
  )
  if(NOT GIT_RESULT EQUAL 0)
    set(RESEQ_GIT_DESCRIBE "v${RESEQ_VERSION}")
  endif()
endif()

# Configure header file to pass CMake settings to code
configure_file(
	"${PROJECT_SOURCE_DIR}/CMakeConfig.h.in"
	"${PROJECT_BINARY_DIR}/CMakeConfig.h"
)

# Build options (declared before add_subdirectory so they are available)
option(RESEQ_BUILD_PYTHON "Build Python bindings (requires SWIG 3+ and python3-dev)" OFF)
option(CODE_COVERAGE "Enable code coverage instrumentation (GCC only)" OFF)

# Enable CTest (must be in root CMakeLists.txt for ctest to discover tests)
enable_testing()

# Subdirectories
add_subdirectory(skewer)
add_subdirectory(reseq)
if(RESEQ_BUILD_PYTHON)
  add_subdirectory(python)
endif()
