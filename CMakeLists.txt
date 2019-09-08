cmake_minimum_required(VERSION 3.1) 
project(reseq)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Adding libraries from LD_LIBRARY_PATH
link_directories( $ENV{LD_LIBRARY_PATH} )

# Testing framework
add_subdirectory(googletest)
SET( GTEST_ROOT "${PROJECT_SOURCE_DIR}/googletest/googletest/" )
include_directories(${GTEST_ROOT}/include)

# Compression packages
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

# Add boost
find_package(Boost 1.48.0 REQUIRED filesystem program_options serialization system)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
SET( BOOST_LIBRARY_FLAGS "-lboost_serialization -lboost_program_options -lboost_filesystem -lboost_system" )
# Suppress warnings of deprecated std::auto_ptr used in boost
SET( BOOST_LIBRARY_FLAGS "${BOOST_LIBRARY_FLAGS} -Wno-deprecated-declarations" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${BOOST_LIBRARY_FLAGS}" )

# Add NLopt
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wdeprecated-declarations" )
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/nlopt/cmake/")
add_subdirectory(nlopt)
include_directories("${PROJECT_SOURCE_DIR}/build/nlopt/src/api/")

# Add SeqAn
SET(SEQAN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/seqan/include")
include_directories(${SEQAN_INCLUDE_PATH})

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake/")
find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Version
set(RESEQ_VERSION_MAJOR 0)
set(RESEQ_VERSION_MINOR 9)

# Configure header file to pass CMake settings to code
configure_file(
	"${PROJECT_SOURCE_DIR}/CMakeConfig.h.in"
	"${PROJECT_BINARY_DIR}/CMakeConfig.h"
)

# Compiler flags
SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++14 -fPIC")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

# Include directories
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/reseq")
include_directories("${PROJECT_SOURCE_DIR}/2016-05-15_ROOTPWA/utilities")

# Subdirectories
add_subdirectory(skewer)
add_subdirectory(reseq)
add_subdirectory(python)
