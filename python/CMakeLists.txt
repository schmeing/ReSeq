# This subdirectory is only added when RESEQ_BUILD_PYTHON=ON.
# Fail clearly if prerequisites are missing rather than silently no-op.

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(SWIG 3 REQUIRED)
include(UseSWIG)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

set(UseSWIG_MODULE_VERSION 2)

set(PYTHON_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/pyMods")
set(PYTHON_SWIG_OUTFILE_DIR "${PROJECT_BINARY_DIR}/python/CMakeFiles/DataStats")
file(MAKE_DIRECTORY "${PYTHON_MODULE_DIRECTORY}")
file(MAKE_DIRECTORY "${PYTHON_SWIG_OUTFILE_DIR}")

set_property(
  SOURCE DataStats.i
  PROPERTY CPLUSPLUS ON
)
set_property(
  SOURCE DataStats.i
  PROPERTY SWIG_MODULE_NAME DataStats
)
set_property(
  SOURCE DataStats.i
  PROPERTY USE_SWIG_DEPENDENCIES TRUE
)
set_property(
  SOURCE DataStats.i
  PROPERTY DEPENDS
           "${PROJECT_SOURCE_DIR}/reseq/DataStatsInterface.h"
           "${PROJECT_SOURCE_DIR}/reseq/Reference.h"
           "${PROJECT_SOURCE_DIR}/reseq/Surrounding.h"
           "${PROJECT_SOURCE_DIR}/reseq/utilities.hpp"
           "${PROJECT_SOURCE_DIR}/reseq/Vect.hpp"
)

swig_add_library(
  DataStats
  TYPE MODULE
  LANGUAGE python
  OUTPUT_DIR "${PYTHON_MODULE_DIRECTORY}"
  OUTFILE_DIR "${PYTHON_SWIG_OUTFILE_DIR}"
  SOURCES DataStats.i
)

target_include_directories(
  DataStats
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/python"
    "${PROJECT_SOURCE_DIR}/reseq"
)
set_target_properties(
  DataStats
  PROPERTIES
    SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE
    SWIG_COMPILE_OPTIONS "-Wall;-modern;-DSWIGWORDSIZE64"
    LIBRARY_OUTPUT_DIRECTORY "${PYTHON_MODULE_DIRECTORY}"
)
target_link_libraries(
  DataStats
  PRIVATE
    Python3::Module
    reseq_lib
)

install(TARGETS DataStats DESTINATION lib)
install(FILES "${PYTHON_MODULE_DIRECTORY}/DataStats.py" DESTINATION lib)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/python/reseq/"
  DESTINATION lib/python/reseq
  FILES_MATCHING PATTERN "*.py"
)

add_test(
  NAME python_script_tests
  COMMAND "${Python3_EXECUTABLE}" -m unittest discover -s "${PROJECT_SOURCE_DIR}/python/tests" -v
)
set_tests_properties(
  python_script_tests
  PROPERTIES
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
)

add_test(
  NAME python_datastats_import
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/python/tests/import_datastats.py"
)
set_tests_properties(
  python_datastats_import
  PROPERTIES
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    ENVIRONMENT "PYTHONPATH=${PYTHON_MODULE_DIRECTORY}"
)
