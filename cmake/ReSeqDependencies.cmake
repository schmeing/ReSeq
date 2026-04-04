include(FetchContent)

# --- Options ---
option(RESEQ_ALLOW_FETCHCONTENT "Allow downloading missing dependencies at configure time" ON)
option(RESEQ_USE_VENDORED_SEQAN "Use the in-tree seqan/ fallback (deprecated)" OFF)

# --- Compression (required, system-provided) ---
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

# --- Boost (required, system-provided) ---
cmake_policy(SET CMP0057 NEW)
find_package(Boost 1.48.0 REQUIRED
  filesystem iostreams math_c99 math_c99f math_c99l
  math_tr1 math_tr1f math_tr1l program_options serialization system)

# --- GoogleTest (dev-only, always FetchContent) ---
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# --- NLopt (find_package first, FetchContent fallback) ---
find_package(NLopt QUIET)
if(NOT NLopt_FOUND)
  FetchContent_Declare(
    nlopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG v2.9.1
    SYSTEM
  )
  set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
  set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
  set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
  set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
  set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(nlopt)
endif()

# --- SeqAn (find_package first, FetchContent fallback, vendored last resort) ---
set(_reseq_seqan_resolved FALSE)

# 1. Try system/conda SeqAn with modern imported target
if(NOT _reseq_seqan_resolved)
  find_package(SeqAn 2.5 CONFIG QUIET)
  if(TARGET seqan::seqan2)
    add_library(reseq_seqan INTERFACE)
    target_link_libraries(reseq_seqan INTERFACE seqan::seqan2)
    set(_reseq_seqan_resolved TRUE)
    message(STATUS "SeqAn: using system package (seqan::seqan2)")
  endif()
endif()

# 2. Try legacy find_package (module mode)
if(NOT _reseq_seqan_resolved)
  find_package(SeqAn 2.5 QUIET)
  if(TARGET seqan::seqan2)
    add_library(reseq_seqan INTERFACE)
    target_link_libraries(reseq_seqan INTERFACE seqan::seqan2)
    set(_reseq_seqan_resolved TRUE)
    message(STATUS "SeqAn: using system package (legacy find_package)")
  endif()
endif()

# 3. FetchContent from upstream GitHub (header-only — download only, do not configure)
if(NOT _reseq_seqan_resolved AND RESEQ_ALLOW_FETCHCONTENT)
  message(STATUS "SeqAn: not found locally, fetching v2.5.2 via FetchContent...")
  FetchContent_Declare(
    seqan2
    GIT_REPOSITORY https://github.com/seqan/seqan.git
    GIT_TAG seqan-v2.5.2
    GIT_SHALLOW TRUE
  )
  # SeqAn2 is header-only. Use Populate (not MakeAvailable) to avoid running
  # SeqAn's own CMakeLists.txt which has internal build-system dependencies.
  FetchContent_GetProperties(seqan2)
  if(NOT seqan2_POPULATED)
    FetchContent_Populate(seqan2)
  endif()
  if(EXISTS "${seqan2_SOURCE_DIR}/include")
    add_library(reseq_seqan INTERFACE)
    target_include_directories(reseq_seqan SYSTEM INTERFACE "${seqan2_SOURCE_DIR}/include")
    target_compile_definitions(reseq_seqan INTERFACE SEQAN_HAS_ZLIB=1 SEQAN_HAS_BZIP2=1)
    target_link_libraries(reseq_seqan INTERFACE ZLIB::ZLIB BZip2::BZip2)
    set(_reseq_seqan_resolved TRUE)
    message(STATUS "SeqAn: using FetchContent (v2.5.2)")
  endif()
endif()

# 4. Vendored fallback (deprecated, for transition only)
if(NOT _reseq_seqan_resolved AND RESEQ_USE_VENDORED_SEQAN)
  if(EXISTS "${PROJECT_SOURCE_DIR}/seqan/include")
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake")
    find_package(SeqAn REQUIRED)
    add_library(reseq_seqan INTERFACE)
    target_include_directories(reseq_seqan SYSTEM INTERFACE
      "${PROJECT_SOURCE_DIR}/seqan/include"
      ${SEQAN_INCLUDE_DIRS}
    )
    target_compile_definitions(reseq_seqan INTERFACE ${SEQAN_DEFINITIONS})
    separate_arguments(SEQAN_CXX_FLAGS_LIST NATIVE_COMMAND "${SEQAN_CXX_FLAGS}")
    target_compile_options(reseq_seqan INTERFACE ${SEQAN_CXX_FLAGS_LIST})
    target_link_libraries(reseq_seqan INTERFACE ${SEQAN_LIBRARIES})
    set(_reseq_seqan_resolved TRUE)
    message(STATUS "SeqAn: using vendored fallback (deprecated)")
  endif()
endif()

if(NOT _reseq_seqan_resolved)
  message(FATAL_ERROR
    "SeqAn >= 2.5 not found. Options:\n"
    "  1. Enable RESEQ_ALLOW_FETCHCONTENT=ON (default) to download automatically\n"
    "  2. Install SeqAn and set CMAKE_PREFIX_PATH\n"
    "  3. Set RESEQ_USE_VENDORED_SEQAN=ON if seqan/ is still in-tree")
endif()
