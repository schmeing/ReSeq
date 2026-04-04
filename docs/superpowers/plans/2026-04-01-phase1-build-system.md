# Phase 1: Build System & Test Architecture — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize the CMake build system to C++20 / CMake 3.16+, extract a core library target, separate tests into a dedicated binary with CTest integration, and switch vendored dependencies to FetchContent where appropriate.

**Architecture:** Incremental migration — each task produces a buildable, testable tree. Global CMake state is replaced one dependency at a time (least entangled first), then a `reseq_lib` static library is extracted from the production sources, a dedicated `reseq_test` binary is created from the test sources, and finally the old `reseq test` subcommand is removed. CI is extended with sanitizer and coverage jobs.

**Tech Stack:** CMake 3.16+, C++20, GoogleTest (FetchContent), NLopt (FetchContent), CTest, gcov/lcov, GitHub Actions

**Spec:** [docs/superpowers/specs/2026-04-01-comprehensive-refactoring-design.md](../specs/2026-04-01-comprehensive-refactoring-design.md) — Phase 1

---

## File Map

| Action | File | Purpose |
|--------|------|---------|
| Modify | `CMakeLists.txt` | Root build: bump minimum, C++20, remove global state, FetchContent deps |
| Modify | `reseq/CMakeLists.txt` | Extract `reseq_lib`, create `reseq_test`, update `reseq` target |
| Modify | `CMakeConfig.h.in` | Add `RESEQ_BINARY_DIR` for test binary to locate `reseq` executable |
| Modify | `reseq/main.cpp` | Remove test registration, gtest include, test subcommand |
| Modify | `Makefile` | Update `test` target to use `ctest` |
| Modify | `.github/workflows/ci.yml` | Switch to ctest, add sanitizer + coverage jobs |
| Create | `reseq/test_main.cpp` | GoogleTest main for `reseq_test` binary |
| Delete | (vendored `googletest/` directory) | Replaced by FetchContent |

---

### Task 1: Capture Pre-Split Test Baseline

**Files:**
- None modified — verification only

Before any build changes, record the current test pass/fail state so we can verify the migration doesn't break anything.

- [ ] **Step 1: Build the current tree**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
```

- [ ] **Step 2: Run the existing test suite and capture output**

```bash
build/bin/reseq test 2>&1 | tee test-baseline.log
```

Expected: All 25 tests pass. Save the count of `[  PASSED  ]` lines.

- [ ] **Step 3: Record baseline checksum**

```bash
grep -c '\[  PASSED  \]' test-baseline.log > test-baseline-count.txt
cat test-baseline-count.txt
```

Expected: A number (e.g., `25` or the total test count). This file is for local verification only — do not commit it.

- [ ] **Step 4: Commit (no changes — just a checkpoint note)**

No files to commit. This is a verification-only task.

---

### Task 2: Bump CMake Minimum to 3.16 and Set C++20

**Files:**
- Modify: `CMakeLists.txt:1` (cmake_minimum_required)
- Modify: `CMakeLists.txt:21` (CMAKE_CXX_STANDARD)

- [ ] **Step 1: Update cmake_minimum_required**

In `CMakeLists.txt`, change line 1:

```cmake
cmake_minimum_required(VERSION 3.16)
```

- [ ] **Step 2: Update C++ standard**

In `CMakeLists.txt`, change line 21 from `set(CMAKE_CXX_STANDARD 14)` to:

```cmake
set(CMAKE_CXX_STANDARD 20)
```

- [ ] **Step 3: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass. If the compiler doesn't support C++20, the configure step will fail with a clear message (due to `CMAKE_CXX_STANDARD_REQUIRED ON` already being set).

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: bump CMake minimum to 3.16 and set C++20 standard"
```

---

### Task 3: Replace Vendored GoogleTest with FetchContent

**Files:**
- Modify: `CMakeLists.txt:28-30` (remove `add_subdirectory(googletest)` and related lines)
- Modify: `CMakeLists.txt` (add FetchContent block)

- [ ] **Step 1: Remove vendored GoogleTest lines from root CMakeLists.txt**

In `CMakeLists.txt`, remove these three lines (28-30):

```cmake
add_subdirectory(googletest)
SET( GTEST_ROOT "${PROJECT_SOURCE_DIR}/googletest/googletest/" )
include_directories(${GTEST_ROOT}/include)
```

Replace with FetchContent block (insert at the same location):

```cmake
# GoogleTest via FetchContent
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
```

- [ ] **Step 2: Update reseq/CMakeLists.txt to use GTest::gtest**

In `reseq/CMakeLists.txt`, change line 20 from:

```cmake
target_link_libraries(reseq gtest)
```

to:

```cmake
target_link_libraries(reseq GTest::gtest)
```

- [ ] **Step 3: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass. GoogleTest is now fetched from GitHub at configure time.

- [ ] **Step 4: Delete vendored googletest directory**

```bash
rm -rf googletest/
```

- [ ] **Step 5: Update .gitmodules if googletest was a submodule**

Check if `googletest` is listed in `.gitmodules`. If so, remove the entry. If `.gitmodules` becomes empty, delete it.

```bash
git ls-files .gitmodules
# If it exists and has a googletest entry:
git rm googletest
```

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt reseq/CMakeLists.txt
git rm -r googletest/ 2>/dev/null || true
git commit -m "build: replace vendored GoogleTest with FetchContent v1.15.2"
```

---

### Task 4: Replace Vendored NLopt with FetchContent

**Files:**
- Modify: `CMakeLists.txt:53-57` (NLopt section)

- [ ] **Step 1: Remove vendored NLopt lines from root CMakeLists.txt**

In `CMakeLists.txt`, remove these lines (53-57):

```cmake
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wdeprecated-declarations" )
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/nlopt/cmake/")
add_subdirectory(nlopt)
include_directories("${CMAKE_BINARY_DIR}/nlopt/src/api/")
```

Replace with FetchContent block:

```cmake
# NLopt via FetchContent (or system install)
find_package(NLopt QUIET)
if(NOT NLopt_FOUND)
  FetchContent_Declare(
    nlopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG v2.9.1
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
```

- [ ] **Step 2: Update reseq/CMakeLists.txt to use nlopt target**

In `reseq/CMakeLists.txt`, verify line 8 already says `target_link_libraries(DataStats nlopt)`. The FetchContent target name is `nlopt`, matching the existing usage. No change needed if it matches.

- [ ] **Step 3: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass.

- [ ] **Step 4: Delete vendored nlopt directory**

```bash
rm -rf nlopt/
```

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt
git rm -r nlopt/ 2>/dev/null || true
git commit -m "build: replace vendored NLopt with FetchContent v2.9.1"
```

---

### Task 5: Wrap SeqAn2 and Skewer as Proper CMake Targets

**Files:**
- Modify: `CMakeLists.txt:59-67` (SeqAn section)
- Modify: `CMakeLists.txt:109` (ROOTPWA include)
- Modify: `CMakeLists.txt:107-108` (project/reseq includes)

- [ ] **Step 1: Replace SeqAn global includes with an INTERFACE target**

In `CMakeLists.txt`, remove lines 59-67:

```cmake
# Add SeqAn
SET(SEQAN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/seqan/include")
include_directories(${SEQAN_INCLUDE_PATH})

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake/")
find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
```

Replace with:

```cmake
# SeqAn2 (vendored) — wrap as INTERFACE target with SYSTEM includes
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake/")
find_package(SeqAn REQUIRED)
add_library(vendored_seqan INTERFACE)
target_include_directories(vendored_seqan SYSTEM INTERFACE
  "${PROJECT_SOURCE_DIR}/seqan/include"
  ${SEQAN_INCLUDE_DIRS}
)
target_compile_definitions(vendored_seqan INTERFACE ${SEQAN_DEFINITIONS})
target_compile_options(vendored_seqan INTERFACE ${SEQAN_CXX_FLAGS})
target_link_libraries(vendored_seqan INTERFACE ${SEQAN_LIBRARIES})
```

- [ ] **Step 2: Wrap ROOTPWA as a temporary INTERFACE target**

In `CMakeLists.txt`, remove line 109:

```cmake
include_directories("${PROJECT_SOURCE_DIR}/2016-05-15_ROOTPWA/utilities")
```

Replace with:

```cmake
# ROOTPWA utilities (vendored, temporary — absorbed into reseq/ in Phase 3)
add_library(vendored_rootpwa INTERFACE)
target_include_directories(vendored_rootpwa SYSTEM INTERFACE
  "${PROJECT_SOURCE_DIR}/2016-05-15_ROOTPWA/utilities"
)
```

- [ ] **Step 3: Update reseq/CMakeLists.txt to use the new targets**

In `reseq/CMakeLists.txt`, update the DataStats target (lines 7-10):

Replace:
```cmake
target_link_libraries(DataStats skewer_matrix)
target_link_libraries(DataStats nlopt)
target_link_libraries(DataStats ${SEQAN_LIBRARIES})
target_link_libraries(DataStats ${Boost_LIBRARIES})
```

With:
```cmake
target_link_libraries(DataStats PRIVATE skewer_matrix nlopt vendored_seqan vendored_rootpwa ${Boost_LIBRARIES})
```

Update the reseq executable target (lines 19-22):

Replace:
```cmake
target_link_libraries(reseq DataStats)
target_link_libraries(reseq GTest::gtest)
target_link_libraries(reseq ${SEQAN_LIBRARIES})
target_link_libraries(reseq ${Boost_LIBRARIES})
```

With:
```cmake
target_link_libraries(reseq PRIVATE DataStats DataStatsInterface GTest::gtest vendored_seqan ${Boost_LIBRARIES})
```

- [ ] **Step 4: Remove remaining global include_directories**

In `CMakeLists.txt`, remove lines 106-108:

```cmake
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/reseq")
```

Replace with target-specific includes. In `reseq/CMakeLists.txt`, add to each target:

```cmake
target_include_directories(DataStats PRIVATE
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/reseq"
)
target_include_directories(DataStatsInterface PRIVATE
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/reseq"
)
target_include_directories(reseq PRIVATE
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/reseq"
)
```

- [ ] **Step 5: Add SYSTEM to Boost include**

In `CMakeLists.txt`, change the Boost include to use the target properly. Remove:

```cmake
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
```

Instead, ensure Boost includes are passed via target. Since `${Boost_LIBRARIES}` includes imported targets on modern CMake, the includes propagate automatically. If they don't, add `target_include_directories(DataStats SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})` to `reseq/CMakeLists.txt`.

- [ ] **Step 6: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass. No global `include_directories()` remain except possibly Boost fallback.

- [ ] **Step 7: Commit**

```bash
git add CMakeLists.txt reseq/CMakeLists.txt
git commit -m "build: wrap vendored deps as proper CMake targets, remove global includes"
```

---

### Task 6: Clean Up Global Compiler Flags

**Files:**
- Modify: `CMakeLists.txt` (remove CMAKE_CXX_FLAGS manipulation)
- Modify: `reseq/CMakeLists.txt` (add target-specific flags)

- [ ] **Step 1: Remove global CMAKE_CXX_FLAGS lines**

In `CMakeLists.txt`, remove all `SET(CMAKE_CXX_FLAGS ...)` lines (approximately lines 42-44, 54, 99-102). These were:

```cmake
SET( BOOST_LIBRARY_FLAGS "-lboost_serialization -lboost_program_options -lboost_filesystem -lboost_system" )
SET( BOOST_LIBRARY_FLAGS "${BOOST_LIBRARY_FLAGS} -Wno-deprecated-declarations" )
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${BOOST_LIBRARY_FLAGS}" )
...
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wdeprecated-declarations" )
...
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
endif()
```

Also remove the global `link_libraries(rt)` on line 103.

- [ ] **Step 2: Add target-specific flags in reseq/CMakeLists.txt**

Update the `PRIVATE_COMPILE_OPTIONS` variable and apply per-target:

```cmake
set(PRIVATE_COMPILE_OPTIONS -Wall -Wextra -Wno-sign-compare -Wno-implicit-fallthrough -fPIC)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND PRIVATE_COMPILE_OPTIONS -fext-numeric-literals -Wno-deprecated-declarations)
endif()
```

Add `rt` to the link libraries for targets that need it:

```cmake
target_link_libraries(DataStats PRIVATE skewer_matrix nlopt vendored_seqan vendored_rootpwa ${Boost_LIBRARIES} rt)
```

And for reseq:

```cmake
target_link_libraries(reseq PRIVATE DataStats DataStatsInterface GTest::gtest vendored_seqan ${Boost_LIBRARIES} rt)
```

- [ ] **Step 3: Remove link_directories for Boost**

In `CMakeLists.txt`, remove the Boost `link_directories()` calls (lines 40, 47-51):

```cmake
link_directories(${Boost_LIBRARY_DIRS})
...
if(NOT $ENV{BOOST_ROOT} STREQUAL "")
    link_directories("$ENV{BOOST_ROOT}/lib/")
endif()
if(NOT $ENV{BOOST_LIBRARYDIR} STREQUAL "")
    link_directories($ENV{BOOST_LIBRARYDIR})
endif()
```

Modern CMake's `find_package(Boost)` provides imported targets. If `BOOST_ROOT` is set, `find_package` finds it automatically.

- [ ] **Step 4: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass. No `CMAKE_CXX_FLAGS`, `include_directories()`, `link_directories()`, or `link_libraries()` remain in root `CMakeLists.txt`.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt reseq/CMakeLists.txt
git commit -m "build: replace global compiler flags with target-specific options"
```

---

### Task 7: Extract `reseq_lib` Static Library

**Files:**
- Modify: `reseq/CMakeLists.txt` (create reseq_lib from DataStats + DataStatsInterface + ProbabilityEstimates + Simulator)

- [ ] **Step 1: Merge DataStats, DataStatsInterface, and remaining source into reseq_lib**

In `reseq/CMakeLists.txt`, replace the existing `DataStats` and `DataStatsInterface` library definitions and the source files in the `reseq` executable with a single `reseq_lib`:

```cmake
set(PRIVATE_COMPILE_OPTIONS -Wall -Wextra -Wno-sign-compare -Wno-implicit-fallthrough -fPIC)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND PRIVATE_COMPILE_OPTIONS -fext-numeric-literals -Wno-deprecated-declarations)
endif()

# Core library — all production source except main.cpp and *Test.*
add_library(reseq_lib STATIC
  AdapterStats.cpp
  CoverageStats.cpp
  DataStats.cpp
  DataStatsInterface.cpp
  ErrorStats.cpp
  FragmentDistributionStats.cpp
  FragmentDuplicationStats.cpp
  ProbabilityEstimates.cpp
  QualityStats.cpp
  Reference.cpp
  Simulator.cpp
  Surrounding.cpp
  TileStats.cpp
)
target_include_directories(reseq_lib PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/reseq"
)
target_link_libraries(reseq_lib PUBLIC
  skewer_matrix nlopt vendored_seqan vendored_rootpwa ${Boost_LIBRARIES} rt
)
target_compile_options(reseq_lib PRIVATE ${PRIVATE_COMPILE_OPTIONS})

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Production executable — just main.cpp
add_executable(reseq main.cpp)
target_link_libraries(reseq PRIVATE reseq_lib)
target_compile_options(reseq PRIVATE ${PRIVATE_COMPILE_OPTIONS})

install(TARGETS reseq
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
```

Note: `ProbabilityEstimates.cpp` and `Simulator.cpp` were previously only in the `reseq` executable target — they now move into `reseq_lib`.

**Temporarily**, also keep the test sources and gtest link on the `reseq` target so `reseq test` still works:

Add after the reseq executable definition:

```cmake
# TEMPORARY: tests still in production binary (migrated in Task 8)
target_sources(reseq PRIVATE
  AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp
  ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp
  FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp
  QualityStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp
  SimulatorTest.cpp SurroundingTest.cpp TileStatsTest.cpp
  utilitiesTest.cpp VectTest.cpp
)
target_link_libraries(reseq PRIVATE GTest::gtest)
```

- [ ] **Step 2: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
build/bin/reseq test
```

Expected: All tests pass. The build now compiles `reseq_lib` as a separate static library.

- [ ] **Step 3: Commit**

```bash
git add reseq/CMakeLists.txt
git commit -m "build: extract reseq_lib static library from production sources"
```

---

### Task 8: Create Separate Test Binary with CTest

**Files:**
- Create: `reseq/test_main.cpp`
- Modify: `reseq/CMakeLists.txt` (add reseq_test target, enable CTest)
- Modify: `CMakeConfig.h.in` (add RESEQ_BINARY_DIR)

- [ ] **Step 1: Add RESEQ_BINARY_DIR to CMakeConfig.h.in**

In `CMakeConfig.h.in`, add a line for the binary dir (needed for regression tests in Phase 2 to find the `reseq` binary):

```c
#ifndef CMAKECONFIG_H
#define CMAKECONFIG_H

#define RESEQ_VERSION_MAJOR @RESEQ_VERSION_MAJOR@
#define RESEQ_VERSION_MINOR @RESEQ_VERSION_MINOR@
#define RESEQ_VERSION_PATCH @RESEQ_VERSION_PATCH@
#define RESEQ_VERSION "@RESEQ_VERSION@"
#define RESEQ_GIT_VERSION "@RESEQ_GIT_DESCRIBE@"
#define PROJECT_SOURCE_DIR "@PROJECT_SOURCE_DIR@"
#define RESEQ_BINARY_DIR "@CMAKE_RUNTIME_OUTPUT_DIRECTORY@"

#endif // CMAKECONFIG_H
```

- [ ] **Step 2: Create test_main.cpp**

Create `reseq/test_main.cpp`:

```cpp
#include <string>

#include "gtest/gtest.h"

#include "reportingUtils.hpp"

namespace reseq {
uint16_t kVerbosityLevel = 2;
bool kNoDebugOutput = false;
} // namespace reseq

#include "AdapterStatsTest.h"
#include "CoverageStatsTest.h"
#include "DataStatsTest.h"
#include "ErrorStatsTest.h"
#include "FragmentDistributionStatsTest.h"
#include "FragmentDuplicationStatsTest.h"
#include "ProbabilityEstimatesTest.h"
#include "QualityStatsTest.h"
#include "ReferenceTest.h"
#include "SeqQualityStatsTest.h"
#include "SimulatorTest.h"
#include "SurroundingTest.h"
#include "TileStatsTest.h"
#include "VectTest.h"
#include "utilitiesTest.h"

#include "BasicTestClass.hpp"
#include "utilities.hpp"

int main(int argc, char** argv) {
    std::string test_dir;
    if (!reseq::BasicTestClass::GetTestDir(test_dir)) {
        std::cerr << "ERROR: Cannot find test data directory. Aborting." << std::endl;
        return 1;
    }

    // Register all test classes
    reseq::AdapterStatsTest::Register();
    reseq::DataStatsTest::Register();
    reseq::FragmentDistributionStatsTest::Register(1); // single-threaded for determinism
    reseq::FragmentDuplicationStatsTest::Register();
    reseq::ProbabilityEstimatesTest::Register();
    reseq::ReferenceTest::Register();
    reseq::SeqQualityStatsTest::Register();
    reseq::SimulatorTest::Register();
    reseq::TileStatsTest::Register();
    reseq::VectTest::Register();
    reseq::utilitiesTest::Register();

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
```

- [ ] **Step 3: Add reseq_test target to reseq/CMakeLists.txt**

Add after the `reseq` executable block:

```cmake
# Test executable — links reseq_lib + GoogleTest
option(RESEQ_BUILD_TESTS "Build test binary" ON)
if(RESEQ_BUILD_TESTS)
  enable_testing()
  add_executable(reseq_test
    test_main.cpp
    AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp
    ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp
    FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp
    QualityStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp
    SimulatorTest.cpp SurroundingTest.cpp TileStatsTest.cpp
    utilitiesTest.cpp VectTest.cpp
  )
  target_link_libraries(reseq_test PRIVATE reseq_lib GTest::gtest)
  target_compile_options(reseq_test PRIVATE ${PRIVATE_COMPILE_OPTIONS})
  include(GoogleTest)
  gtest_discover_tests(reseq_test)
endif()
```

- [ ] **Step 4: Build and run both test paths**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)

# Run via new test binary
cd build && ctest --output-on-failure && cd ..

# Run via old path (should still work)
build/bin/reseq test
```

Expected: Both paths pass all tests. CTest discovers and runs the same tests.

- [ ] **Step 5: Verify test count matches baseline**

```bash
cd build && ctest --output-on-failure 2>&1 | grep -c "Test #" && cd ..
```

Compare the count against the baseline from Task 1.

- [ ] **Step 6: Commit**

```bash
git add reseq/test_main.cpp reseq/CMakeLists.txt CMakeConfig.h.in
git commit -m "build: create separate reseq_test binary with CTest integration"
```

---

### Task 9: Remove Test Code from Production Binary

**Files:**
- Modify: `reseq/CMakeLists.txt` (remove test sources and gtest from reseq target)
- Modify: `reseq/main.cpp` (remove test subcommand, gtest include, test registration)

- [ ] **Step 1: Remove test sources from reseq target**

In `reseq/CMakeLists.txt`, remove the `# TEMPORARY` block added in Task 7:

```cmake
# TEMPORARY: tests still in production binary (migrated in Task 8)
target_sources(reseq PRIVATE
  AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp
  ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp
  FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp
  QualityStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp
  SimulatorTest.cpp SurroundingTest.cpp TileStatsTest.cpp
  utilitiesTest.cpp VectTest.cpp
)
target_link_libraries(reseq PRIVATE GTest::gtest)
```

Delete those lines entirely.

- [ ] **Step 2: Remove gtest include from main.cpp**

In `reseq/main.cpp`, remove line 19:

```cpp
#include "gtest/gtest.h"
```

- [ ] **Step 3: Remove test header includes from main.cpp**

Search for and remove all test header includes from `main.cpp`. Look for includes like:

```cpp
#include "AdapterStatsTest.h"
#include "CoverageStatsTest.h"
```

etc. (Check if there are any — the current `main.cpp` includes them via the test `.cpp` files, not directly. If none found, skip this step.)

- [ ] **Step 4: Remove the test subcommand from main.cpp**

In `reseq/main.cpp`, remove the entire `} else if ("test" == unrecognized_opts.at(0)) {` block (lines 1179-1239). This includes:
- The test mode entry check
- Test registration calls (Register)
- GoogleTest initialization
- RunGoogleTests calls
- The "tests have failed" error message

After removal, the `else` chain should go directly from the last valid command to the "Unrecognized command" error:

```cpp
                }
            }
        } else {
            if (2 < kVerbosityLevel) {
                cerr << std::endl;
            }
            printErr << "Unrecognized command: '" << unrecognized_opts.at(0) << "'" << std::endl;
```

- [ ] **Step 5: Remove RunGoogleTests helper function**

Search for the `RunGoogleTests` function definition in `main.cpp` and remove it. It's a helper used only by the test subcommand.

- [ ] **Step 6: Remove test-related includes that are no longer needed**

After removing the test block, some includes in `main.cpp` may become unused. Check if these are still needed:
- `#include "gtest/gtest.h"` (already removed in Step 2)
- Any `*Test.h` includes

The `#include "BasicTestClass.hpp"` include was likely not in main.cpp directly (tests pull it in), but verify and remove if present.

- [ ] **Step 7: Build and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)

# New binary should NOT have test subcommand
build/bin/reseq test 2>&1 | head -5
# Expected: "Unrecognized command: 'test'"

# Tests run via ctest only
cd build && ctest --output-on-failure && cd ..
```

Expected: `reseq test` now prints "Unrecognized command". All tests pass via `ctest`.

- [ ] **Step 8: Verify production binary size decreased**

```bash
ls -la build/bin/reseq
```

The binary should be noticeably smaller without GoogleTest and 15 test files linked in.

- [ ] **Step 9: Commit**

```bash
git add reseq/CMakeLists.txt reseq/main.cpp
git commit -m "build: remove test code from production binary

Tests now run exclusively via the reseq_test binary and ctest.
The reseq binary no longer links GoogleTest or compiles test sources."
```

---

### Task 10: Add CODE_COVERAGE CMake Option

**Files:**
- Modify: `reseq/CMakeLists.txt` (add coverage flags)
- Modify: `CMakeLists.txt` (add option)

- [ ] **Step 1: Add the option to root CMakeLists.txt**

In `CMakeLists.txt`, add near the existing `RESEQ_BUILD_PYTHON` option:

```cmake
option(RESEQ_BUILD_PYTHON "Build Python bindings (requires SWIG 3+ and python3-dev)" OFF)
option(CODE_COVERAGE "Enable code coverage instrumentation (GCC only)" OFF)
```

- [ ] **Step 2: Add coverage flags to reseq_lib in reseq/CMakeLists.txt**

In `reseq/CMakeLists.txt`, after the `reseq_lib` target definition, add:

```cmake
if(CODE_COVERAGE)
  target_compile_options(reseq_lib PRIVATE --coverage)
  target_link_options(reseq_lib PRIVATE --coverage)
  target_compile_options(reseq_test PRIVATE --coverage)
  target_link_options(reseq_test PRIVATE --coverage)
endif()
```

- [ ] **Step 3: Build with coverage disabled (default) and test**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure && cd ..
```

Expected: All tests pass, no coverage files generated.

- [ ] **Step 4: Build with coverage enabled and verify instrumentation**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCODE_COVERAGE=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure && cd ..
find build -name "*.gcda" | head -5
```

Expected: `.gcda` files exist after running tests, proving coverage instrumentation works.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt reseq/CMakeLists.txt
git commit -m "build: add CODE_COVERAGE CMake option for gcov instrumentation"
```

---

### Task 11: Update Makefile

**Files:**
- Modify: `Makefile`

- [ ] **Step 1: Update the test target to use ctest**

In `Makefile`, change the `test` target from:

```makefile
test: build
	$(BUILD_DIR)/bin/reseq test
```

to:

```makefile
test: build
	cd $(BUILD_DIR) && ctest --output-on-failure -j$$(nproc)
```

- [ ] **Step 2: Add a coverage target**

Add a new target:

```makefile
coverage:
	cmake -S . -B $(BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCODE_COVERAGE=ON \
		$(CMAKE_FLAGS)
	cmake --build $(BUILD_DIR) -j$$(nproc)
	cd $(BUILD_DIR) && ctest --output-on-failure
	lcov --capture --directory $(BUILD_DIR) --output-file $(BUILD_DIR)/coverage.info \
		--ignore-errors mismatch
	lcov --remove $(BUILD_DIR)/coverage.info \
		'*/googletest/*' '*/seqan/*' '*/nlopt/*' '*/skewer/*' '*/2016-05-15_ROOTPWA/*' \
		'/usr/*' '*/build/*' \
		--output-file $(BUILD_DIR)/coverage.info --ignore-errors unused
	genhtml $(BUILD_DIR)/coverage.info --output-directory $(BUILD_DIR)/coverage-report
	@echo "Coverage report: $(BUILD_DIR)/coverage-report/index.html"
```

- [ ] **Step 3: Update help text**

Add the coverage target to the help output:

```makefile
help:
	@echo "Targets:"
	@echo "  build        Configure and build (default)"
	@echo "  test         Build and run unit tests via ctest"
	@echo "  format       Format C++ and Python files in-place"
	@echo "  format-check Dry-run format check (CI use)"
	@echo "  lint         Run clang-tidy and ruff"
	@echo "  coverage     Build with coverage, run tests, generate HTML report"
	@echo "  changelog    Generate CHANGELOG.md from git history"
	@echo "  clean        Remove build directory"
	@echo "  install      Build and install"
```

- [ ] **Step 4: Test the updated Makefile**

```bash
make clean
make test
```

Expected: Builds and runs tests via ctest.

- [ ] **Step 5: Commit**

```bash
git add Makefile
git commit -m "build: update Makefile to use ctest and add coverage target"
```

---

### Task 12: Update CI Pipeline

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Update the Test step to use ctest**

In `.github/workflows/ci.yml`, change the Test step (line 44-45) from:

```yaml
      - name: Test
        run: build/bin/reseq test
```

to:

```yaml
      - name: Test
        run: cd build && ctest --output-on-failure -j$(nproc)
```

- [ ] **Step 2: Add a sanitizer job**

Add a new job after `build-and-test`:

```yaml
  sanitizers:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y gcc-13 g++-13 \
            libboost-all-dev libbz2-dev zlib1g-dev

      - name: Configure with sanitizers
        env:
          CC: gcc-13
          CXX: g++-13
        run: |
          cmake -S . -B build \
            -DCMAKE_BUILD_TYPE=Debug \
            -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"

      - name: Build
        run: cmake --build build -j$(nproc)

      - name: Test with sanitizers
        run: cd build && ctest --output-on-failure -j$(nproc)
```

- [ ] **Step 3: Add a coverage job**

Add another job:

```yaml
  coverage:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y gcc-13 g++-13 lcov \
            libboost-all-dev libbz2-dev zlib1g-dev

      - name: Configure with coverage
        env:
          CC: gcc-13
          CXX: g++-13
        run: cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCODE_COVERAGE=ON

      - name: Build
        run: cmake --build build -j$(nproc)

      - name: Test
        run: cd build && ctest --output-on-failure

      - name: Generate coverage report
        run: |
          lcov --capture --directory build --output-file build/coverage.info \
            --ignore-errors mismatch
          lcov --remove build/coverage.info \
            '*/googletest/*' '*/seqan/*' '*/nlopt/*' '*/skewer/*' \
            '*/2016-05-15_ROOTPWA/*' '/usr/*' '*/build/*' \
            --output-file build/coverage.info --ignore-errors unused
          lcov --list build/coverage.info

      - name: Upload coverage artifact
        uses: actions/upload-artifact@v4
        with:
          name: coverage-report
          path: build/coverage.info
```

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: switch to ctest, add sanitizer and coverage jobs"
```

---

### Task 13: Final Verification

**Files:**
- None modified — verification only

- [ ] **Step 1: Clean build from scratch**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
```

Expected: Builds cleanly with no warnings from reseq sources.

- [ ] **Step 2: Run all tests via ctest**

```bash
cd build && ctest --output-on-failure -j$(nproc) && cd ..
```

Expected: All tests pass.

- [ ] **Step 3: Verify production binary has no test code**

```bash
# Binary should not contain GoogleTest symbols
nm build/bin/reseq | grep -i gtest | head -5
```

Expected: No output (no GoogleTest symbols).

```bash
# Test subcommand should be gone
build/bin/reseq test 2>&1 | head -1
```

Expected: "Unrecognized command: 'test'" or similar error.

- [ ] **Step 4: Verify test binary exists and works**

```bash
build/bin/reseq_test --gtest_list_tests | head -20
```

Expected: Lists all test cases.

- [ ] **Step 5: Verify ctest discovers tests**

```bash
cd build && ctest -N && cd ..
```

Expected: Lists all discovered tests with `Test #N:` prefix.

- [ ] **Step 6: Verify coverage instrumentation**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCODE_COVERAGE=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure && cd ..
find build -name "*.gcda" | wc -l
```

Expected: Non-zero count of `.gcda` files.

- [ ] **Step 7: Verify format check still passes**

```bash
make format-check
```

Expected: No formatting issues.

- [ ] **Step 8: Compare test count to baseline**

```bash
cd build && ctest -N 2>&1 | tail -1 && cd ..
cat test-baseline-count.txt
```

Expected: Same number of tests as the baseline captured in Task 1.

- [ ] **Step 9: Clean up baseline files**

```bash
rm -f test-baseline.log test-baseline-count.txt
```

---

## Summary

After completing all 13 tasks, the build system state will be:

| Aspect | Before | After |
|--------|--------|-------|
| CMake minimum | 3.1 | 3.16 |
| C++ standard | C++14 | C++20 |
| GoogleTest | Vendored subdirectory | FetchContent v1.15.2 |
| NLopt | Vendored subdirectory | FetchContent v2.9.1 |
| SeqAn2 | Global includes | `vendored_seqan` INTERFACE target with SYSTEM |
| Skewer | Unchanged | Unchanged (already a proper target) |
| ROOTPWA | Global includes | `vendored_rootpwa` INTERFACE target (temporary) |
| Library targets | DataStats, DataStatsInterface | `reseq_lib` (single static library) |
| Test location | Inside `reseq` binary | Separate `reseq_test` binary |
| Test runner | `reseq test` subcommand | `ctest` |
| Global CMake state | 8+ global commands | Zero |
| CI sanitizers | None | ASan + UBSan job |
| CI coverage | None | gcov/lcov job |
| `CODE_COVERAGE` option | N/A | Available (OFF by default) |
