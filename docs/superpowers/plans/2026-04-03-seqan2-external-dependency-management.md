# SeqAn2 External Dependency Management Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade ReSeq from vendored SeqAn `2.4.0` to external SeqAn `2.5.2`, stop carrying SeqAn source in-tree, and establish a reusable dependency policy for the repository's other third-party libraries.

**Architecture:** Keep ReSeq's build centered on CMake-native `find_package()` and imported targets. Production dependencies should be consumed through provider-neutral package configs so the same CMake code works with system packages, conda environments, Conan-generated configs, and vcpkg toolchains. `FetchContent` remains allowed for development-only or temporary-gap dependencies; vendoring is reserved for patched forks like `skewer`.

**Tech Stack:** CMake 3.16+, C++20, SeqAn 2.5.2, imported CMake targets, CMake Presets, conda/bioconda, Conan-compatible package configs, vcpkg-compatible toolchain workflows

**Research Basis:**
- CMake `find_package()`: https://cmake.org/cmake/help/latest/command/find_package.html
- CMake `FetchContent`: https://cmake.org/cmake/help/latest/module/FetchContent.html
- SeqAn `2.5.2` release: https://github.com/seqan/seqan/releases
- SeqAn imported-target guidance (`seqan::seqan2`): https://github.com/seqan/seqan/releases/tag/seqan-v2.5.1
- Bioconda `seqan` package (`2.5.2`): https://anaconda.org/bioconda/seqan
- Conan/CMakeDeps model: https://docs.conan.io/2/reference/tools/cmake/cmakedeps.html
- vcpkg/CMake integration: https://learn.microsoft.com/en-gb/vcpkg/users/buildsystems/cmake-integration

---

## File Map

| Action | File | Purpose |
|--------|------|---------|
| Modify | `CMakeLists.txt` | Remove hard-wired vendored SeqAn setup and delegate dependency resolution |
| Create | `cmake/ReSeqDependencies.cmake` | Centralize dependency policy and imported-target wrappers |
| Create | `CMakePresets.json` | Add provider-friendly configure presets |
| Modify | `reseq/CMakeLists.txt` | Link ReSeq via project wrapper targets instead of raw SeqAn variables |
| Modify | `python/CMakeLists.txt` | Link Python module through project wrapper target instead of `${SEQAN_LIBRARIES}` |
| Modify | `README.md` | Replace vendored SeqAn install assumptions with external dependency instructions |
| Modify | `CLAUDE.md` | Update repository guidance and dependency policy |
| Modify | `conda/reseq/meta.yaml` | Consume external SeqAn package in the existing conda workflow |
| Delete | `seqan/` | Remove vendored SeqAn after external path is proven in CI |

---

## Dependency Policy

Use this policy for SeqAn and future non-patched dependencies:

1. Prefer `find_package()` with imported targets.
2. Keep project code provider-neutral.
3. Allow `FetchContent` only when no acceptable package path exists or for dev-only/test-only libraries.
4. Vendor only patched forks or sources that the project intentionally owns.

This matches the current repository direction:
- already uses `find_package()` for `Boost`, `ZLIB`, `BZip2`, `SWIG`
- already uses `FetchContent` for `googletest` and `NLopt`
- already vendors only the difficult cases (`seqan`, patched `skewer`)

---

## Dependency Classification

Apply the new dependency model repo-wide with these buckets:

### External-required

These should be consumed via `find_package()` and imported targets, with no vendored source tree in the repository:

- `SeqAn`
- `Boost`
- `ZLIB`
- `BZip2`
- `Python`
- `SWIG`

### External-with-fallback

These may use `find_package()` first and fall back to `FetchContent` when no suitable package is available:

- `NLopt`

### Vendored-by-exception

These remain vendored until their local patch burden is removed or moved to a maintained fork:

- `skewer`

### Development-only FetchContent

These are acceptable as `FetchContent` dependencies because they are test/developer tooling, not production runtime dependencies:

- `GoogleTest`

---

### Task 1: Capture the Current SeqAn Baseline

**Files:**
- None modified — verification only

- [ ] **Step 1: Record current vendored SeqAn version**

Run:

```bash
sed -n '40,46p' seqan/include/seqan/version.h
```

Expected: `SEQAN_VERSION_MAJOR 2`, `SEQAN_VERSION_MINOR 4`, `SEQAN_VERSION_PATCH 0`.

- [ ] **Step 2: Record the current build/test baseline**

Run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: configure succeeds, build succeeds, `reseq_tests` passes.

- [ ] **Step 3: Record the external target expectation from upstream**

Open the SeqAn `2.5.1/2.5.2` release notes and capture the guidance that consumers should link `seqan::seqan2` instead of manually applying `SEQAN_*` variables.

Expected: this note is copied into the implementation PR description or issue for rationale.

- [ ] **Step 4: Commit**

No code changes in this task. Use it as a checkpoint only.

---

### Task 2: Centralize Dependency Resolution in `cmake/ReSeqDependencies.cmake`

**Files:**
- Create: `cmake/ReSeqDependencies.cmake`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Create the dependency policy module**

Create `cmake/ReSeqDependencies.cmake` with this initial structure:

```cmake
include(FetchContent)

option(RESEQ_ALLOW_FETCHCONTENT "Allow downloading missing dependencies at configure time" ON)
option(RESEQ_USE_VENDORED_SEQAN "Use the in-tree SeqAn fallback" OFF)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Boost 1.48.0 REQUIRED
  COMPONENTS filesystem iostreams math_c99 math_c99f math_c99l math_tr1 math_tr1f math_tr1l
             program_options serialization system)

add_library(reseq_zlib INTERFACE)
target_link_libraries(reseq_zlib INTERFACE ZLIB::ZLIB)

add_library(reseq_bzip2 INTERFACE)
target_link_libraries(reseq_bzip2 INTERFACE BZip2::BZip2)

add_library(reseq_boost INTERFACE)
target_link_libraries(reseq_boost INTERFACE ${Boost_LIBRARIES})
```

- [ ] **Step 2: Move GoogleTest and NLopt setup behind the module**

Extend `cmake/ReSeqDependencies.cmake` with the existing `googletest` and `NLopt` resolution logic from the root `CMakeLists.txt`.

Use this shape for `googletest`:

```cmake
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
```

Use this shape for `NLopt`:

```cmake
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

- [ ] **Step 3: Load the dependency module from the root project**

In `CMakeLists.txt`, remove the direct dependency setup block and replace it with:

```cmake
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(ReSeqDependencies)
```

- [ ] **Step 4: Reconfigure and verify no behavior changed yet**

Run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: same pass result as Task 1.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cmake/ReSeqDependencies.cmake
git commit -m "build: centralize dependency resolution"
```

---

### Task 2A: Modernize Existing External Dependencies to Match the Policy

**Files:**
- Modify: `cmake/ReSeqDependencies.cmake`
- Modify: `reseq/CMakeLists.txt`
- Modify: `python/CMakeLists.txt`

- [ ] **Step 1: Wrap Boost in a project-level imported-interface target**

In `cmake/ReSeqDependencies.cmake`, prefer a project wrapper target over leaking raw `${Boost_LIBRARIES}` through the tree:

```cmake
add_library(reseq_boost INTERFACE)
target_link_libraries(reseq_boost INTERFACE ${Boost_LIBRARIES})
```

Then in `reseq/CMakeLists.txt`, replace direct `${Boost_LIBRARIES}` usage with `reseq_boost`.

- [ ] **Step 2: Keep ZLIB and BZip2 consumed through imported targets**

Ensure production targets link through the wrapper targets or directly through:

```cmake
ZLIB::ZLIB
BZip2::BZip2
```

Expected: no manual include directories or link directories are added for these libraries.

- [ ] **Step 3: Keep NLopt as a hybrid dependency**

Do not vendor `NLopt`. Preserve the pattern:

```cmake
find_package(NLopt QUIET)
```

with `FetchContent` fallback when it is not present.

Expected: this dependency remains provider-neutral and consistent with the new policy.

- [ ] **Step 4: Leave GoogleTest as FetchContent**

Do not convert `GoogleTest` to a required system dependency.

Expected: test-only dependencies remain easy to bootstrap and do not complicate production packaging.

- [ ] **Step 5: Leave `skewer` vendored**

Do not de-vendor `skewer` in this migration. Add a short comment in `CLAUDE.md` and `README.md` that its vendored status is intentional because the project carries local modifications.

- [ ] **Step 6: Verify the build still works after dependency-target cleanup**

Run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: no behavior change beyond cleaner dependency boundaries.

- [ ] **Step 7: Commit**

```bash
git add cmake/ReSeqDependencies.cmake reseq/CMakeLists.txt python/CMakeLists.txt README.md CLAUDE.md
git commit -m "build: classify and normalize external dependency usage"
```

---

### Task 3: Replace Vendored SeqAn Wiring with an External-First Imported Target

**Files:**
- Modify: `cmake/ReSeqDependencies.cmake`
- Modify: `CMakeLists.txt`
- Modify: `reseq/CMakeLists.txt`
- Modify: `python/CMakeLists.txt`

- [ ] **Step 1: Add external-first SeqAn resolution**

In `cmake/ReSeqDependencies.cmake`, add this logic after the common dependency block:

```cmake
set(_reseq_seqan_target "")

find_package(SeqAn 2.5 CONFIG QUIET)
if(TARGET seqan::seqan2)
  set(_reseq_seqan_target seqan::seqan2)
endif()

if(NOT _reseq_seqan_target)
  find_package(SeqAn 2.5 QUIET)
  if(TARGET seqan::seqan2)
    set(_reseq_seqan_target seqan::seqan2)
  endif()
endif()

if(NOT _reseq_seqan_target AND RESEQ_ALLOW_FETCHCONTENT)
  FetchContent_Declare(
    seqan2
    GIT_REPOSITORY https://github.com/seqan/seqan.git
    GIT_TAG seqan-v2.5.2
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(seqan2)
  if(EXISTS "${seqan2_SOURCE_DIR}/include")
    add_library(reseq_seqan INTERFACE)
    target_include_directories(reseq_seqan SYSTEM INTERFACE "${seqan2_SOURCE_DIR}/include")
    set(_reseq_seqan_target reseq_seqan)
  endif()
endif()

if(NOT _reseq_seqan_target AND RESEQ_USE_VENDORED_SEQAN)
  add_library(reseq_seqan INTERFACE)
  target_include_directories(reseq_seqan SYSTEM INTERFACE "${PROJECT_SOURCE_DIR}/seqan/include")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake")
  find_package(SeqAn REQUIRED)
  target_include_directories(reseq_seqan SYSTEM INTERFACE ${SEQAN_INCLUDE_DIRS})
  target_compile_definitions(reseq_seqan INTERFACE ${SEQAN_DEFINITIONS})
  separate_arguments(SEQAN_CXX_FLAGS_LIST NATIVE_COMMAND "${SEQAN_CXX_FLAGS}")
  target_compile_options(reseq_seqan INTERFACE ${SEQAN_CXX_FLAGS_LIST})
  target_link_libraries(reseq_seqan INTERFACE ${SEQAN_LIBRARIES})
  set(_reseq_seqan_target reseq_seqan)
endif()

if(NOT _reseq_seqan_target)
  message(FATAL_ERROR
    "SeqAn >= 2.5 was not found. Provide it via package manager / CMAKE_PREFIX_PATH, "
    "enable RESEQ_ALLOW_FETCHCONTENT, or explicitly set RESEQ_USE_VENDORED_SEQAN=ON.")
endif()

if(NOT TARGET reseq_seqan)
  add_library(reseq_seqan INTERFACE)
  target_link_libraries(reseq_seqan INTERFACE ${_reseq_seqan_target})
endif()
```

- [ ] **Step 2: Remove the root-level vendored SeqAn block**

Delete this block from `CMakeLists.txt`:

```cmake
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/seqan/util/cmake/")
find_package(SeqAn REQUIRED)
add_library(vendored_seqan INTERFACE)
target_include_directories(vendored_seqan SYSTEM INTERFACE
  "${PROJECT_SOURCE_DIR}/seqan/include"
  ${SEQAN_INCLUDE_DIRS}
)
target_compile_definitions(vendored_seqan INTERFACE ${SEQAN_DEFINITIONS})
separate_arguments(SEQAN_CXX_FLAGS_LIST NATIVE_COMMAND "${SEQAN_CXX_FLAGS}")
target_compile_options(vendored_seqan INTERFACE ${SEQAN_CXX_FLAGS_LIST})
target_link_libraries(vendored_seqan INTERFACE ${SEQAN_LIBRARIES})
```

- [ ] **Step 3: Relink production code through the project target**

In `reseq/CMakeLists.txt`, replace `vendored_seqan` with `reseq_seqan`:

```cmake
target_link_libraries(reseq_lib PUBLIC
  skewer_matrix nlopt reseq_seqan GTest::gtest
  ${Boost_LIBRARIES} ZLIB::ZLIB BZip2::BZip2 rt
)
```

- [ ] **Step 4: Relink the Python module through the project target**

In `python/CMakeLists.txt`, replace:

```cmake
SWIG_LINK_LIBRARIES(DataStats ${SEQAN_LIBRARIES})
```

with:

```cmake
SWIG_LINK_LIBRARIES(DataStats reseq_seqan)
```

- [ ] **Step 5: Verify the vendored fallback still works**

Run:

```bash
cmake -S . -B build -DRESEQ_USE_VENDORED_SEQAN=ON -DRESEQ_ALLOW_FETCHCONTENT=OFF
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: build and tests still pass through the compatibility path.

- [ ] **Step 6: Verify the external path works**

Run in a conda environment that contains `seqan=2.5.2`:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DRESEQ_USE_VENDORED_SEQAN=OFF \
  -DRESEQ_ALLOW_FETCHCONTENT=OFF \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: configure succeeds without reading `seqan/` from the repository.

- [ ] **Step 7: Commit**

```bash
git add CMakeLists.txt cmake/ReSeqDependencies.cmake reseq/CMakeLists.txt python/CMakeLists.txt
git commit -m "build: consume SeqAn 2.5 externally via imported target"
```

---

### Task 4: Add Provider-Friendly Configure Presets

**Files:**
- Create: `CMakePresets.json`

- [ ] **Step 1: Create repository presets for system and conda workflows**

Create `CMakePresets.json`:

```json
{
  "version": 6,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23,
    "patch": 0
  },
  "configurePresets": [
    {
      "name": "base",
      "hidden": true,
      "generator": "Unix Makefiles",
      "binaryDir": "${sourceDir}/build",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "RelWithDebInfo",
        "RESEQ_USE_VENDORED_SEQAN": "OFF",
        "RESEQ_ALLOW_FETCHCONTENT": "ON"
      }
    },
    {
      "name": "system-release",
      "inherits": "base"
    },
    {
      "name": "conda-release",
      "inherits": "base",
      "cacheVariables": {
        "CMAKE_PREFIX_PATH": "$env{CONDA_PREFIX}"
      }
    }
  ],
  "buildPresets": [
    {
      "name": "build",
      "configurePreset": "system-release"
    }
  ],
  "testPresets": [
    {
      "name": "test",
      "configurePreset": "system-release",
      "output": {
        "outputOnFailure": true
      }
    }
  ]
}
```

- [ ] **Step 2: Document provider compatibility**

Do not hard-code Conan or vcpkg into project CMake yet. Instead, verify in docs that the same `find_package()`-based project build works when users provide:
- `CMAKE_PREFIX_PATH` from conda or system package roots
- Conan-generated package configs
- vcpkg toolchain files

Expected: no provider-specific code appears in `CMakeLists.txt`.

- [ ] **Step 3: Verify presets**

Run:

```bash
cmake --preset system-release
cmake --build build -j"$(nproc)"
ctest --preset test
```

If testing the conda path:

```bash
cmake --preset conda-release
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: preset-driven configure/build/test succeeds.

- [ ] **Step 4: Commit**

```bash
git add CMakePresets.json
git commit -m "build: add provider-friendly CMake presets"
```

---

### Task 5: Update Docs and Packaging Around the External Dependency Model

**Files:**
- Modify: `README.md`
- Modify: `CLAUDE.md`
- Modify: `conda/reseq/meta.yaml`

- [ ] **Step 1: Update README dependency policy**

In `README.md`, replace the implicit vendored-SeqAn assumption with text like:

```md
ReSeq consumes SeqAn externally and expects a CMake-discoverable SeqAn 2.5.x installation.
Supported provider workflows:

- system packages exposing `find_package(SeqAn ...)`
- bioconda environments via `CMAKE_PREFIX_PATH=$CONDA_PREFIX`
- Conan-generated CMake package configs
- vcpkg toolchain integration

Vendored source trees are reserved for patched forks such as `skewer`.
```

- [ ] **Step 2: Add concrete install examples**

Add these examples to `README.md`:

```bash
# system / manually installed packages
cmake --preset system-release
cmake --build build -j$(nproc)

# conda / bioconda
micromamba create -n reseq-build -c conda-forge -c bioconda \
  cmake make compilers boost zlib bzip2 seqan
micromamba activate reseq-build
cmake --preset conda-release
cmake --build build -j$(nproc)
```

- [ ] **Step 3: Update repository guidance in `CLAUDE.md`**

Replace the dependency summary so it reflects:
- SeqAn is no longer listed as vendored
- `skewer` remains vendored because it carries local modifications
- `find_package()` plus imported targets is the default pattern for production dependencies

- [ ] **Step 4: Update the conda recipe to consume external SeqAn**

In `conda/reseq/meta.yaml`, add `seqan >=2.5.2,<2.6` to the build and host requirements:

```yaml
requirements:
  build:
    - cmake
    - make
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - boost
    - swig
    - seqan >=2.5.2,<2.6
  host:
    - python
    - zlib
    - bzip2
    - boost
    - seqan >=2.5.2,<2.6
```

- [ ] **Step 5: Validate docs and recipe**

Run:

```bash
rg -n "vendored SeqAn|seqan/" README.md CLAUDE.md conda/reseq/meta.yaml
```

Expected: references to SeqAn vendoring are gone except in migration notes.

- [ ] **Step 6: Commit**

```bash
git add README.md CLAUDE.md conda/reseq/meta.yaml
git commit -m "docs: document external SeqAn dependency workflow"
```

---

### Task 6: Remove the Vendored SeqAn Tree After External CI Is Green

**Files:**
- Delete: `seqan/`
- Modify: `README.md`
- Modify: `CLAUDE.md`
- Modify: any formatting/lint exclusions that still mention `seqan/`

- [ ] **Step 1: Add one CI path that forbids vendored SeqAn**

Whichever CI workflow owns the primary build should configure with:

```bash
cmake -S . -B build \
  -DRESEQ_USE_VENDORED_SEQAN=OFF \
  -DRESEQ_ALLOW_FETCHCONTENT=OFF \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: the project proves it can build without `seqan/` in-tree.

- [ ] **Step 2: Delete the vendored tree**

Run:

```bash
rm -rf seqan
```

Expected: `git status` shows the full directory removed.

- [ ] **Step 3: Remove any stale build references**

Run:

```bash
rg -n "seqan/" . \
  --glob '!docs/superpowers/plans/**' \
  --glob '!docs/superpowers/specs/**'
```

Expected: only intentional historical references remain.

- [ ] **Step 4: Rebuild and retest from a clean tree**

Run:

```bash
rm -rf build
cmake --preset conda-release
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

Expected: clean external build passes with no vendored SeqAn directory present.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "build: remove vendored SeqAn tree"
```

---

## Acceptance Criteria

- ReSeq no longer requires `seqan/` to be present in the repository.
- The default CMake code path uses `find_package()` and imported targets.
- SeqAn is upgraded to `2.5.2`.
- The build still passes with current `g++`/`cmake` toolchain and existing tests.
- The same project CMake remains compatible with system packages, bioconda, Conan-generated CMake configs, and vcpkg toolchain workflows.
- The project keeps vendoring only the dependencies it truly owns or patches.

---

## Rollout Notes

- Land the central dependency module first.
- Normalize existing external dependencies (`Boost`, `ZLIB`, `BZip2`, `NLopt`, `GoogleTest`) under the same policy before removing vendored SeqAn.
- Keep `RESEQ_USE_VENDORED_SEQAN=ON` available during one transition cycle.
- Remove the vendored tree only after one external-only CI path is stable.
- Do not apply the same de-vendoring policy to `skewer` until its local patch set is audited and either upstreamed or split into a maintained fork.
