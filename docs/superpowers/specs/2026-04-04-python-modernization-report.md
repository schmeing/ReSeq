# Python Modernization Report

**Date:** 2026-04-04

## Scope

This report covers the Python-facing parts of ReSeq:

- the optional SWIG bindings in [python/DataStats.i](/home/bernt-popp/development/ReSeq/python/DataStats.i)
- the Python build integration in [python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt)
- the helper scripts [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py) and [python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py)

The goal is to answer three questions:

1. Is Python still needed here?
2. What is outdated or fragile?
3. What should be modernized, and in what order?

## Executive Summary

The Python layer is optional and not required for the core `reseq` C++ binary. It exists for two purposes:

- a SWIG-generated `DataStats` Python module
- helper scripts for plotting and FASTQ name normalization

There is real modernization work worth doing, but it should be treated as a focused cleanup project, not as part of the SeqAn migration itself.

The most important finding is that [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py) still uses `time.clock()`, which was removed in Python 3.8. That makes the script incompatible with modern Python runtimes unless it is already patched locally or simply unused.

The build integration is also dated: [python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt) uses legacy `PythonLibs` and old-style SWIG commands, and still links raw `${Boost_LIBRARIES}` / `${SEQAN_LIBRARIES}` instead of modern targets.

## What Exists Today

### Python Build Surface

The Python subdirectory is only enabled when `RESEQ_BUILD_PYTHON=ON` in [CMakeLists.txt](/home/bernt-popp/development/ReSeq/CMakeLists.txt). The current implementation:

- locates an interpreter via `find_program(PYTHON ...)`
- requires `SWIG 3`
- uses `FIND_PACKAGE(PythonLibs 3 REQUIRED)`
- builds a SWIG module from [python/DataStats.i](/home/bernt-popp/development/ReSeq/python/DataStats.i)
- installs helper scripts to `bin`

This means Python is a build-time optional feature, not a mandatory project dependency.

### Python Code Surface

There are two very different Python maintenance targets:

- **Bindings:** [python/DataStats.i](/home/bernt-popp/development/ReSeq/python/DataStats.i)
- **Scripts:** [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py), [python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py)

The scripts are likely easier to keep alive than the SWIG bindings. The SWIG module is tightly coupled to the C++ data model and build system.

## Findings

### High Priority

#### 1. `plotDataStats.py` is incompatible with modern Python

[python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py) imports:

```python
from time import clock
```

and still calls `clock()` at multiple sites, including:

- [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py#L1005)
- [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py#L1019)
- [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py#L1236)

According to Python’s official documentation, `time.clock()` was deprecated in Python 3.3 and removed in Python 3.8. The documented replacements are `time.perf_counter()` or `time.process_time()`, depending on intent.

Impact:

- the plotting script is broken on current Python versions
- this is the clearest concrete modernization requirement in the Python layer

#### 2. The CMake Python integration is legacy

[python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt) still uses:

- `FIND_PACKAGE(PythonLibs 3 REQUIRED)` at line 14
- `SWIG_ADD_MODULE(...)` at line 31
- repeated `SWIG_LINK_LIBRARIES(...)` calls at lines 33-36

This is outdated compared with modern CMake practice:

- `find_package(Python3 COMPONENTS Interpreter Development REQUIRED)` is the current standard
- `swig_add_library(...)` is the modern `UseSWIG` interface
- target-based linking should be used instead of raw library variables

Impact:

- more fragile configuration behavior
- harder future maintenance
- unnecessary inconsistency with the rest of the repository’s ongoing CMake modernization

### Medium Priority

#### 3. The Python scripts are old-style but serviceable

[python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py) is not obviously broken, but it is dated:

- uses `getopt` instead of `argparse`
- detects gzip only by filename suffix
- has minimal validation and no tests

This script is small and easy to modernize. The main risk is not complexity, but lack of automated verification.

#### 4. The SWIG binding is tightly coupled and likely untested

[python/DataStats.i](/home/bernt-popp/development/ReSeq/python/DataStats.i) exposes C++ types directly, includes project headers in the wrapper body, and depends on generated code plus a custom “swigtrick” target from [python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt#L45).

Nothing in the current repository structure suggests this binding is exercised in CI.

Impact:

- breakage risk is high if C++ headers change
- build rot can accumulate silently if Python support is rarely enabled

### Low Priority

#### 5. Packaging and dependency boundaries could be clearer

Today the repo implicitly mixes:

- core C++ build dependencies
- optional Python binding dependencies
- script runtime dependencies

This is manageable, but the documentation should explicitly distinguish:

- core ReSeq requirements
- optional Python-binding requirements
- plotting-script requirements like `matplotlib` and `numpy`

## Recommendations

### Recommendation 1: Keep Python support, but narrow the modernization scope

Do not try to redesign the Python layer during SeqAn work. Treat Python modernization as a separate, bounded cleanup:

1. make the scripts work on modern Python
2. modernize the Python CMake/SWIG build
3. add smoke tests
4. decide whether the SWIG module is still worth keeping

This gives value quickly without mixing too many concerns.

### Recommendation 2: Fix the scripts first

Start with the helper scripts because they are the lowest-risk, highest-certainty improvements.

Priority changes:

- replace `time.clock()` with `time.perf_counter()` in [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py)
- replace `getopt` with `argparse` in [python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py)
- centralize file opening logic so `.gz` vs plain text handling is explicit and testable

### Recommendation 3: Modernize the CMake Python integration

The Python build should be brought in line with the repository’s target-based dependency direction.

Concrete changes:

- replace `find_program(PYTHON ...)` + `PythonLibs` with `find_package(Python3 COMPONENTS Interpreter Development REQUIRED)`
- replace `SWIG_ADD_MODULE` with `swig_add_library`
- link the generated module against project-level wrapper targets, not raw `${Boost_LIBRARIES}` / `${SEQAN_LIBRARIES}`
- keep Python support behind `RESEQ_BUILD_PYTHON`

### Recommendation 4: Add minimum verification if Python support is kept

If the repository continues shipping Python support, it should have at least:

- one import smoke test for the generated `DataStats` module
- one CLI smoke test for `reseq-prepare-names.py`
- one script launch test for `plotDataStats.py --help` or equivalent

Without this, the Python path will continue to drift.

### Recommendation 5: Reassess whether the SWIG module is still worth maintaining

The scripts are clearly useful. The SWIG bindings are less obviously justified.

A pragmatic decision point:

- if the `DataStats` Python module is actively used, modernize and test it
- if it is not actively used, deprecate or remove it and keep only the helper scripts

That decision should be based on actual project usage, not sentiment.

## Proposed Modernization Backlog

### Phase 1: Compatibility Fixes

- replace `time.clock()` in [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py)
- run the plotting script under a modern Python interpreter
- add at least one smoke test for both scripts

### Phase 2: Build-System Cleanup

- migrate to modern `Python3` CMake discovery
- migrate to modern `UseSWIG` commands
- switch Python linking to project wrapper targets

### Phase 3: Maintenance Decision

- confirm whether the SWIG bindings are still used
- either:
  - keep and test them properly
  - or deprecate/remove them

## Suggested Policy

Use this policy going forward:

- Python is an **optional feature**, not a core dependency of ReSeq
- helper scripts should remain lightweight and compatible with current Python
- the SWIG binding must either be tested in CI or be removed
- Python build configuration should follow the same target-based dependency rules as the C++ build

## Conclusion

Yes, there is meaningful Python modernization work to do.

The most urgent issue is not architectural, it is basic compatibility: [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py) still depends on an API removed in Python 3.8. After that, the next best investment is modernizing [python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt) so Python support is built with current CMake and linked through proper targets.

The scripts are worth keeping. The SWIG module should be kept only if it is actually used and tested.

## Sources

- Python `time.clock()` removal: https://docs.python.org/3/whatsnew/3.8.html
- Python `argparse`: https://docs.python.org/3/library/argparse.html
- CMake `FindPython3`: https://cmake.org/cmake/help/latest/module/FindPython3.html
- CMake `UseSWIG`: https://cmake.org/cmake/help/latest/module/UseSWIG.html
