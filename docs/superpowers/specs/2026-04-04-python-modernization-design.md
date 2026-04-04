# Python Modernization Design

**Date:** 2026-04-04
**Status:** Approved
**Input:** [docs/superpowers/specs/2026-04-04-python-modernization-report.md](/home/bernt-popp/development/ReSeq/docs/superpowers/specs/2026-04-04-python-modernization-report.md)

## Goal

Modernize ReSeq's optional Python surface so it works with current Python and current CMake/SWIG practices without changing the project's core dependency model.

This work keeps Python optional behind `RESEQ_BUILD_PYTHON`, preserves the current helper-script entry points, and keeps the `DataStats` SWIG module as a supported feature with build and import verification.

## Scope

In scope:

- modernize [python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt) to use modern Python3 discovery and current `UseSWIG` APIs
- keep the SWIG `DataStats` module buildable and importable as an optional feature
- replace removed Python APIs in [python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py)
- modernize [python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py) to current CLI/file-handling practices
- add automated test coverage for scripts and the SWIG module

Out of scope:

- packaging ReSeq as a Python package
- redesigning the `DataStats` Python API
- changing the installed script names or their basic stdout-oriented behavior

## Architecture

The Python layer remains split into two responsibilities:

- standalone helper scripts under `python/`
- an optional SWIG-generated `DataStats` extension module

The main build continues to gate Python support with `RESEQ_BUILD_PYTHON`. When that option is enabled, `python/CMakeLists.txt` becomes responsible for:

- discovering `Python3::Interpreter` and `Python3::Development`
- generating the SWIG wrapper with `swig_add_library`
- linking the module through repository targets instead of raw legacy library variables where possible
- registering Python smoke tests with CTest

The generated Python artifacts continue to live in the build tree under `build/pyMods` so the existing script lookup behavior remains intact.

## Components

### Python Build Integration

[python/CMakeLists.txt](/home/bernt-popp/development/ReSeq/python/CMakeLists.txt) will move from `find_program(PYTHON)` and `PythonLibs` to `find_package(Python3 COMPONENTS Interpreter Development REQUIRED)`.

The SWIG build will move from `SWIG_ADD_MODULE` and repeated `SWIG_LINK_LIBRARIES` calls to `swig_add_library(...)` and target-based `target_link_libraries(...)`.

The SWIG target should reuse the existing project implementation through `reseq_lib`, while still linking Python runtime requirements through the imported Python target. The output layout should continue to place `DataStats.py` and the compiled module in predictable build/install locations used by the scripts.

### Helper Scripts

[python/plotDataStats.py](/home/bernt-popp/development/ReSeq/python/plotDataStats.py) keeps its current command-line interface and plotting role, but updates removed timing APIs to modern replacements.

[python/reseq-prepare-names.py](/home/bernt-popp/development/ReSeq/python/reseq-prepare-names.py) keeps the same positional-argument workflow and stdout output, but moves from `getopt` to `argparse` and uses one explicit helper for opening plain text versus gzip-compressed files.

### Test Surface

Python verification will be integrated through CTest so it can run alongside the rest of the repository's tests. The minimum supported Python checks are:

- CLI smoke test for `reseq-prepare-names.py --help`
- CLI smoke test for `plotDataStats.py --help` or another shallow invocation that does not require real data
- import smoke test for the generated `DataStats` module from the build output directory when `RESEQ_BUILD_PYTHON=ON`
- behavior tests for `reseq-prepare-names.py`, including plain text and gzip input handling

## Data Flow

`reseq-prepare-names.py` reads the first record header from both FASTQ inputs to infer the shared name prefix, then streams the first file to stdout with rewritten names. That behavior stays the same; only the CLI parsing and file-opening boundary become more explicit and testable.

`plotDataStats.py` continues to locate `DataStats` from `RESEQ_PYMODS`, `build/pyMods`, or installed `lib` paths. The modernization preserves that search order.

The SWIG module continues to expose `DataStats` directly to Python. The new build wiring changes how the module is configured, not what it exports.

## Error Handling

- If `RESEQ_BUILD_PYTHON=ON` and Python development headers or SWIG are missing, configure/build should fail clearly.
- Script argument errors should come from standard `argparse` usage errors instead of custom `getopt` branches.
- Plain text and gzip opening should be explicit and centralized so unsupported inputs fail at a single boundary.
- The SWIG import smoke test should fail fast when the built module cannot be loaded from the build tree.

## Testing Strategy

Testing is part of the modernization, not a follow-up:

- add focused Python tests for script behavior
- add CTest-registered smoke tests for Python CLI/help flows
- add a build-tree import smoke test for the SWIG module
- keep tests shallow enough that they do not require external datasets or installed Python packages beyond the optional plotting runtime for plotting-specific checks

## Acceptance Criteria

- `plotDataStats.py` no longer depends on `time.clock()`
- `reseq-prepare-names.py` uses modern CLI parsing and has automated behavior coverage
- `RESEQ_BUILD_PYTHON=ON` configures through modern `Python3` CMake discovery
- the SWIG module is built through modern `UseSWIG` APIs and can be imported from the build tree
- Python smoke/behavior tests are runnable via CTest
