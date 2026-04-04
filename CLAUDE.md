# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository

This is a fork (`berntpopp/ReSeq`) of upstream (`schmeing/ReSeq`). The upstream remote has been removed.

**All work stays in the fork.** Never add the upstream remote, create PRs on `schmeing/ReSeq`, or push to it. Use `--repo berntpopp/ReSeq` with all `gh` commands.

## Build and Test Commands

```bash
make build              # Configure (CMake) and build
make test               # Build and run all unit tests via ctest
make coverage           # Build with coverage, run tests, generate HTML report
make format             # Format C++ (clang-format) and Python (ruff) in-place
make format-check       # Dry-run format check (used in CI)
make lint               # Run clang-tidy and ruff
make clean              # Remove build directory
make changelog          # Generate CHANGELOG.md via git-cliff
pre-commit run --all-files  # Run all pre-commit hooks
```

Tests run via a separate `reseq_test` binary using CTest. GoogleTest filter syntax works: `build/bin/reseq_test --gtest_filter="SimulatorTest.*"`. Tests must run sequentially (not parallel) due to inter-test dependencies.

The build requires: C++20 compiler (GCC 10+, Clang 12+), CMake 3.16+, Boost 1.48+ (serialization, program_options, filesystem, system, math), ZLIB, BZip2. SeqAn 2.5.2, GoogleTest, and NLopt are fetched automatically via FetchContent (or found via `find_package()` if installed). Python bindings are OFF by default (`-DRESEQ_BUILD_PYTHON=ON` requires SWIG 3+ and python3-dev).

## Architecture

ReSeq is a bioinformatics tool that learns error/quality profiles from real Illumina paired-end sequencing data and uses them to simulate realistic reads.

### Commands (entry point: `reseq/main.cpp`)

| Command | Purpose |
|---------|---------|
| `illuminaPE` | Full pipeline: collect stats from BAM → estimate probabilities via IPF → simulate reads |
| `seqToIllumina` | Apply Illumina error/quality model to input sequences (no coverage model) |
| `queryProfile` | Extract info from `.reseq` stats files (fragment length bias, ref seq bias, etc.) |
| `replaceN` | Replace ambiguous bases (N) in reference sequences |

### Core Components (in `reseq/`)

One CMake static library + two executables:

**`reseq_lib`** (static) — all production source:
- `DataStats` — top-level aggregator that orchestrates all sub-stats
- `AdapterStats`, `CoverageStats`, `ErrorStats`, `FragmentDistributionStats`, `FragmentDuplicationStats`, `QualityStats`, `TileStats` — each models a specific aspect of sequencing
- `Reference` — reference genome loading, surrounding context, excluded regions
- `Surrounding` / `SurroundingBase` — sequence context modeling for error patterns
- `Vect` — offset vector (indexed from non-zero starting position), used pervasively
- `types.hpp` — type aliases, VectorAtomic, SeqAn compat (extracted from utilities.hpp)
- `utilities.hpp` — helper functions, classes (includes types.hpp)
- `ProbabilityEstimates` — Iterative Proportional Fitting (IPF) for multi-dimensional probability tables
- `Simulator` — block-based read simulation engine with threading support

**`reseq`** — thin CLI executable linking `reseq_lib`

**`reseq_test`** — test executable linking `reseq_lib` + GoogleTest

### External Dependencies

SeqAn 2.5.2 (bioinformatics, header-only) is resolved via `find_package()` or `FetchContent` — see `cmake/ReSeqDependencies.cmake`. GoogleTest and NLopt are also fetched via FetchContent. `skewer/` remains vendored (adapter trimming, carries local modifications — see `skewer/MODIFICATIONS.md`).

### Test Structure

Each component has a `*Test.cpp` / `*Test.h` pair inheriting from `BasicTestClass.hpp` (which extends `::testing::Test`). Test data lives in `test/` (E. coli and Drosophila references, BAMs, adapters). Tests are compiled into the `reseq_test` binary and run via CTest (`make test` or `ctest --output-on-failure`). GoogleTest filter syntax: `build/bin/reseq_test --gtest_filter="SimulatorTest.*"`.

### Versioning

Single source of truth: `VERSION` file (currently `1.1.0`). CMake reads it at configure time. `CMakeConfig.h.in` generates version macros including `RESEQ_GIT_VERSION` from `git describe`.

## Code Style

- C++: clang-format (LLVM-based, 120 column limit, 4-space indent, C++20)
- Python: ruff (py39 target, 120 line length, rules: E/F/W/I/B/SIM); package at `python/reseq/`
- Commits: conventional commits (`feat:`, `fix:`, `build:`, `style:`, `test:`, `ci:`)
- Pre-commit hooks enforce formatting on `reseq/` and `python/` only

## CI

GitHub Actions (`.github/workflows/ci.yml`): builds with GCC 13 and Clang 17 on Ubuntu 24.04, runs format checks, verifies version tags match `VERSION` file, runs ASan+UBSan sanitizer checks, and generates code coverage reports.
