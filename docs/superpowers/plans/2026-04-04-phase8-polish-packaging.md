# Phase 8: Polish & Packaging — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Final modernization: C++20 language cleanup, type extraction, Python packaging, Vect hardening, and documentation.

**Architecture:** Seven sequential tasks. Tasks 1-2 are C++ modernization (language cleanup then type extraction). Task 3 is Vect hardening. Task 4 is Python packaging. Task 5 is Skewer docs. Task 6 is documentation refresh. Task 7 is final verification. Each task independently committable.

**Tech Stack:** C++20, CMake 3.16+, GoogleTest, Python 3.9+, SWIG, ruff

**Design spec:** `docs/superpowers/specs/2026-04-04-phase8-polish-packaging-design.md`

**Prerequisite:** Phase 7 and SeqAn de-vendoring merged. `make build && make test` passes on master.

---

## Critical Invariants

1. **All existing tests must pass.** `make build && make test` after every task.
2. **Conventional commits** with `(8)` scope tag and sub-task letter.
3. **No behavioral changes** — this is purely mechanical cleanup and packaging.

---

## Task 1: C++ Language Cleanup (8a)

**Goal:** Replace `<stdint.h>` with `<cstdint>`, `typedef` with `using`, `static const` with `constexpr`, add `[[nodiscard]]`, and use structured bindings where appropriate.

**Files modified:** ~30 files in `reseq/`

### Step-by-step

- [ ] **Step 1.1: Replace `<stdint.h>` with `<cstdint>` across all files**

Run this replacement across all production and test source files in `reseq/`:

```bash
sed -i 's/#include <stdint.h>/#include <cstdint>/' reseq/*.cpp reseq/*.h reseq/*.hpp
```

Verify no `<stdint.h>` remains:

```bash
grep -rn '<stdint.h>' reseq/
```

Expected: empty output.

- [ ] **Step 1.2: Convert `typedef` to `using` in `utilities.hpp`**

In `reseq/utilities.hpp`, replace all `typedef` declarations (lines 39-89) with `using` syntax. The pattern is:

```cpp
// Before:
typedef int8_t intQualDiff;
// After:
using intQualDiff = int8_t;
```

Also convert the SeqAn type aliases in the `namespace utilities` block (lines 94-110):

```cpp
// Before:
typedef seqan::String<seqan::CigarElement<>> CigarString;
// After:
using CigarString = seqan::String<seqan::CigarElement<>>;
```

There are approximately 48 typedef declarations to convert. Use `sed` for the simple cases:

```bash
cd reseq && sed -i -E 's/^typedef ([^ ]+) ([^;]+);/using \2 = \1;/' utilities.hpp
```

Then manually fix the multi-line and `const` typedefs (lines 97-110) which the sed won't handle correctly. These have the pattern `typedef const seqan::ModifiedString<...> Name;` which should become `using Name = const seqan::ModifiedString<...>;`.

- [ ] **Step 1.3: Convert `typedef` in other files**

Search for remaining typedefs in other production files:

```bash
grep -rn "typedef " reseq/*.h reseq/*.hpp reseq/*.cpp | grep -v Test | grep -v "seqan/" | grep -v "//"
```

Convert any found to `using` syntax.

- [ ] **Step 1.4: Convert `static const` to `static constexpr`**

In the following files, replace `static const` with `static constexpr` for integral and bool constants:

- `reseq/AdapterStats.h:28` — `kKmerLength`
- `reseq/CoverageStats.h:136-137` — `kBlockSize`, `kMaxCoverage`
- `reseq/FragmentDistributionStats.h:70-99,250,315` — `kMaxLikelihoodCalculations`, `kSplinePrecisionFactor`, `kGCSplineDf`, `kMaxKnotShift`, `kPercentGCSitesForNormalization`, `kNumFitsInsertLength`, `kSurMult`, `kGCExp`, `kMaxDuplications`, `kDispersionBinSize`, `kMaxChi2Calculations`, `kMaxBinsQueuedForBiasCalc`
- `reseq/FragmentDuplicationStats.h:22` — `kMaxDuplication`
- `reseq/ProbabilityEstimates.h:57,110,375,631` — `kNumMargins` (4 occurrences in template classes)
- `reseq/QualityStats.h` — `kSqFragmentLengthBinSize`, `kWriteOut4dMatrixCsvs`

Use:

```bash
sed -i 's/static const \(uint\|int\|bool\)/static constexpr \1/' reseq/AdapterStats.h reseq/CoverageStats.h reseq/FragmentDistributionStats.h reseq/FragmentDuplicationStats.h reseq/ProbabilityEstimates.h reseq/QualityStats.h
```

Do NOT convert `static constexpr const char*` patterns (e.g., `kAdapterSearchInfoFile`) — `constexpr` pointer-to-const requires `constexpr const char* const` which changes semantics.

Do NOT convert `static const T dummy_` in `Vect.hpp` — it requires out-of-line definition.

- [ ] **Step 1.5: Add `[[nodiscard]]` to key bool-returning functions**

In `reseq/DataStats.h`, add `[[nodiscard]]` to:

```cpp
[[nodiscard]] bool IsValidRecord(const seqan::BamAlignmentRecord& record);
[[nodiscard]] bool ReadBam(const char* bam_file, ...);
[[nodiscard]] bool Load(const char* archive_file);
[[nodiscard]] bool Save(const char* archive_file, bool text_format = false) const;
```

In `reseq/DataStatsInterface.h`:

```cpp
[[nodiscard]] bool Load(const char* archive_file);
```

In `reseq/AdapterStats.h`, add to `LoadAdapters`, `Detect`, `PredictAdapters`:

```cpp
[[nodiscard]] bool LoadAdapters(const char* adapter_file, const char* adapter_matrix);
[[nodiscard]] bool PredictAdapters();
[[nodiscard]] bool Detect(...);
```

In `reseq/Reference.h`, add to `Load`, `ReplaceN`, and any other bool-returning public methods.

In `reseq/BamIngestionEngine.h`:

```cpp
[[nodiscard]] bool Run(...);
```

- [ ] **Step 1.6: Use structured bindings where they improve readability**

Search for range-for loops iterating over containers of pairs or maps:

```bash
grep -rn "for.*auto&.*:.*\." reseq/*.cpp | grep -v Test | head -20
```

Convert clear cases where `.first` / `.second` are used immediately. Example pattern:

```cpp
// Before:
for (auto& entry : map) {
    use(entry.first, entry.second);
}
// After:
for (auto& [key, value] : map) {
    use(key, value);
}
```

Only convert where the result is unambiguously clearer. Skip any case where the pair members are used many lines apart or where the loop body is complex.

- [ ] **Step 1.7: Build and test**

```bash
make build && make test
```

- [ ] **Step 1.8: Format and commit**

```bash
make format
git add -A reseq/
git commit -m "style(8a): C++20 language cleanup

Replace <stdint.h> with <cstdint> across 29 files.
Convert 48 typedef to using declarations.
Convert static const to static constexpr for integral constants.
Add [[nodiscard]] to bool-returning I/O and validation functions.
Use structured bindings where they improve readability."
```

---

## Task 2: Extract types.hpp from utilities.hpp (8b)

**Goal:** Extract type aliases and VectorAtomic into a standalone header for reuse.

**Files created:**
- `reseq/types.hpp`

**Files modified:**
- `reseq/utilities.hpp`

### Step-by-step

- [ ] **Step 2.1: Create `reseq/types.hpp`**

Extract from `utilities.hpp` into a new file `reseq/types.hpp`:

1. The include guard and `<cstdint>` include
2. The SeqAn includes + namespace compat alias (lines 18-26)
3. The `UNUSED` macro (lines 31-35)
4. All type aliases in `namespace reseq` (lines 38-89, now `using` declarations)
5. The `#ifndef SWIG` block with SeqAn type aliases (lines 91-110)
6. The `VectorAtomic<T>` struct (around line 349 in current file)

The header should be self-contained — include everything it needs.

- [ ] **Step 2.2: Update `utilities.hpp` to include `types.hpp`**

Replace the extracted sections in `utilities.hpp` with:

```cpp
#include "types.hpp"
```

Remove the now-redundant includes (`<cstdint>`, SeqAn headers, `CMakeConfig.h`) that are provided by `types.hpp`. Keep everything else in `utilities.hpp` (helper functions, classes, `UNUSED` macro if used outside types).

- [ ] **Step 2.3: Build and test**

```bash
make build && make test
```

No other files should need changes — they include `utilities.hpp` which transitively provides everything from `types.hpp`.

- [ ] **Step 2.4: Format and commit**

```bash
make format
git add reseq/types.hpp reseq/utilities.hpp
git commit -m "refactor(8b): extract types.hpp from utilities.hpp

Move 48 type aliases, VectorAtomic<T>, and SeqAn compat alias
into standalone types.hpp. utilities.hpp includes types.hpp for
backward compatibility — no consumer changes needed."
```

---

## Task 3: Vect Debug Bounds Checking + Documentation (8f)

**Goal:** Add debug assertions and documentation to Vect.

**Files modified:**
- `reseq/Vect.hpp`

### Step-by-step

- [ ] **Step 3.1: Add class-level documentation**

Add a documentation comment before the class definition in `reseq/Vect.hpp`:

```cpp
/// Offset vector: a std::vector with a non-zero starting index.
///
/// Vect<T> stores a pair (offset, data) where element n is at data[n - offset].
/// Key semantics:
///   - from() returns the index of the first element (the offset)
///   - to() returns one past the index of the last element
///   - operator[] (non-const) auto-extends the vector to accommodate the index
///   - operator[] (const) returns a zero-valued dummy for out-of-range access
///   - at() throws std::out_of_range for invalid indices
///   - Shrink() removes leading/trailing zeros and frees unused capacity
template <typename T>
class Vect {
```

- [ ] **Step 3.2: Add debug assertion to non-const operator[]**

Add `#include <cassert>` at the top of Vect.hpp. Then in the non-const `operator[]`, add an assertion before `Ensure()`:

```cpp
T& operator[](typename std::vector<T>::size_type n) {
    assert(n < std::numeric_limits<typename std::vector<T>::size_type>::max() - 1024 &&
           "Vect::operator[] called with suspiciously large index");
    Ensure(n);
    return vec_.second[n - vec_.first];
}
```

This catches clearly invalid indices (e.g., negative values cast to size_t producing huge numbers) in debug builds without affecting release performance.

- [ ] **Step 3.3: Document const operator[] behavior**

Add a comment to the const `operator[]`:

```cpp
/// Returns element at index n, or a zero-valued dummy if n is out of range.
/// This silent-default behavior is intentional — Vect is used in contexts
/// where out-of-range reads should return zero (e.g., sparse histograms).
const T& operator[](typename std::vector<T>::size_type n) const {
```

- [ ] **Step 3.4: Build and test**

```bash
make build && make test
```

- [ ] **Step 3.5: Format and commit**

```bash
make format
git add reseq/Vect.hpp
git commit -m "refactor(8f): add Vect documentation and debug bounds checking

Add class-level documentation explaining offset semantics.
Add debug assertion in non-const operator[] to catch invalid indices.
Document const operator[] silent-default behavior."
```

---

## Task 4: Python Packaging (8c)

**Goal:** Restructure Python code as a proper installable package.

**Files created:**
- `python/reseq/__init__.py`
- `python/reseq/plot_data_stats.py` (moved from `plotDataStats.py`)
- `python/reseq/prepare_names.py` (moved from `reseq-prepare-names.py`)

**Files modified:**
- `pyproject.toml`
- `python/CMakeLists.txt`
- `conda/reseq/meta.yaml`

### Step-by-step

- [ ] **Step 4.1: Create package directory and move scripts**

```bash
mkdir -p python/reseq
touch python/reseq/__init__.py
cp python/plotDataStats.py python/reseq/plot_data_stats.py
cp python/reseq-prepare-names.py python/reseq/prepare_names.py
```

Keep the originals in place for now (remove in a later step after wiring is verified).

- [ ] **Step 4.2: Remove sys.path hacks from `python/reseq/plot_data_stats.py`**

In `python/reseq/plot_data_stats.py`, replace the sys.path hack block (lines 20-26):

```python
if "RESEQ_PYMODS" in os.environ and os.path.isdir(os.path.realpath(os.environ["RESEQ_PYMODS"])):
    sys.path.append(os.path.realpath(os.environ["RESEQ_PYMODS"]))
elif os.path.isdir(os.path.dirname(os.path.realpath(__file__)) + "/../build/pyMods/"):
    sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../build/pyMods/")
elif os.path.isdir(os.path.dirname(os.path.realpath(__file__)) + "/../lib/"):
    sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../lib/")
import DataStats
```

with:

```python
try:
    import DataStats
except ImportError:
    # Fall back to common build locations when SWIG module is not installed
    _script_dir = os.path.dirname(os.path.realpath(__file__))
    for _candidate in [
        os.environ.get("RESEQ_PYMODS", ""),
        os.path.join(_script_dir, "..", "..", "build", "pyMods"),
        os.path.join(_script_dir, "..", "..", "lib"),
    ]:
        if _candidate and os.path.isdir(_candidate):
            sys.path.append(os.path.realpath(_candidate))
            try:
                import DataStats
                break
            except ImportError:
                sys.path.pop()
                continue
    else:
        raise ImportError(
            "Cannot find DataStats SWIG module. Build with -DRESEQ_BUILD_PYTHON=ON "
            "or set RESEQ_PYMODS to the module directory."
        )
```

- [ ] **Step 4.3: Update pyproject.toml with project metadata**

Add project metadata to `pyproject.toml`:

```toml
[project]
name = "reseq"
version = "1.1.0"
description = "Plotting tools for ReSeq sequencing simulator"
requires-python = ">=3.9"
license = {text = "MIT"}
dependencies = [
    "matplotlib",
    "numpy",
]

[project.scripts]
reseq-plot = "reseq.plot_data_stats:main"
reseq-prepare-names = "reseq.prepare_names:main"

[build-system]
requires = ["setuptools>=64"]
build-backend = "setuptools.backends._legacy:_Backend"

[tool.ruff]
target-version = "py39"
line-length = 120

[tool.ruff.lint]
select = ["E", "F", "W", "I", "B", "SIM"]

[tool.ruff.lint.per-file-ignores]
"python/reseq/plot_data_stats.py" = ["E501", "SIM105"]
"python/plotDataStats.py" = ["E501", "SIM105"]

[tool.mypy]
python_version = "3.9"
check_untyped_defs = true
ignore_missing_imports = true
warn_redundant_casts = true
warn_unused_ignores = true
no_implicit_optional = true
exclude_gitignore = true
files = ["python/reseq", "python/tests"]
```

- [ ] **Step 4.4: Add `main()` entry points to the moved scripts**

In `python/reseq/plot_data_stats.py`, wrap the existing `if __name__ == "__main__"` block into a `main()` function:

```python
def main():
    # existing argparse + plotting code
    ...

if __name__ == "__main__":
    main()
```

Do the same for `python/reseq/prepare_names.py`.

- [ ] **Step 4.5: Update `python/reseq/__init__.py`**

```python
"""ReSeq plotting and preparation tools."""
```

- [ ] **Step 4.6: Update `python/CMakeLists.txt` install paths**

Update the install section to include the new package:

```cmake
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/python/reseq/"
  DESTINATION lib/python/reseq
  FILES_MATCHING PATTERN "*.py"
)
```

Keep the existing install targets for backward compatibility.

- [ ] **Step 4.7: Update ruff config paths**

In the Makefile, ensure `make format` and `make lint` cover `python/reseq/` in addition to `python/`.

- [ ] **Step 4.8: Remove original script files**

```bash
git rm python/plotDataStats.py python/reseq-prepare-names.py
```

Update any references in `python/DataStats.i`, test files, or CI to point to new locations.

- [ ] **Step 4.9: Build and test**

```bash
make build && make test
make format-check
```

If Python bindings are testable: `cmake -S . -B build -DRESEQ_BUILD_PYTHON=ON && cmake --build build -j$(nproc)`

- [ ] **Step 4.10: Format and commit**

```bash
make format
git add -A python/ pyproject.toml conda/reseq/meta.yaml
git commit -m "feat(8c): restructure Python as installable package

Create python/reseq/ package with proper __init__.py.
Move plotDataStats.py → reseq/plot_data_stats.py with main() entry point.
Move reseq-prepare-names.py → reseq/prepare_names.py.
Replace sys.path hacks with try/except import fallback.
Update pyproject.toml with project metadata, py39+ target, entry points."
```

---

## Task 5: Skewer MODIFICATIONS.md (8d)

**Goal:** Document the local modifications to the vendored skewer code.

**Files created:**
- `skewer/MODIFICATIONS.md`

### Step-by-step

- [ ] **Step 5.1: Create `skewer/MODIFICATIONS.md`**

```markdown
# Skewer Local Modifications

Vendored from: [relipmoc/skewer](https://github.com/relipmoc/skewer) v0.2.2 (2016)
License: MIT

## Why vendored

ReSeq uses only the core bit-masked k-difference adapter matching algorithm
(~1,200 lines from `matrix.cpp/h` and `fastq.cpp/h`). The vendored copy
carries local modifications that are not suitable for upstreaming.

## Modifications

1. **Namespace wrapping** — All code wrapped in `namespace skewer {}` to
   avoid global namespace pollution when linking as a static library.

2. **Removed unused code** (~4,700 lines) — Deleted `main.cpp` (CLI tool),
   `parameter.cpp/h` (argument parsing), and all code not used by ReSeq's
   `AdapterStats::Detect()` pathway.

3. **CMakeLists.txt** — Created minimal CMake build config producing a
   static library (`skewer_matrix`) from `matrix.cpp` and `fastq.cpp`.

4. **GCC 15 compatibility** — Added `const` qualifier to
   `ElementComparator::operator()` for C++20 compliance.
```

- [ ] **Step 5.2: Commit**

```bash
git add skewer/MODIFICATIONS.md
git commit -m "docs(8d): document skewer vendored modifications"
```

---

## Task 6: Documentation Refresh (8e)

**Goal:** Final pass on README.md and CLAUDE.md to reflect completed refactoring.

**Files modified:**
- `README.md`
- `CLAUDE.md`

### Step-by-step

- [ ] **Step 6.1: Update README.md**

Refresh the Requirements section to reflect:
- C++20 compiler (GCC 10+, Clang 12+)
- CMake 3.16+
- Boost 1.48+ (serialization, program_options, filesystem, system, math)
- ZLIB, BZip2
- SeqAn 2.5.2, GoogleTest, NLopt fetched automatically via FetchContent
- Python 3.9+ for optional plotting tools (requires SWIG 3+ and python3-dev)

Update the Installation section to show:
```bash
cmake -S . -B build
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```

Update the skewer attribution to note it carries local modifications (see `skewer/MODIFICATIONS.md`).

Remove any remaining references to vendored SeqAn.

- [ ] **Step 6.2: Update CLAUDE.md**

Final pass to reflect:
- External SeqAn 2.5.2 dependency model (already partially updated)
- `cmake/ReSeqDependencies.cmake` for dependency policy
- New `reseq/types.hpp` header
- Python package structure (`python/reseq/`)
- Completed Phase 1-8 status
- Skewer vendored with `MODIFICATIONS.md`

- [ ] **Step 6.3: Build and test**

```bash
make build && make test
```

- [ ] **Step 6.4: Format and commit**

```bash
make format
git add README.md CLAUDE.md
git commit -m "docs(8e): final documentation refresh

Update README.md with C++20 requirements, FetchContent dependencies,
new test commands, and Python package structure.
Update CLAUDE.md with final architecture after Phases 1-8."
```

---

## Task 7: Final Verification (8g)

**Goal:** Run the complete verification gate.

### Step-by-step

- [ ] **Step 7.1: Full build and test**

```bash
make build && make test
```

- [ ] **Step 7.2: Format and lint**

```bash
make format-check
make lint
```

- [ ] **Step 7.3: Pre-commit hooks**

```bash
pre-commit run --all-files
```

- [ ] **Step 7.4: Regression tests**

```bash
build/bin/reseq_test --gtest_filter="RegressionTest.*"
```

- [ ] **Step 7.5: Verify cleanup targets**

```bash
grep -rn '<stdint.h>' reseq/*.cpp reseq/*.h reseq/*.hpp
# Expected: empty

grep -c 'typedef' reseq/types.hpp reseq/utilities.hpp
# Expected: 0 for both
```

- [ ] **Step 7.6: Git log review**

```bash
git log --oneline master..HEAD
```

Expected: clean conventional commits with `(8)` scope tags.

---

## Verification Checklist

- [ ] `make build` succeeds (GCC and Clang)
- [ ] `make test` passes all tests
- [ ] `make format-check` clean
- [ ] `make lint` passes
- [ ] `pre-commit run --all-files` passes
- [ ] Regression tests pass
- [ ] No `<stdint.h>` in production source
- [ ] No `typedef` in `types.hpp` or `utilities.hpp`
- [ ] `types.hpp` exists and is included by `utilities.hpp`
- [ ] `skewer/MODIFICATIONS.md` exists
- [ ] Python package structure: `python/reseq/__init__.py` exists
