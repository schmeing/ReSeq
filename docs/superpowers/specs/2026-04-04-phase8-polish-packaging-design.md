# Phase 8: Polish & Packaging — Design Spec

**Goal:** Final modernization pass — C++20 language cleanup, type extraction, Python packaging, documentation refresh, and Vect container hardening.

**Prerequisite:** Phase 7 merged, SeqAn de-vendored. `make build && make test` passes on master.

---

## Current State

- **48 `typedef`** declarations in `utilities.hpp` (should be `using`)
- **29 `<stdint.h>`** includes across production source (should be `<cstdint>`)
- **29 `static const`** declarations that could be `constexpr`
- **618-line `utilities.hpp`** combining type aliases, helper functions, and SeqAn compat
- **Python bindings** OFF by default, `plotDataStats.py` has `sys.path` hacks, no `pyproject.toml`
- **Skewer** vendored with undocumented local modifications
- **README/CLAUDE.md** partially updated but not reflecting final architecture
- **Vect** has no debug bounds checking, offset semantics undocumented

---

## Sub-tasks

### 8a: C++ Language Cleanup

Mechanical replacements across ~30 production source files:

| Change | Count | Files |
|--------|-------|-------|
| `<stdint.h>` → `<cstdint>` | ~29 | All `.cpp`/`.h`/`.hpp` in `reseq/` |
| `typedef X Y` → `using Y = X` | ~48 | `utilities.hpp` + scattered |
| `static const` → `constexpr` | ~29 | Stats classes, BiasCalculationVectors, etc. |
| `[[nodiscard]]` | ~10 | `Load()`, `Save()`, `ReadBam()`, `IsValidRecord()`, etc. |
| Structured bindings | Opportunistic | Range-for with pairs, structured returns |

**`[[nodiscard]]` policy:** Apply to functions where ignoring the return value is always a bug — primarily bool-returning I/O functions (`Load`, `Save`, `ReadBam`) and validation functions (`IsValidRecord`). Do not apply to getters or functions with legitimate fire-and-forget usage.

**Structured bindings policy:** Use where they improve readability (e.g., iterating over `std::pair`, unpacking `std::tuple` returns). Do not force them where the existing code is already clear.

---

### 8b: Extract types.hpp

Extract from `utilities.hpp` into `reseq/types.hpp`:

- All 48 type aliases (`intQualDiff`, `uintBaseCall`, `uintQual`, `uintTempSeq`, etc.)
- `VectorAtomic<T>` struct
- SeqAn namespace compat alias (`namespace seqan = seqan2`)
- SeqAn includes needed by the above (`<seqan/version.h>`, `<seqan/bam_io.h>`, `<seqan/modifier.h>`, `<seqan/sequence.h>`)
- The `SWIG` guard (`#ifndef SWIG`) for the SeqAn-dependent portion

`utilities.hpp` retains:
- Helper functions (`SafePercent`, `Divide`, `SetToMin`, `SetToMax`, etc.)
- `DominantBase`, `DominantBaseWithMemory` classes
- `Complement` struct
- `KmerCount` template
- `FragmentSite` struct

`utilities.hpp` adds `#include "types.hpp"` at the top for full backward compatibility — no other files need their includes updated.

---

### 8c: Python Packaging

**Goal:** Make the Python plotting tools installable as a proper package.

- Create `pyproject.toml` at repository root with project metadata, py39+ target, ruff config
- Create `python/reseq/__init__.py` package structure
- Move `plotDataStats.py` and `reseq-prepare-names.py` into `python/reseq/` as modules
- Remove `sys.path` hacks from `plotDataStats.py`
- Update `python/CMakeLists.txt` to install the SWIG module into the package directory
- Update `python/DataStats.i` path references if needed
- Keep `RESEQ_BUILD_PYTHON=OFF` as default
- Update `conda/reseq/meta.yaml` to reflect new package structure

---

### 8d: Skewer MODIFICATIONS.md

Create `skewer/MODIFICATIONS.md`:

```markdown
# Skewer Local Modifications

Vendored from: [relipmoc/skewer](https://github.com/relipmoc/skewer) v0.2.2 (2016)
License: MIT

## Why vendored

ReSeq uses only the core bit-masked k-difference adapter matching algorithm
(~1,200 lines from matrix.cpp/h and fastq.cpp/h). The vendored copy carries
local modifications that are not suitable for upstreaming.

## Modifications

1. **Namespace wrapping** — All code wrapped in `namespace skewer {}` to avoid
   global namespace pollution when linking as a static library.

2. **Removed unused code** (~4,700 lines) — Deleted `main.cpp` (CLI tool),
   `parameter.cpp/h` (argument parsing), and all code not used by ReSeq's
   `AdapterStats::Detect()` pathway.

3. **CMakeLists.txt** — Created minimal CMake build config producing a static
   library (`skewer_matrix`) from `matrix.cpp` and `fastq.cpp`.

4. **GCC 15 compatibility** — Added `const` qualifier to
   `ElementComparator::operator()` for C++20 compliance.
```

---

### 8e: Documentation Refresh

**README.md** — Final pass to reflect completed refactoring:
- Build requirements: C++20, CMake 3.16+, GCC 10+/Clang 12+
- Dependencies: SeqAn 2.5.2 (FetchContent), GoogleTest (FetchContent), NLopt (FetchContent), Boost, ZLIB, BZip2
- Skewer remains vendored (with rationale)
- Test commands: `make test`, `make coverage`, GoogleTest filter syntax
- Codecov badge (already present)
- Python bindings section updated for new package structure

**CLAUDE.md** — Final pass:
- Remove references to `2016-05-15_ROOTPWA/` (absorbed in Phase 3)
- Update dependency section for external SeqAn model
- Update architecture section for decomposed classes (ReadSequenceStats, BamIngestionEngine, CLI modules)
- Reference `cmake/ReSeqDependencies.cmake` for dependency policy

---

### 8f: Vect Debug Bounds Checking

Add debug-only bounds checking to `Vect<T>`:

- In non-const `operator[]`: add `assert(n >= from() || vec_.second.empty())` before `Ensure(n)` to catch clearly invalid access patterns in debug builds
- In const `operator[]`: document the silent-default behavior (returns `dummy_` for out-of-range)
- Add class-level documentation comment explaining:
  - Vect is a vector with an offset (`from()`) — element `n` is stored at `vec_.second[n - from()]`
  - `operator[]` auto-extends the vector (non-const) or returns zero (const)
  - `at()` throws for out-of-range access
  - `from()` / `to()` give the valid index range

No functional changes — assertions only fire in debug builds (`-DCMAKE_BUILD_TYPE=Debug`).

---

## Verification Gate

After all sub-tasks complete:

- [ ] `make build` succeeds (GCC and Clang)
- [ ] `make test` passes all tests
- [ ] `make format-check` clean
- [ ] `make lint` passes (clang-tidy clean)
- [ ] `pre-commit run --all-files` passes
- [ ] Regression tests pass (`build/bin/reseq_test --gtest_filter="RegressionTest.*"`)
- [ ] Python bindings build and import (`-DRESEQ_BUILD_PYTHON=ON`)
- [ ] No `<stdint.h>` in production source: `grep -rn '<stdint.h>' reseq/*.cpp reseq/*.h reseq/*.hpp` empty
- [ ] No `typedef` in `utilities.hpp`/`types.hpp`: `grep -c 'typedef' reseq/types.hpp reseq/utilities.hpp` = 0

---

## Files Created/Modified

| File | Action |
|------|--------|
| `reseq/types.hpp` | New — type aliases, VectorAtomic, SeqAn compat |
| `reseq/utilities.hpp` | Slim down, include types.hpp |
| `reseq/*.cpp`, `reseq/*.h`, `reseq/*.hpp` | stdint.h→cstdint, typedef→using, static const→constexpr, [[nodiscard]], structured bindings |
| `reseq/Vect.hpp` | Debug bounds checking, documentation |
| `pyproject.toml` | New — Python project metadata |
| `python/reseq/__init__.py` | New — package init |
| `python/reseq/plotDataStats.py` | Moved + cleaned up |
| `python/reseq/reseq_prepare_names.py` | Moved + cleaned up |
| `python/CMakeLists.txt` | Updated paths |
| `python/DataStats.i` | Updated paths if needed |
| `conda/reseq/meta.yaml` | Updated for new package structure |
| `skewer/MODIFICATIONS.md` | New — documents local patches |
| `README.md` | Final refresh |
| `CLAUDE.md` | Final refresh |

---

## What This Phase Does NOT Do

- No god object decomposition (Phase 5c-5f permanently deferred)
- No Vect functional changes (debug assertions only)
- No SeqAn API migration (namespace alias handles it)
- No new test classes (test infrastructure complete from Phase 7)
