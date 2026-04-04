# Phase 6: DRY & Interface Cleanup — Design Specification

**Date:** 2026-04-04
**Status:** Draft
**Scope:** Named constants, container type aliases, FRIEND_TEST removal, DataStatsInterface slimming
**Prerequisite:** Phase 5b (DataStats decomposition) merged. `make build && make test` passes on master.
**Parent spec:** `docs/superpowers/specs/2026-04-01-comprehensive-refactoring-design.md` (Phase 6)

---

## Scope Decisions

The original Phase 6 spec defines 5 sub-tasks (6a-6e). This implementation **defers 6a (StatVariable\<T\> registration)** because it touches serialization across all 8 Stats classes — high risk immediately after the Phase 5b serialization-preserving decomposition. The remaining 4 sub-tasks are lower risk and largely mechanical.

**In scope:** 6b, 6c, 6d, 6e
**Deferred:** 6a (can be revisited as a standalone phase)

---

## Execution Order

Bottom-up: foundations first, broadest change last.

| Step | Sub-task | What | Risk |
|------|----------|------|------|
| 1 | 6b | Named constants (`constants.hpp`) | LOW |
| 2 | 6c | Container type aliases (`container_types.hpp`) | LOW |
| 3 | 6e | FRIEND_TEST removal + gtest cleanup | LOW |
| 4 | 6d | Slim DataStatsInterface | MEDIUM |

Each step is merged to master before the next begins.

---

## Critical Invariants

1. **All existing tests pass** after each step.
2. **Binary serialization compatibility preserved** — no data layout changes, no field reordering.
3. **No public DataStats API changes** — DataStatsInterface changes are internal to its consumers.
4. **Conventional commits** with `(6)` scope tag and sub-task letter (e.g., `refactor(6b):`).

---

## Step 1: Named Constants (6b)

### Goal

Replace magic numbers with named domain constants. Improves readability and reduces error risk when the same value appears in array template parameters, loop bounds, and conditionals across 15+ files.

### Design

Create `reseq/constants.hpp`:

```cpp
#ifndef RESEQ_CONSTANTS_H
#define RESEQ_CONSTANTS_H

#include "utilities.hpp" // For type aliases (uintBaseCall, uintTempSeq, uintQual)

namespace reseq {

inline constexpr uintBaseCall kNumBases = 4;        // A, C, G, T
inline constexpr uintBaseCall kNumBasesN = 5;        // A, C, G, T, N
inline constexpr uintTempSeq kTemplateSegments = 2;  // first, second read
inline constexpr uintTempSeq kStrands = 2;           // forward, reverse
inline constexpr uintQual kPhredSangerOffset = 33;
inline constexpr uintQual kPhredIlluminaOffset = 64;
inline constexpr uint16_t kGCBins = 101;             // 0-100% inclusive

} // namespace reseq

#endif // RESEQ_CONSTANTS_H
```

### Replacement Rules

| Magic | Constant | Context | ~Count |
|-------|----------|---------|--------|
| `4` | `kNumBases` | Array template params, loop bounds for DNA bases | ~60 |
| `5` | `kNumBasesN` | Array template params, loop bounds for ACGTN | ~50 |
| `2` | `kTemplateSegments` | Array template params for first/second read | ~30 |
| `2` | `kStrands` | Array template params for forward/reverse | ~10 |
| `33` | `kPhredSangerOffset` | Quality offset comparison | ~3 |
| `64` | `kPhredIlluminaOffset` | Quality offset comparison | ~2 |
| `101` | `kGCBins` | GC content bin array sizes | ~12 |

**Not replaced:** Loop counters that count down (`for (auto base = 5; base--;)`) keep their numeric form — replacing with a constant would obscure the countdown pattern. Only array declarations and explicit comparisons are updated.

**Disambiguation:** Not every `4` or `5` is a base count. Each replacement must be verified in context. For example, `4` in `Surrounding::Length()` is unrelated to DNA bases.

### Verification

`make build && make test` — purely mechanical, no logic changes.

---

## Step 2: Container Type Aliases (6c)

### Goal

Replace deeply nested container declarations (up to 7 levels in ErrorStats/QualityStats) with readable type aliases. The aliases use constants from Step 1.

### Design

Create `reseq/container_types.hpp`:

```cpp
#ifndef RESEQ_CONTAINER_TYPES_H
#define RESEQ_CONTAINER_TYPES_H

#include <array>
#include <vector>

#include "constants.hpp"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq {

// Atomic accumulator vector (thread-safe, filled during multi-threaded BAM pass)
template <typename T>
using AtomicVec = std::vector<utilities::VectorAtomic<T>>;

// Domain-indexed array aliases (use constants for sizes)
template <typename T>
using PerSegment = std::array<T, kTemplateSegments>;

template <typename T>
using PerStrand = std::array<T, kStrands>;

template <typename T>
using PerBase = std::array<T, kNumBases>;

template <typename T>
using PerBaseN = std::array<T, kNumBasesN>;

template <typename T>
using PerSegmentPerBase = PerSegment<PerBase<T>>;

template <typename T>
using PerSegmentPerBaseN = PerSegment<PerBaseN<T>>;

template <typename T>
using PerSegmentPerStrand = PerSegment<PerStrand<T>>;

template <typename T>
using PerSegmentPerStrandPerBase = PerSegment<PerStrand<PerBase<T>>>;

} // namespace reseq

#endif // RESEQ_CONTAINER_TYPES_H
```

### Application Scope

Apply to **header declarations only** — the 35-40 deeply nested type declarations in Stats class headers. Examples:

**Before (ErrorStats.h):**
```cpp
std::array<std::array<std::array<std::vector<std::vector<std::vector<
    utilities::VectorAtomic<uintNucCount>>>>, 5>, 4>, 2>
    tmp_called_bases_by_base_quality_per_tile_;
```

**After:**
```cpp
PerSegmentPerBase<PerBaseN<std::vector<std::vector<AtomicVec<uintNucCount>>>>>
    tmp_called_bases_by_base_quality_per_tile_;
```

**Not changed:**
- `.cpp` file local variables and loop structures — too much churn for marginal readability gain
- Containers where the array dimension doesn't correspond to a domain concept (e.g., arbitrary size arrays)
- Serialization code — field order must remain identical

### Verification

`make build && make test` — type aliases are transparent to the compiler.

---

## Step 3: FRIEND_TEST Removal & gtest Cleanup (6e)

### Goal

Remove all 6 `FRIEND_TEST` macros from production headers. Remove `gtest/gtest.h` include from `Vect.hpp`. The existing `friend class *Test;` declarations stay — they already provide the private access tests need. Only introduce Peer classes if a specific test still needs narrower private access after removing FRIEND_TEST.

### Strategy Per Macro

| Location | Macro | Action |
|----------|-------|--------|
| `Vect.hpp` | `FRIEND_TEST(VectTest, BasicFunctionality)` | **Eliminate.** Rewrite test to use public API (`.std()`, iterators, `size()`, `from()`, `to()`). Remove `#include "gtest/gtest.h"` from Vect.hpp. |
| `Vect.hpp` | `FRIEND_TEST(VectTest, CopyAndClear)` | **Eliminate.** Same — public API is sufficient. |
| `DataStats.h` | `FRIEND_TEST(DataStatsTest, Construction)` | **Eliminate.** `friend class DataStatsTest;` already exists (line 136). The FRIEND_TEST is redundant — the Construction test body runs inside `DataStatsTest` fixture, which is already a friend. Rewrite the test to use public getters (`ReadLengths(0).size() == 0`, `SequenceContent(0, 0).size() == 0`) so the friend is not exercised. |
| `ReadSequenceStats.h` | `FRIEND_TEST(DataStatsTest, Construction)` | **Eliminate.** Same — `friend class DataStatsTest;` already exists (line 57). Rewrite test to use public API. |
| `FragmentDistributionStats.h` | `FRIEND_TEST(FragmentDistributionStatsTest, UpdateRefSeqBias)` | **Eliminate.** `friend class FragmentDistributionStatsTest;` already exists (line 517). The FRIEND_TEST is redundant. No test rewrite needed — existing friend class covers access. |
| `FragmentDuplicationStats.h` | `FRIEND_TEST(FragmentDuplicationStatsTest, DispersionCalculation)` | **Eliminate.** `friend class FragmentDuplicationStatsTest;` already exists (line 40). Redundant. |
| `FragmentDuplicationStats.h` | `FRIEND_TEST(FragmentDistributionStatsTest, BiasBinningAndFragmentCounts)` | **Eliminate.** Cross-class test access — `friend class FragmentDistributionStatsTest;` does NOT exist in FragmentDuplicationStats.h. This FRIEND_TEST grants access to FragmentDistributionStatsTest for one specific test. Replace with `friend class FragmentDuplicationStatsPeer;` only if the test genuinely needs private access after analysis; otherwise refactor the test to access via the existing FragmentDuplicationStatsTest friend. |

### gtest Cleanup

After removing all `FRIEND_TEST` from `Vect.hpp`, remove the `#include "gtest/gtest.h"` line. This eliminates the transitive gtest dependency from all 8 Stats class headers that include Vect.hpp.

Note: Other production headers that use `FRIEND_TEST` (DataStats.h, ReadSequenceStats.h, FragmentDistributionStats.h, FragmentDuplicationStats.h) get gtest transitively through `utilities.hpp` → `Vect.hpp`. Once Vect.hpp drops the include, verify these headers still compile (they should — `friend class` does not require gtest).

### Verification

`make build && make test` — all tests must pass. No test coverage loss. The Vect tests (BasicFunctionality, CopyAndClear) must be rewritten to use public API before removing FRIEND_TEST.

---

## Step 4: Prune DataStatsInterface (6d)

### Goal

Remove unused methods from DataStatsInterface. The only consumer of DataStatsInterface is the Python/SWIG binding (`plotDataStats.py` via `python/DataStats.i`). The `queryProfile` CLI command uses `DataStats` directly — it is NOT a DataStatsInterface consumer.

### Design

#### 4a: Audit Python Usage

Grep `plotDataStats.py` and `python/DataStats.i` for all DataStatsInterface method calls. Remove every method from DataStatsInterface.h/.cpp that is not called by either file.

#### 4b: Keep `.std()` Wrapping (No Binding Contract Change)

The current `.std()` wrapping converts `Vect<T>` to `pair<size_type, vector<T>>&`, which matches the SWIG `%template` declarations in `python/DataStats.i`. Changing this return type would require redesigning the SWIG bindings — a separate concern with its own risk profile.

**Decision:** Keep the `.std()` wrapping pattern for this phase. A SWIG binding redesign (removing `.std()`, exposing `Vect<T>` to Python) can be a follow-up phase if desired.

#### 4c: Simplify Where Possible

For remaining methods, look for opportunities to consolidate or simplify forwarding chains, but do not change return types or the SWIG contract.

### Expected Outcome

- DataStatsInterface: reduced to only Python-consumed methods, ~200-400 LOC (down from 679)
- `.std()` wrapping unchanged — SWIG bindings remain stable
- `queryProfile` CLI unaffected (it does not use DataStatsInterface)

### Verification

- `make build && make test`
- Run `plotDataStats.py` on a real `.reseq` profile to verify Python output is correct

---

## Files Created

| File | Step | Purpose |
|------|------|---------|
| `reseq/constants.hpp` | 6b | Named domain constants |
| `reseq/container_types.hpp` | 6c | Template aliases for nested containers |

## Files Modified (Summary)

| Step | Files | Nature of Change |
|------|-------|------------------|
| 6b | ~15 Stats/utility files | Replace magic numbers with constants |
| 6c | ~8 Stats headers | Replace nested type declarations with aliases |
| 6e | 5 production headers + 2-3 test files | Remove FRIEND_TEST macros, remove gtest include from Vect.hpp, rewrite Vect/DataStats tests to use public API |
| 6d | DataStatsInterface.h/.cpp | Audit Python usage, remove uncalled methods |

## Risk Mitigations

| Risk | Mitigation |
|------|------------|
| Wrong magic number replaced | Each replacement verified in context; `4` in non-base contexts left as-is |
| Container alias breaks serialization | Aliases are transparent to compiler; serialization order unchanged |
| FRIEND_TEST removal loses test coverage | Rewrite Vect/DataStats tests to use public API; other tests keep existing `friend class` access |
| DataStatsInterface pruning breaks Python | Test plotDataStats.py against real profile after changes; SWIG contract (.std() wrapping) unchanged |
| Type alias confuses IDE/tooling | Standard C++ using declarations; all major IDEs resolve them |
| gtest include removal breaks compilation | Verify all production headers compile without transitive gtest from Vect.hpp |
