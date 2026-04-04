# Phase 6: DRY & Interface Cleanup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace magic numbers with named constants, deeply nested container types with readable aliases, remove FRIEND_TEST from production headers, and prune unused DataStatsInterface methods.

**Architecture:** Four sequential sub-tasks (6b, 6c, 6e, 6d), each independently mergeable. No logic changes — purely mechanical cleanup. All changes are transparent to the compiler and preserve binary serialization compatibility.

**Tech Stack:** C++20, GoogleTest, Boost.Serialization, SWIG/Python

**Design spec:** `docs/superpowers/specs/2026-04-04-phase6-dry-interface-cleanup-design.md`

**Prerequisite:** Phase 5b merged. `make build && make test` passes on master.

---

## Critical Invariants

1. **All existing tests must pass.** `make build && make test` after every task.
2. **Binary serialization compatibility preserved.** Type aliases are transparent — no data layout changes.
3. **Conventional commits** with `(6)` scope tag and sub-task letter.
4. **No public DataStats API changes.**

---

## Task 1: Create Named Constants (6b)

**Goal:** Create `reseq/constants.hpp` with domain constants that replace magic numbers.

**Files created:**
- `reseq/constants.hpp`

### Step-by-step

- [ ] **Step 1.1: Create `reseq/constants.hpp`**

```cpp
#ifndef RESEQ_CONSTANTS_H
#define RESEQ_CONSTANTS_H

#include "utilities.hpp"

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

- [ ] **Step 1.2: Build and test**

```bash
make build && make test
```

Header-only file, no CMake changes needed. Just verify it compiles when included.

- [ ] **Step 1.3: Format and commit**

```bash
make format
git add reseq/constants.hpp
git commit -m "refactor(6b): add named domain constants header"
```

---

## Task 2: Apply Named Constants to Codebase (6b)

**Goal:** Replace magic numbers with constants from `constants.hpp` across all production source files. This is a large mechanical replacement — do it file by file, building after each file to catch errors immediately.

**Files modified:** All `.h` and `.cpp` files in `reseq/` that use the magic numbers (approximately 15-20 files).

**Replacement rules:**
- Replace `4` with `kNumBases` in array template params (`std::array<..., 4>`) and loop bounds (`for (auto base = 4; base--;)`) where `4` represents DNA bases
- Replace `5` with `kNumBasesN` where it represents ACGTN
- Replace `2` with `kTemplateSegments` where it represents first/second read (identifiable by variable name `template_segment`)
- Replace `2` with `kStrands` where it represents forward/reverse strand (identifiable by variable name `strand`)
- Replace `33` with `kPhredSangerOffset` in quality offset code
- Replace `64` with `kPhredIlluminaOffset` in quality offset code
- Replace `101` with `kGCBins` in GC content bin array declarations and loop bounds

**IMPORTANT disambiguation:** Not every `4`, `5`, or `2` is a domain constant. Only replace when the number provably represents the domain concept. Verify each replacement in context. Examples of FALSE positives to leave alone:
- `Surrounding::Length()` uses `4` but it means surrounding context length, not bases
- `cigar_element.count` comparisons
- Array indices like `.at(0)`, `.at(1)` 
- Constants in `ErrorStats.h` line 53 where `6` represents `kNumBasesN + 1` (previous called base including "none")

### Step-by-step

- [ ] **Step 2.1: Replace in header files**

Add `#include "constants.hpp"` to each header that gains a constant reference. Apply replacements in these headers:
- `ErrorStats.h` — replace `4` (kNumBases), `5` (kNumBasesN) in array template params
- `QualityStats.h` — replace `4`, `5` in array template params
- `CoverageStats.h` — replace `4`, `5` in array template params and inline loop bounds
- `DataStats.h` — replace `4` (kNumBases) in `tmp_sequence_content_reference_` and `sequence_content_reference_`
- `ReadSequenceStats.h` — replace `5` (kNumBasesN) in array template params
- `AdapterStats.h` — replace `5` (kNumBasesN) in array template params
- `FragmentDistributionStats.h` — replace `101` (kGCBins) in array template params
- `utilities.hpp` — replace `5` in `SeqQualityStats` content array
- `Simulator.h` — replace `2` (kTemplateSegments) in loop bounds

Build after all headers are done: `make build`

- [ ] **Step 2.2: Replace in `.cpp` implementation files**

Apply replacements in these files:
- `ErrorStats.cpp` — replace `4` (kNumBases), `5` (kNumBasesN) in loop bounds
- `QualityStats.cpp` — replace `4`, `5`, `2` (kTemplateSegments) in loop bounds
- `CoverageStats.cpp` — replace `4`, `5`, `2` in loop bounds
- `BamIngestionEngine.cpp` — replace `4` (kNumBases), `33` (kPhredSangerOffset), `64` (kPhredIlluminaOffset), `2` (kTemplateSegments/kStrands)
- `DataStats.cpp` — replace `2` (kTemplateSegments) in constructor loop
- `ReadSequenceStats.cpp` — replace `5` (kNumBasesN), `2` (kTemplateSegments) in loops
- `Simulator.cpp` — replace `2` (kTemplateSegments) in loops
- `AdapterStats.cpp` — replace `2` (kTemplateSegments) in loops
- `FragmentDistributionStats.cpp` — replace `4` (kNumBases), `101` (kGCBins) in loops and array declarations
- `ProbabilityEstimates.cpp` — replace `2` (kTemplateSegments) in loops
- `Reference.cpp` — replace `4` (kNumBases) in loops

Build after all files are done: `make build`

- [ ] **Step 2.3: Build and test**

```bash
make build && make test
```

All 95 tests must pass.

- [ ] **Step 2.4: Format and commit**

```bash
make format
git add -A reseq/
git commit -m "refactor(6b): replace magic numbers with named constants

Replace magic 4/5/2/33/64/101 with kNumBases, kNumBasesN,
kTemplateSegments, kStrands, kPhredSangerOffset, kPhredIlluminaOffset,
kGCBins across ~20 source files. Each replacement verified in context."
```

---

## Task 3: Create Container Type Aliases (6c)

**Goal:** Create `reseq/container_types.hpp` with template aliases for deeply nested container declarations.

**Files created:**
- `reseq/container_types.hpp`

### Step-by-step

- [ ] **Step 3.1: Create `reseq/container_types.hpp`**

```cpp
#ifndef RESEQ_CONTAINER_TYPES_H
#define RESEQ_CONTAINER_TYPES_H

#include <array>
#include <vector>

#include "constants.hpp"
#include "utilities.hpp"

namespace reseq {

// Atomic accumulator vector (thread-safe, filled during multi-threaded BAM pass)
template <typename T>
using AtomicVec = std::vector<utilities::VectorAtomic<T>>;

// Domain-indexed array aliases
template <typename T>
using PerSegment = std::array<T, kTemplateSegments>;

template <typename T>
using PerStrand = std::array<T, kStrands>;

template <typename T>
using PerBase = std::array<T, kNumBases>;

template <typename T>
using PerBaseN = std::array<T, kNumBasesN>;

// Composite aliases
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

- [ ] **Step 3.2: Build and test**

```bash
make build && make test
```

- [ ] **Step 3.3: Commit**

```bash
make format
git add reseq/container_types.hpp
git commit -m "refactor(6c): add container type alias header"
```

---

## Task 4: Apply Container Type Aliases to Headers (6c)

**Goal:** Replace deeply nested `std::array` declarations in Stats class headers with the new type aliases. Apply to header declarations only — not `.cpp` files.

**Files modified:**
- `reseq/ErrorStats.h`
- `reseq/QualityStats.h`
- `reseq/CoverageStats.h`
- `reseq/DataStats.h`
- `reseq/ReadSequenceStats.h`
- `reseq/AdapterStats.h`
- `reseq/FragmentDistributionStats.h`

### Step-by-step

- [ ] **Step 4.1: Replace in ErrorStats.h**

Add `#include "container_types.hpp"` at top. Replace the 7 identical `tmp_*_per_tile_` declarations (lines 20-40) that follow the pattern:
```cpp
std::array<std::array<std::array<std::vector<std::vector<std::vector<
    utilities::VectorAtomic<uintNucCount>>>>, kNumBasesN>, kNumBases>, kTemplateSegments>
```
with:
```cpp
PerSegmentPerBase<PerBaseN<std::vector<std::vector<AtomicVec<uintNucCount>>>>>
```

Also replace the corresponding final histogram declarations (lines 57-77) and the `called_bases_by_base_quality_per_previous_called_base_` declarations.

Build: `make build`

- [ ] **Step 4.2: Replace in QualityStats.h**

Add `#include "container_types.hpp"`. Replace the ~26 deeply nested declarations. The dominant patterns are:
- `PerSegmentPerBase<PerBaseN<...>>` for per-tile-per-error-reference vars
- `PerSegment<PerBase<...>>` for per-tile quality vars
- `PerSegmentPerBaseN<...>` for per-sequence vars

Build: `make build`

- [ ] **Step 4.3: Replace in CoverageStats.h**

Add `#include "container_types.hpp"`. Replace the ~12 deeply nested declarations. The dominant pattern here involves `PerBase<PerBaseN<PerBase<...>>>` for error rate tables.

Build: `make build`

- [ ] **Step 4.4: Replace in DataStats.h, ReadSequenceStats.h, AdapterStats.h, FragmentDistributionStats.h**

Apply aliases to the remaining files (2-3 declarations each). These are the simplest replacements.

Build: `make build`

- [ ] **Step 4.5: Build and test**

```bash
make build && make test
```

All 95 tests must pass. Type aliases are transparent to the compiler — if it builds, it's correct.

- [ ] **Step 4.6: Format and commit**

```bash
make format
git add reseq/ErrorStats.h reseq/QualityStats.h reseq/CoverageStats.h reseq/DataStats.h reseq/ReadSequenceStats.h reseq/AdapterStats.h reseq/FragmentDistributionStats.h
git commit -m "refactor(6c): apply container type aliases to Stats headers

Replace 35+ deeply nested std::array/vector declarations with
PerSegment, PerBase, PerBaseN, AtomicVec and composite aliases.
Header declarations only — no .cpp changes."
```

---

## Task 5: Remove FRIEND_TEST and gtest from Production Headers (6e)

**Goal:** Remove all 6 `FRIEND_TEST` macros from production headers. Remove `#include "gtest/gtest.h"` from `Vect.hpp`.

**Files modified:**
- `reseq/Vect.hpp`
- `reseq/VectTest.cpp`
- `reseq/DataStats.h`
- `reseq/ReadSequenceStats.h`
- `reseq/DataStatsTest.cpp`
- `reseq/FragmentDistributionStats.h`
- `reseq/FragmentDuplicationStats.h`

### Step-by-step

- [ ] **Step 5.1: Rewrite VectTest to use public API**

In `reseq/VectTest.cpp`, rewrite `BasicFunctionality` and `CopyAndClear` tests to use public methods instead of private `vec_` member:

For `BasicFunctionality`: Replace `test.vec_.first` with `test.from()`, `test.vec_.second.empty()` with `test.empty()`, `test.vec_.second.size()` with `test.size()`.

For `CopyAndClear`: The test sets `vec_.first` and `vec_.second.push_back()` directly. Replace with a Vect constructor that takes offset + initializer, or use the `Acquire()` method to set up test state through public API.

Build and test: `make build && make test`

- [ ] **Step 5.2: Remove FRIEND_TEST and gtest from Vect.hpp**

Remove these lines from `reseq/Vect.hpp`:
```cpp
#include "gtest/gtest.h"
```
and:
```cpp
    FRIEND_TEST(VectTest, BasicFunctionality);
    FRIEND_TEST(VectTest, CopyAndClear);
```

Build: `make build` — verify no compilation errors from missing gtest in transitive includes.

- [ ] **Step 5.3: Rewrite DataStatsTest Construction test**

In `reseq/DataStatsTest.cpp`, rewrite the `Construction` test (lines 487-506) to use public getters:

```cpp
// Old: EXPECT_TRUE(test_->read_lengths_.at(0).empty())
// New: EXPECT_EQ(0, test_->ReadLengths(0).size())

// Old: EXPECT_TRUE(test_->read_sequence_stats_.sequence_content_.at(i).at(j).empty())
// New: EXPECT_EQ(0, test_->SequenceContent(i, j).size())
```

The `Shrink()` call stays — it's a public method on DataStats (it was kept there in Phase 5b).

Build and test: `make build && make test`

- [ ] **Step 5.4: Remove FRIEND_TEST from DataStats.h and ReadSequenceStats.h**

Remove from `reseq/DataStats.h`:
```cpp
    FRIEND_TEST(DataStatsTest, Construction);
```

Remove from `reseq/ReadSequenceStats.h`:
```cpp
    FRIEND_TEST(DataStatsTest, Construction);
```

Build: `make build`

- [ ] **Step 5.5: Remove FRIEND_TEST from FragmentDistributionStats.h and FragmentDuplicationStats.h**

Remove from `reseq/FragmentDistributionStats.h`:
```cpp
    FRIEND_TEST(FragmentDistributionStatsTest, UpdateRefSeqBias);
```
The existing `friend class FragmentDistributionStatsTest;` already covers this test's access.

Remove from `reseq/FragmentDuplicationStats.h`:
```cpp
    FRIEND_TEST(FragmentDuplicationStatsTest, DispersionCalculation);
```
The existing `friend class FragmentDuplicationStatsTest;` already covers this.

For the cross-class FRIEND_TEST:
```cpp
    FRIEND_TEST(FragmentDistributionStatsTest, BiasBinningAndFragmentCounts);
```
Investigate whether this test actually accesses FragmentDuplicationStats private members. If yes, add `friend class FragmentDistributionStatsTest;` to FragmentDuplicationStats.h (replacing the FRIEND_TEST). If no, just remove the FRIEND_TEST.

Build: `make build`

- [ ] **Step 5.6: Build and test**

```bash
make build && make test
```

All 95 tests must pass. Verify no gtest include remains in any production header:
```bash
grep -rn "gtest/gtest.h" reseq/*.hpp reseq/*.h | grep -v Test
```
Should return empty.

- [ ] **Step 5.7: Format and commit**

```bash
make format
git add reseq/Vect.hpp reseq/VectTest.cpp reseq/DataStats.h reseq/ReadSequenceStats.h reseq/DataStatsTest.cpp reseq/FragmentDistributionStats.h reseq/FragmentDuplicationStats.h
git commit -m "refactor(6e): remove FRIEND_TEST and gtest from production headers

Remove all 6 FRIEND_TEST macros from production headers. Remove
gtest/gtest.h include from Vect.hpp, eliminating transitive gtest
dependency from all Stats headers. Rewrite Vect and DataStats
construction tests to use public API."
```

---

## Task 6: Prune DataStatsInterface (6d)

**Goal:** Remove unused methods from DataStatsInterface. Keep `.std()` wrapping and SWIG contract unchanged.

**Files modified:**
- `reseq/DataStatsInterface.h`
- `reseq/DataStatsInterface.cpp`

### Step-by-step

- [ ] **Step 6.1: Audit Python consumer usage**

Run these commands to find which methods are actually called:

```bash
# Methods called in plotDataStats.py
grep -oP '\.\K[A-Z][a-zA-Z]+(?=\()' python/plotDataStats.py | sort -u

# Methods declared in DataStatsInterface.h
grep -oP '^\s+\S+\s+\K[A-Z][a-zA-Z]+(?=\()' reseq/DataStatsInterface.h | sort -u

# Diff to find unused
comm -23 <(grep ... DataStatsInterface.h) <(grep ... plotDataStats.py)
```

Based on the audit, the following methods are NOT called by `plotDataStats.py`:
- `AdapterOverrunBases`
- `AdapterPolyATailLength`
- `BaseQualityStats` (raw stats, not the derived min/max/median/quartile)
- `BaseQualityStatsReference` (raw stats)
- `PhredQualityOffset`
- `ReadBam`
- `Save`
- `TotalNumberReads`

Verify each candidate is truly unused by also checking `python/DataStats.i` for any direct references.

- [ ] **Step 6.2: Remove unused methods**

Remove each confirmed-unused method from both `DataStatsInterface.h` (declaration) and `DataStatsInterface.cpp` (implementation). Remove methods one at a time or in small groups, building between each group.

- [ ] **Step 6.3: Build and test**

```bash
make build && make test
```

- [ ] **Step 6.4: Verify Python bindings still work**

If Python bindings are buildable (`-DRESEQ_BUILD_PYTHON=ON`), build and verify. If not, at minimum verify that the SWIG interface file `python/DataStats.i` still references valid methods.

```bash
grep -oP '%include\s+"[^"]+"|%template\s+' python/DataStats.i
```

- [ ] **Step 6.5: Format and commit**

```bash
make format
git add reseq/DataStatsInterface.h reseq/DataStatsInterface.cpp
git commit -m "refactor(6d): prune unused DataStatsInterface methods

Remove ~10 methods not called by plotDataStats.py (the sole consumer).
SWIG binding contract (.std() wrapping) unchanged."
```

---

## Verification Checklist

After all 6 tasks are complete, verify:

- [ ] `make build` succeeds with no warnings in changed files
- [ ] `make test` passes all 95 tests (0 skipped)
- [ ] `make format-check` is clean
- [ ] `pre-commit run --all-files` passes
- [ ] No `FRIEND_TEST` in production headers: `grep -rn FRIEND_TEST reseq/*.h reseq/*.hpp | grep -v Test`
- [ ] No `gtest/gtest.h` in production headers: `grep -rn "gtest/gtest.h" reseq/*.h reseq/*.hpp | grep -v Test`
- [ ] Magic numbers reduced: `grep -rn "std::array<.*,\s*[45]\s*>" reseq/*.h | wc -l` should be ~0
- [ ] Binary serialization round-trip (DataStatsTest::Ecoli save/load/verify passes)
