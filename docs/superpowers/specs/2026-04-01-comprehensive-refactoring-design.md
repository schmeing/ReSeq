# ReSeq Comprehensive Refactoring — Design Specification

**Date:** 2026-04-01
**Status:** Approved (rev 2 — incorporates code review findings)
**Scope:** Full modernization of the ReSeq codebase (C++20, modern CMake, safety, concurrency, decomposition, testing, polish)
**Approach:** 8 sequential phases, each merged to master before the next begins
**Input:** [CODE_REVIEW_REPORT.md](../../../.planning/CODE_REVIEW_REPORT.md)

---

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| C++ standard | C++20 | Unlocks `std::jthread`, concepts, `<format>`, ranges |
| CMake minimum | 3.16 | `FetchContent`, `gtest_discover_tests`, modern target-based |
| Compiler floor | GCC 10+, Clang 12+ | C++20 support |
| Vendored deps | Hybrid — FetchContent for GoogleTest/NLopt, keep SeqAn2/Skewer vendored | SeqAn3 is a full API rewrite; Skewer has local mods |
| Test architecture | Separate `reseq_test` binary + CTest, drop `reseq test` subcommand | Tests must not ship in production binary |
| Regression testing | Golden file tests with fixed seed + single thread | Safety net for all subsequent phases |
| Phasing | Strict sequential (A merges before B starts) | Clean bisect, clean master at every step |
| Test coverage | Fill gaps + E2E + coverage instrumentation (gcov/lcov) | Visibility, no hard threshold |
| ROOTPWA | Absorb into `reseq/logging.hpp` + `reseq/format_utils.hpp`, delete vendored dir | Single 404-line file, used in 17 files, no reason to keep vendored |
| Multi-thread determinism | Single-thread golden files only; multi-thread tests use bounded invariants | Cross-thread-count byte-identical output is not achievable without algorithmic redesign (see Phase 7 notes) |

---

## Phase 1: Build System & Test Architecture

**Goal:** Modern CMake with C++20, separate test binary, CTest, FetchContent.

### CMake Migration Strategy

The current build uses global-state CMake (`include_directories()`, `link_directories()`, `link_libraries()`, `CMAKE_CXX_FLAGS` string manipulation). Converting to target-based CMake must be **incremental** — never remove a global command until all consumers have explicit `target_link_libraries()` to a wrapper target.

**Step order (each step builds and passes tests):**

1. **Bump minimum to 3.16**, no functional changes — unlocks all required features
2. **Set `CMAKE_CXX_STANDARD 20`** with `CMAKE_CXX_STANDARD_REQUIRED ON`
3. **Wrap vendored deps as INTERFACE targets** (one at a time, least entangled first):
   - GoogleTest → `FetchContent` (remove vendored `googletest/`)
   - ROOTPWA → `vendored::rootpwa` INTERFACE target (temporary — absorbed in Phase 3)
   - NLopt → `FetchContent` with `find_package()` fallback
   - Skewer → `vendored::skewer` INTERFACE target with `SYSTEM` includes
   - SeqAn2 → `vendored::seqan` INTERFACE target with `SYSTEM` includes (last — most pervasive)
4. **Convert system deps** to imported targets: `Boost::serialization`, `Boost::program_options`, `Boost::filesystem`, `Boost::system`, `ZLIB::ZLIB`, `BZip2::BZip2`
5. **Extract `reseq_lib` static library** from production sources (all `reseq/*.cpp` except `main.cpp` and `*Test.*`)
6. **Convert compiler flags** from `CMAKE_CXX_FLAGS` to `target_compile_options()` / `target_compile_features()`
7. **Remove all global state** (`include_directories()`, `link_directories()`, `link_libraries()`) — only after every target has explicit dependencies

### CMake Targets (final state)

- `reseq_lib` — static library, all core source (`reseq/` except `main.cpp` and `*Test.*`)
- `reseq` — thin CLI executable linking `reseq_lib`
- `reseq_test` — test executable linking `reseq_lib` + GoogleTest

### CMake Options

| Option | Default | Purpose |
|--------|---------|---------|
| `RESEQ_BUILD_TESTS` | ON | Build test binary |
| `RESEQ_BUILD_PYTHON` | OFF | Build Python/SWIG bindings |
| `CODE_COVERAGE` | OFF | Enable gcov instrumentation |

### Test Harness Migration

The current test harness is deeply coupled to the production binary:
- `reseq/CMakeLists.txt:18` compiles all 15 `*Test.cpp` files into the `reseq` executable
- `reseq/main.cpp:1179-1226` contains test registration, sequencing, and execution
- Tests access private state via `FRIEND_TEST` (6 macros) and direct member access
- Test data discovered via `utilities::GetReSeqDir()` which uses `PROJECT_SOURCE_DIR` from `CMakeConfig.h`

**Migration path:**

1. **Capture pre-split baseline checksums** — before any build changes, run `reseq test` and record pass/fail state + any output checksums. This is a manual verification step, not automated golden files (those come in Phase 2).
2. **Create `reseq_test` executable** — link `reseq_lib` + GoogleTest + all `*Test.cpp` files
3. **Replicate test registration** — extract the sequential registration logic from `main.cpp:1196-1206` into a test main (or use `gtest_discover_tests` for automatic registration)
4. **Propagate `PROJECT_SOURCE_DIR`** to the test binary via `CMakeConfig.h` (same mechanism as today) so `GetReSeqDir()` finds `test/` data
5. **Verify identical test behavior** — run both `reseq test` (old) and `ctest` (new), confirm identical pass/fail results
6. **Remove test code from production binary** — delete `*Test.cpp` from `reseq` target, remove GoogleTest link, remove `reseq test` subcommand from `main.cpp`
7. **Handle `FRIEND_TEST` temporarily** — keep the 6 `FRIEND_TEST` macros in production headers for now (they compile harmlessly without GoogleTest). Full removal happens in Phase 6e after decomposition creates cleaner public interfaces.

### CTest Integration

- `enable_testing()` + `gtest_discover_tests(reseq_test)`
- Tests run via `ctest -j$(nproc)` or `make test`
- Per-test temp directories via `WORKING_DIRECTORY` property

### CI Updates

The repo already has GitHub Actions CI (`.github/workflows/ci.yml`) with GCC 13 and Clang 17 on Ubuntu 24.04. This phase **extends** existing CI, not creates it from scratch:

- Update compiler floor requirements in CI documentation (GCC 10+, Clang 12+)
- Switch test execution from `build/bin/reseq test` to `ctest --output-on-failure`
- Add sanitizer job: `-fsanitize=address,undefined` as separate CI matrix entry
- Add coverage job: `CODE_COVERAGE=ON` + gcov/lcov, upload as artifact
- Keep existing format-check and version-check jobs

---

## Phase 2: Golden File Regression Tests

**Goal:** Capture current output as baseline. These become the safety net for all subsequent phases.

### Strategy

- Single-threaded, fixed-seed execution (eliminates scheduling variance)
- Small test data from `test/` (E. coli reference + BAM files)
- Simple file-diff helper in GoogleTest (no external framework)
- Tests run as part of the new `reseq_test` binary (built in Phase 1), invoking the `reseq` CLI binary as a subprocess

### Regression Tests

| Test | Command | Golden File | Comparison |
|------|---------|-------------|------------|
| `illuminaPE` stats | Run on test BAM → `.reseq` profile | `test/expected/ecoli.reseq` | Binary compare |
| `illuminaPE` simulation | Profile → FASTQ (seed=42, threads=1) | `test/expected/ecoli_sim_R{1,2}.fq.gz` | Decompress + diff |
| `seqToIllumina` | Input sequences → FASTQ (seed=42, threads=1) | `test/expected/seq_to_illumina_R{1,2}.fq.gz` | Decompress + diff |
| `queryProfile` | Extract from `.reseq` file | `test/expected/query_profile.txt` | Text diff |
| `replaceN` | Replace Ns in reference | `test/expected/replaceN.fa` | Text diff |

### Test Architecture

These regression tests are **subprocess-based** — they invoke the `reseq` CLI binary and compare its file output. This is distinct from the unit tests which link against `reseq_lib` and test internal APIs. The subprocess pattern:

1. Locate the `reseq` binary via CMake-configured path (`RESEQ_BINARY_DIR`)
2. Run it with fixed arguments (seed, threads=1, test data paths)
3. Compare output files against committed golden files
4. Each test uses a per-test temp directory (cleaned up on success, preserved on failure for debugging)

```cpp
void ExpectFilesEqual(const std::filesystem::path& expected,
                      const std::filesystem::path& actual);
void ExpectGzFilesEqual(const std::filesystem::path& expected,
                        const std::filesystem::path& actual);
```

### Workflow

1. Build with Phase 1 changes (new test binary + production binary)
2. Run each `reseq` command, capture output as golden files
3. Commit golden files to `test/expected/`
4. Write GoogleTest fixtures in `reseq_test` that invoke `reseq` binary and diff against golden files
5. Verify all pass before proceeding

### CI Gate

Every subsequent phase PR must pass regression tests. If a phase intentionally changes output, golden files are updated with justification in the commit message.

### Key Constraint: Cross-Platform Determinism

`std::uniform_real_distribution` output is not portable across stdlib implementations (libstdc++ vs libc++). Solution: **golden files are committed for GCC/libstdc++ only** (the primary CI platform). The Clang CI job runs regression tests in "generate-and-self-check" mode — it generates its own golden files and verifies they don't change across runs, but does not compare against the GCC golden files. This avoids maintaining duplicate golden files while still catching regressions on both compilers.

---

## Phase 3: Safety Stabilization

**Goal:** Memory safety, null handling, exception-safe locking, ROOTPWA absorption.

### 3a: Absorb `2016-05-15_ROOTPWA/`

Single file `reportingUtils.hpp` (404 lines, used in 17 files). Split into:

| New File | Contents |
|----------|----------|
| `reseq/logging.hpp` | `NullStream`, verbosity globals, terminal color detection, VT100 escape codes, logging macros (`printErr`/`printWarn`/`printInfo`/`printSucc`/`printDebug`) |
| `reseq/format_utils.hpp` | `maxPrecision`, `maxPrecisionAlign`, `maxPrecisionDouble`, `indent`, stream operators for `std::pair`/`std::vector`, `nmbOfDigits` |

Changes during absorption:
- Move `rpwa::` namespace contents to `reseq::`
- `<stdlib.h>` → `<cstdlib>`, `<stdint.h>` → `<cstdint>`
- `enum` → `enum class` for `vt100EscapeCodesEnum`
- `sizeof(array)/sizeof(element)` → `std::size()`
- Update all 17 `#include "reportingUtils.hpp"` sites
- Delete `2016-05-15_ROOTPWA/` directory
- Remove `vendored::rootpwa` CMake target (created temporarily in Phase 1)

### 3b: `NULL` → `nullptr`

82 sites across 31 files. Mechanical replacement.

### 3c: `mutex.lock()`/`unlock()` → `std::scoped_lock`

~20 manual lock/unlock pairs in `ProbabilityEstimates.h` and `FragmentDistributionStats.cpp`. Replace with `std::scoped_lock` (handles multiple mutexes, exception-safe).

### 3d: Raw `new`/`delete` → `std::unique_ptr`

124 raw allocation sites, prioritized by risk:

| Priority | Location | Pattern |
|----------|----------|---------|
| HIGH | `DataStats` → `CoverageStats` (`FullRecord`) | `unique_ptr` with ownership transfer |
| HIGH | `Simulator` (`SimUnit`/`SimBlock` linked lists) | `unique_ptr` chains or `std::vector<unique_ptr>` |
| HIGH | `CoverageStats` (`CoverageBlock` linked list) | `std::vector<unique_ptr>` |
| MEDIUM | `Simulator` (`StringSet`) | `unique_ptr` member |
| LOW | Remaining local `new`/`delete` pairs | `make_unique` |

**Critical note for Phase 4 prerequisite:** The `unique_ptr` work on Simulator's `SimBlock`/`SimUnit` linked lists and CoverageStats' `CoverageBlock` chains is not just a safety improvement — it is a **required prerequisite** for Phase 4's concurrency modernization. The current intrusive linked lists with mixed atomic/raw pointers cannot be safely converted to `std::jthread` without first normalizing ownership. Specifically:
- `SimBlock::next_block_` is `std::atomic<SimBlock*>` but `partner_block_`, `first_block_`, `last_block_` are raw — mixed atomicity
- `CoverageStats::reusable_blocks_` is accessed by worker threads without mutex protection
- `Simulator::GetNextBlock()` traverses `first_unit_->first_block_->next_block_` without locks while teardown deletes the same nodes

Phase 3d must convert these to `std::vector<std::unique_ptr<>>` or similar owning containers, eliminating the intrusive linked-list patterns entirely. This is the container/ownership normalization that makes Phase 4's thread-primitive swap safe.

### 3e: Minor Modernizations

- 3 plain `enum` → `enum class` (ErrorStats.h, FragmentDistributionStats.h, ProbabilityEstimates.h)
- `CoverageStats.h:63` throwing copy constructor → `= delete`
- `Simulator.cpp:1295` — make failure state explicit (don't continue through failed block creation)

### Verification Gate

All 57 existing tests + all regression golden file tests pass after each sub-step.

---

## Phase 4: Concurrency Modernization

**Goal:** Replace hand-rolled schedulers with C++20 standard primitives.

**Prerequisite:** Phase 3d's ownership normalization must be complete. Phase 4 is **not** a thread-primitive swap — it is a data-structure redesign followed by thread-primitive adoption. The investigation revealed that the current architecture mixes intrusive linked structures, atomic pointer choreography, ownership transfer, and lifecycle teardown in ways that make simple `std::thread` → `std::jthread` replacement unsafe.

### 4a: Data Structure Redesign (before thread primitives)

Each concurrent subsystem needs its container normalized before thread primitives can change:

| Subsystem | Current Structure | Problem | Target Structure |
|-----------|------------------|---------|-----------------|
| `Simulator` SimBlock chain | Singly-linked list via `atomic<SimBlock*> next_block_` + raw `partner_block_`, `first_block_`, `last_block_` | Mixed atomicity; `GetNextBlock()` does unlocked traversal while teardown deletes nodes | `std::vector<std::unique_ptr<SimBlock>>` indexed by block ID; work queue via `std::queue<size_t>` protected by mutex+condvar |
| `CoverageStats` CoverageBlock | Doubly-linked list via `atomic<CoverageBlock*> next_block_` + raw `previous_block_`; `reusable_blocks_` vector | Reuse pool accessed by workers without mutex; concurrent modification corrupts iterator | `std::vector<std::unique_ptr<CoverageBlock>>` with explicit pool; mutex-guarded reuse |
| `FragmentDistributionStats` bias queue | `std::deque<pair<atomic_flag, BiasCalculationVectors>>` | Multiple threads modify deque elements via `AcquireBiases()` without locking | Mutex-guarded work queue with `std::condition_variable` for signaling |
| `ProbabilityEstimates` fit state | `LogArrayCalc` nested vectors (`dim2_`) accessed during parallel fitting | No explicit mutex around shared fit state | Explicit `std::mutex` guarding `dim2_` access, or per-thread local copies |

### 4b: Thread Primitive Adoption (after container normalization)

| Location | Current | Replacement |
|----------|---------|-------------|
| `Simulator` threading | `std::thread` + manual join + manual teardown | `std::jthread` + `std::stop_token` for cancellation; RAII teardown via owning containers |
| `CoverageStats` block processing | Manual thread lifecycle | `std::jthread` + `std::condition_variable` for producer/consumer |
| `FragmentDistributionStats` bias calc | `atomic_flag` polling, spin waits | `std::jthread` + `std::stop_token` + `std::condition_variable` |
| `ProbabilityEstimates` fit threads | `std::thread` + manual join | `std::jthread`; fit coordination via `std::condition_variable` |

### 4c: Fix Global Mutable State

`kVerbosityLevel` (`main.cpp:22`): non-atomic `uint16_t` → `std::atomic<uint16_t>`. Write-once-read-many, zero performance cost.

### 4d: Centralize Teardown

With owning containers from 4a and `std::jthread` from 4b, all manual cleanup paths become automatic:
- `Simulator.cpp:2993-3015` manual traversal+delete → implicit via vector/unique_ptr destruction
- `CoverageStats` reuse pool cleanup → implicit via vector destruction
- Thread join → implicit via `std::jthread` destructor

### Design Principle

No custom synchronization primitives. Everything uses `std::mutex` + `std::scoped_lock` + `std::condition_variable` + `std::jthread`. If a pattern can't be expressed with these, it's a design smell that needs decomposition.

### Verification Gate

All tests + regression suite. Run under `-fsanitize=thread` (TSan) in CI to verify no data races remain.

---

## Phase 5: God Object Decomposition

**Goal:** Break 5 god objects into focused, single-responsibility components.

### 5a: `main.cpp` (808 lines → ~80 lines)

```
main.cpp (slim dispatcher)
├── cli/cli_parser.hpp/cpp          — Boost.program_options setup, validation
├── cli/illumina_pe_cmd.hpp/cpp     — illuminaPE workflow
├── cli/seq_to_illumina_cmd.hpp/cpp — seqToIllumina workflow
├── cli/query_profile_cmd.hpp/cpp   — queryProfile workflow
└── cli/replace_n_cmd.hpp/cpp       — replaceN workflow
```

### 5b: `DataStats` (1,344 lines)

| New Class | Responsibility |
|-----------|---------------|
| `BamReader` | BAM iteration, record pairing, filtering |
| `RecordClassifier` | Read classification logic |
| `StatisticsAggregator` | Orchestrates sub-stats, owns lifecycle (thin shell of old `DataStats`) |

### 5c: `Simulator` (3,039 lines)

| New Class | Responsibility |
|-----------|---------------|
| `SimulationOutput` | FASTQ/BAM file writing, compression |
| `BlockManager` | SimBlock/SimUnit lifecycle, work distribution (uses normalized containers from Phase 4a) |
| `ReadGenerator` | Error injection, quality simulation, read construction |
| `Simulator` (thin) | Orchestration: configure components, run pipeline |

### 5d: `ProbabilityEstimates` (1,559-line header)

The review finding identified this as more than a mechanical file split. Investigation reveals **four concrete architectural seams**:

| New Class | Responsibility | Key Members |
|-----------|---------------|-------------|
| `IpfEngine` | IPF computation: margin definitions, iterative proportional fitting algorithm, convergence logic | `quality_`, `sequence_quality_`, `base_call_`, `dom_error_`, `error_rate_`, `indels_` (LogIPF arrays); `DefineMarginsQuality/BaseCall/...`, `IterativeProportionalFitting` |
| `FitResult` | Result materialization, post-processing, probability queries | `quality_result_`, `sequence_quality_result_`, ..., `_result_` arrays (LogArrayResult); `PrepareResult`, getter methods, `ChangeErrorRate`, `RemoveSubstitutionErrors` |
| `FitCoordinator` | Thread coordination: work-stealing loop, progress tracking, error handling | `print_mutex_`, `current_param_`, `error_during_fitting_`, `precision_improved_`; `IPFThread` dispatch |
| `FitPersistence` | Serialization, cache validation, file I/O | `stats_creation_time_`; `Load`, `Save`, boost::serialization template |

The critical seam is between `IpfEngine` and `FitResult`: `PrepareResult()` explicitly clears computation arrays (`quality_`, etc.) when materializing results, proving logical separation. Result getters should be encapsulated to avoid leaking array subscripting details.

Move from header-only to `.hpp`/`.cpp` pairs for all extracted classes.

### 5e: `FragmentDistributionStats` (3,693 lines)

| New Class | Responsibility |
|-----------|---------------|
| `AbundanceModel` | Abundance estimation, dispersion fitting |
| `InsertSizeModel` | Insert length distribution, GC bias |
| `SurroundingBiasModel` | Sequence context bias modeling |
| `FragmentPlotter` | All plotting/output logic |

### 5f: `Reference`

| New Class | Responsibility |
|-----------|---------------|
| `ReferenceGenome` | Immutable sequence data, loading, context lookup |
| `VariantSource` | VCF parsing, variant caching |
| `ExclusionRegions` | Excluded region management |

### Decomposition Strategy

Per god object:
1. Identify extraction seam (cohesive group of members + methods)
2. Extract to new class with clean public interface
3. Original class holds new class by value or `unique_ptr`
4. Update internal callers
5. Run full test + regression suite
6. Repeat for next extraction

Each extraction is one commit.

### Verification Gate

All tests + regression suite. Each new class gets at least one unit test.

---

## Phase 6: DRY & Interface Cleanup

**Goal:** Eliminate boilerplate, tighten interfaces.

### 6a: `StatVariable<T>` Registration

Replace repeated Prepare/Finalize/Shrink boilerplate with a registration pattern:

```cpp
class StatsBase {
    std::vector<StatVariableBase*> registered_vars_;
public:
    void Prepare()  { for (auto* v : registered_vars_) v->Prepare(); }
    void Finalize() { for (auto* v : registered_vars_) v->Finalize(); }
    void Shrink()   { for (auto* v : registered_vars_) v->Shrink(); }
};

template<typename T>
class StatVariable : public StatVariableBase {
    Vect<T> final_;
    Vect<T> tmp_;
};
```

~30% LOC reduction across Stats classes.

### 6b: Named Constants (`reseq/constants.hpp`)

| Constant | Value | Replaces |
|----------|-------|----------|
| `kNumBases` | 4 | Magic `4` in ~15 files |
| `kNumBasesWithN` | 5 | Magic `5` in ~10 files |
| `kPhredSangerOffset` | 33 | Magic `33` |

### 6c: Named Multi-Dimensional Container Types

Replace 5-7 level nested `std::array`/`std::vector` with typed aliases:

```cpp
using ErrorRateTable = NDArray<double, kTemplateSegments, kNumBases, kNumBases>;
```

Thin wrapper or type alias — readability goal, not a new container library.

### 6d: Narrow `DataStatsInterface`

Replace 407-line forwarding class with narrow view APIs:

| View | Purpose |
|------|---------|
| `FragmentSizeView` | Read-only fragment length distribution |
| `CoverageView` | Read-only coverage stats |
| `ErrorProfileView` | Read-only error model |

### 6e: Remove `FRIEND_TEST` from Production Headers

After Phase 5 decomposition, most of the 6 `FRIEND_TEST` macros are unnecessary — decomposed classes have cleaner public interfaces. For remaining cases:
- Prefer extracting tested logic into public methods
- If internal access is truly needed, use `#ifdef RESEQ_TESTING` guarded accessors (compile test binary with `-DRESEQ_TESTING`)
- Remove `gtest/gtest.h` include from `Vect.hpp` and any other production header

### Verification Gate

All tests + regression suite. LOC should decrease measurably.

---

## Phase 7: Test Coverage & CI Instrumentation

**Goal:** Fill test gaps, add E2E tests, add coverage visibility.

### 7a: Fill Unit Test Gaps

| Class | What to Add |
|-------|-------------|
| `ErrorStats` | Error rate calculation, base-specific profiles, edge cases |
| `QualityStats` | Distribution building, mean/median, score boundaries |
| `FragmentDuplicationStats` | Duplicate detection, duplication rate |
| New view APIs | Each narrow view interface from Phase 6 |

### 7b: Tests for Decomposed Classes (Phase 5)

At least one unit test per new class: `BamReader`, `BlockManager`, `ReadGenerator`, `SimulationOutput`, `IpfEngine`, `FitResult`, `FitCoordinator`, `AbundanceModel`, CLI subcommands, etc.

### 7c: End-to-End Pipeline Tests

| Test | Verifies |
|------|----------|
| `illuminaPE` round-trip | Stats from simulated reads → re-simulate → similar stats |
| Error handling | Invalid BAM, missing reference, corrupt `.reseq` → clean errors, non-zero exit |
| Edge cases | Empty BAM, single read, short reference, all-N reference |

### 7d: Multi-Thread Invariant Tests

**Important:** Byte-identical output across different thread counts is **not achievable** with the current simulation architecture without deep algorithmic redesign. The investigation found:

- `block_seed_gen_()` is consumed sequentially through `GetNextBlock()`, but thread scheduling changes the rate/order of consumption across thread counts
- `UpdateRefSeqBias()` consumes variable RNG state before block seeding begins
- `try_lock` paths in variant/methylation loading cause variable RNG consumption based on contention

**What we test instead:**

| Test | Assertion | Rationale |
|------|-----------|-----------|
| Single-thread determinism | seed=42, threads=1 → byte-identical to golden files | Baseline reproducibility |
| Same-config determinism | seed=42, threads=4 run twice → identical output | Same config = same result |
| Cross-thread-count invariants | threads=1 vs 2 vs 4: same read count, same mean quality, GC% within 0.1%, error rate within 0.1% | Statistical equivalence without requiring byte identity |

If byte-identical multi-thread output is desired in the future, it would require: (1) pre-generating all block seeds deterministically via `hash(initial_seed, block_id)` instead of sequential RNG consumption, (2) eliminating variable RNG consumption in `UpdateRefSeqBias()`, (3) removing `try_lock` paths. This is a separate feature, not part of this refactoring.

### 7e: Coverage Instrumentation

```cmake
option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE)
    target_compile_options(reseq_lib PRIVATE --coverage)
    target_link_options(reseq_lib PRIVATE --coverage)
endif()
```

CI: build with `CODE_COVERAGE=ON`, generate lcov report, upload as artifact or to Codecov. No hard threshold.

### 7f: Test Isolation Verification

`ctest -j$(nproc) --schedule-random` in CI — proves no inter-test dependencies.

### Verification Gate

Coverage report generated. All tests pass in random order + parallel.

---

## Phase 8: Polish & Packaging

**Goal:** Final modernization pass.

### 8a: C++ Language Cleanup

| Change | Scope |
|--------|-------|
| `<stdint.h>` → `<cstdint>` etc. | 32 files |
| `typedef` → `using` | ~30 in `utilities.hpp` + scattered |
| `static const` → `constexpr` | ~20 sites |
| `[[nodiscard]]` on key APIs | Selective |
| Structured bindings | Opportunistic |

### 8b: Split `utilities.hpp` (606 lines)

| New Header | Contents |
|------------|----------|
| `reseq/types.hpp` | Type aliases, atomic wrappers |
| `reseq/math_utils.hpp` | `SafePercent`, rounding, statistical helpers |
| `reseq/filesystem_utils.hpp` | Path handling, file checks |
| `reseq/domain_types.hpp` | Domain-specific types |

### 8c: Python Packaging

- Remove `sys.path` hacks from `plotDataStats.py`
- Proper package structure: `python/reseq/__init__.py` + `pyproject.toml` entry points
- Target py39+
- Modern CMake SWIG module (3.16+)

### 8d: Skewer `MODIFICATIONS.md`

- Diff vendored code against upstream release
- Document each modification with rationale
- Note version vendored from

### 8e: Documentation

- Refresh `README.md`: C++20, CMake 3.16+, GCC 10+/Clang 12+, `ctest`, new structure
- Brief API docs on public interfaces of decomposed classes
- Update `CLAUDE.md` to reflect new architecture

### 8f: `Vect` Container Redesign

- Document or rename offset semantics
- Add bounds checking in debug builds
- Remove `gtest/gtest.h` include (no longer needed after Phase 6 `FRIEND_TEST` removal)

### Verification Gate

Full test suite + regression tests. `make format-check` and `make lint` pass. CI green on GCC and Clang.

---

## Phase Summary

| Phase | Focus | Key Deliverable | Depends On |
|-------|-------|-----------------|------------|
| 1 | Build system & test architecture | CMake 3.16+/C++20, separate test binary, CTest, FetchContent, incremental CMake migration | — |
| 2 | Golden file regression tests | Baseline output captured via subprocess tests, safety net in place | Phase 1 |
| 3 | Safety stabilization | Smart pointers, nullptr, scoped_lock, ROOTPWA absorbed, **container/ownership normalization for concurrency** | Phase 2 |
| 4 | Concurrency modernization | Data structure redesign → std::jthread adoption, TSan-clean | Phase 3 (specifically 3d container normalization) |
| 5 | God object decomposition | 5 god objects → ~20 focused classes (ProbabilityEstimates along 4 architectural seams) | Phase 4 |
| 6 | DRY & interface cleanup | StatVariable registration, narrow APIs, FRIEND_TEST removal | Phase 5 |
| 7 | Test coverage & CI | Gap-fill tests, E2E tests, bounded multi-thread invariants, coverage instrumentation | Phase 6 |
| 8 | Polish & packaging | Language cleanup, Python fix, docs, utilities split | Phase 7 |

## Cross-Cutting Concerns

- **Every phase PR must pass:** all unit tests + golden file regression tests
- **Every phase merges to master** before the next begins
- **Commit style:** conventional commits (`feat:`, `fix:`, `refactor:`, `build:`, `test:`)
- **PRs target `berntpopp/ReSeq`**, never upstream `schmeing/ReSeq`
- **Existing CI is extended**, not replaced — current GCC 13 / Clang 17 matrix preserved

## Revision History

| Rev | Date | Changes |
|-----|------|---------|
| 1 | 2026-04-01 | Initial design |
| 2 | 2026-04-01 | Incorporated 6 code review findings: (1) explicit test harness migration path in Phase 1, (2) Phase 4 reframed as data-structure redesign + thread adoption with explicit Phase 3d prerequisite, (3) incremental CMake migration strategy with step ordering, (4) ProbabilityEstimates decomposition along 4 concrete architectural seams, (5) multi-thread determinism reframed as bounded invariants not byte-identity, (6) CI updates framed as extending existing infrastructure |
