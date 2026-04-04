# Phase 7: Test Coverage & CI Instrumentation — Design Spec

**Goal:** Establish coverage visibility, fix the CI coverage job, fill critical unit test gaps, and verify test isolation. No hard coverage threshold — use auto-ratchet (Codecov `target: auto`) to prevent regression.

**Prerequisite:** Phase 6 merged. `make build && make test` passes on master.

---

## Current State

- **141 tests** (116 unit + 25 regression) across 16 test files
- **3 core classes with no direct fixture-based tests:** ErrorStats, QualityStats, FragmentDuplicationStats — each has substantial assertion helpers (e.g., `TestSrr490124Equality`, `TestDuplicates`) exercised indirectly from DataStatsTest, but no standalone `TEST_F` targeting the class directly
- **2 Phase 5 classes with no tests:** BamIngestionEngine, ReadSequenceStats — no test files exist
- **CI coverage job broken:** `geninfo: ERROR: Unexpected negative count` from non-atomic gcov counters in threaded code
- **No Codecov integration:** coverage report uploaded as artifact only, no PR annotations or badges
- **Tests run serially** (`ctest -j1`) — no verification that tests are order-independent

---

## Sub-tasks

### 7e: Fix Coverage CI + Codecov Integration

**Priority:** First — enables visibility for all subsequent work.

#### 7e.1: Fix gcov negative count error

Add `-fprofile-update=atomic` to coverage compile options in `reseq/CMakeLists.txt`:

```cmake
if(CODE_COVERAGE)
    target_compile_options(reseq_lib PRIVATE --coverage -fprofile-update=atomic)
    target_link_options(reseq_lib PRIVATE --coverage)
    target_compile_options(reseq_test PRIVATE --coverage -fprofile-update=atomic)
    target_link_options(reseq_test PRIVATE --coverage)
endif()
```

This makes gcov counter increments atomic, eliminating the race condition that produces negative counts in multi-threaded code.

#### 7e.2: Fix lcov version mismatch

Add `--gcov-tool gcov-13` to lcov commands in `.github/workflows/ci.yml` coverage job to ensure gcov version matches the GCC 13 compiler used for the coverage build.

#### 7e.3: Integrate Codecov

Add `codecov/codecov-action@v5` step after lcov report generation in `.github/workflows/ci.yml`:

```yaml
- name: Upload coverage to Codecov
  uses: codecov/codecov-action@v5
  with:
    files: build/coverage.info
    fail_ci_if_error: false
    token: ${{ secrets.CODECOV_TOKEN }}
```

**Token policy:** Add `CODECOV_TOKEN` as a repository secret. Set `fail_ci_if_error: false` so that upload failures (e.g., fork PRs where secrets are unavailable, Codecov outages) do not block CI. Coverage is informational — build and test results are the hard gates. Fork PRs will not get coverage annotations but will still pass CI.

#### 7e.4: Configure Codecov auto-ratchet

Create `codecov.yml` at repository root:

```yaml
coverage:
  status:
    project:
      default:
        target: auto      # must not decrease from base branch
        threshold: 1%     # tolerance for measurement flakiness
    patch:
      default:
        target: 80%       # new/changed code must be 80% covered
```

This enforces:
- **Project-level:** coverage cannot decrease (ratchet)
- **Patch-level:** new code in PRs must be at least 80% covered

#### 7e.5: Add coverage badge to README

Add Codecov badge after the existing project description in `README.md`.

---

### 7f: Test Isolation Verification

**Priority:** Quick CI change, run in parallel with 7e.

**Context:** The repo explicitly documents inter-test dependencies in `reseq/CMakeLists.txt:64` — the `Register()` pattern in test classes creates shared state that requires tests to run together in a single binary invocation.

**Approach:** Add a **separate, non-blocking CI job** called `test-isolation` that runs:

```bash
build/bin/reseq_test --gtest_shuffle --gtest_random_seed=0
```

Seed `0` means "use current time" — each CI run uses a different test order. This job uses `continue-on-error: true` so it does not block PRs initially. Not using `ctest --schedule-random` because all GoogleTest cases run as a single CTest test.

**Acceptance criteria:**
- If the isolation job passes: mark as complete, consider making it blocking in a future phase.
- If the isolation job fails: file the failures as known issues. Do NOT attempt to fix `Register()` dependencies in this phase — that is a separate refactoring concern. The job stays non-blocking until dependencies are resolved.

**New tests written in this phase** (7a, 7b, 7c) MUST be order-independent by design: use GoogleTest fixtures with proper `SetUp()`/`TearDown()`, no reliance on `Register()` or global state.

---

### 7a: Add Direct Fixture-Based Tests

**Priority:** High-impact — these classes have indirect coverage via DataStatsTest helpers but no standalone `TEST_F` exercising their public API directly. Goal: add direct tests without duplicating existing golden checks.

#### 7a.1: ErrorStats tests

Add `TEST_F` tests to `reseq/ErrorStatsTest.cpp` (file exists with assertion helpers but no direct test fixtures):

| Test | What it verifies |
|------|------------------|
| `Prepare` | Allocates tmp vectors to correct dimensions for given tile/quality/position sizes |
| `AddBaseAndFinalize` | Accumulate bases via `AddBase()`, finalize, verify histograms are non-empty and consistent |
| `AddInDelAndFinalize` | Accumulate indels via `AddInDel()`, finalize, verify indel histograms |
| `Shrink` | After finalize, shrink removes trailing zeros without losing data |
| `PreparePlotting` | Derived plotting variables (called_bases_by_base_quality_, etc.) computed correctly |

Use synthetic setup: create ErrorStats, call `Prepare()` with small dimensions, inject known data via `AddBase()`/`AddInDel()`, then verify output.

#### 7a.2: QualityStats tests

Add `TEST_F` tests to `reseq/QualityStatsTest.cpp` (file exists with assertion helpers but no direct test fixtures):

| Test | What it verifies |
|------|------------------|
| `Prepare` | Allocates accumulators to correct dimensions |
| `AddQualityAndFinalize` | Accumulate quality data, finalize, verify per-tile and per-strand histograms |
| `Shrink` | Shrink preserves data integrity |
| `SummaryStatistics` | Mean, median, quartiles, min, max computed correctly from known input |

#### 7a.3: FragmentDuplicationStats tests

Add `TEST_F` tests to `reseq/FragmentDuplicationStatsTest.cpp` (file exists with assertion helpers but no direct test fixtures):

| Test | What it verifies |
|------|------------------|
| `AddDuplication` | Counts below kMaxDuplication go to correct bin, above goes to overflow bin |
| `FinalizeDuplicationVector` | Tmp vector correctly finalized into duplication_number_ Vect |
| `DuplicationNumber` | Public getter returns correct data after finalization |

These are simpler tests — FragmentDuplicationStats has a small public API.

---

### 7b: Tests for Phase 5 Decomposed Classes

New test files require test harness wiring:
1. Add `.cpp` files to the `add_executable(reseq_test ...)` source list in `reseq/CMakeLists.txt:52`
2. Add `#include` for each new `*Test.h` in `reseq/test_main.cpp:13`
3. Add `Register()` call in `reseq/test_main.cpp:44` for each new test class

#### 7b.1: BamIngestionEngine tests

Create `reseq/BamIngestionEngineTest.cpp` + `reseq/BamIngestionEngineTest.h` (new files):

| Test | What it verifies |
|------|------------------|
| `ProcessEcoliBam` | Process ecoli-SRR490124-4pairs.bam, verify record counts match expected |
| `MappingQualityFilter` | Records below minimum mapping quality are excluded |
| `ProperPairDetection` | Proper vs improper pairs classified correctly |

Uses existing ecoli BAM test data. Test fixture follows `BasicTestClassWithReference` pattern.

#### 7b.2: ReadSequenceStats tests

Create `reseq/ReadSequenceStatsTest.cpp` + `reseq/ReadSequenceStatsTest.h` (new files):

| Test | What it verifies |
|------|------------------|
| `PrepareAndFinalize` | Accumulators allocated, finalized into histograms |
| `MappingQualityAccumulation` | `IncrementMappingQuality()` correctly populates proper/improper/single histograms |
| `Shrink` | Shrink preserves data integrity |

---

### 7c: E2E Edge-Case Tests

Add to `reseq/RegressionTest.cpp` using subprocess execution (same pattern as existing regression tests — invoke `reseq` binary, check exit code and stderr):

| Test | Input | Expected behavior |
|------|-------|-------------------|
| `InvalidBamPath` | `illuminaPE` with non-existent BAM path | Non-zero exit code, stderr contains error message |
| `MissingReference` | `illuminaPE` with valid BAM, non-existent reference path | Non-zero exit code, stderr contains error message |
| `EmptyBam` | `illuminaPE` with valid BAM header but zero records | Non-zero exit code — empty input is an error, not silent success |
| `CorruptReseqFile` | `queryProfile` with truncated/garbage `.reseq` file as input | Non-zero exit code, stderr contains error message, no crash/segfault |

**Note:** The `EmptyBam` contract is pinned here: the tool MUST report an error for empty input, not silently produce empty output. If the current code silently succeeds on empty BAM, the test should document that as a known issue (not fail).

---

### 7d: Multi-Thread Invariant Tests — DEFERRED

The original Phase 7 spec included multi-thread invariant tests (byte-identity across thread counts, statistical equivalence). This is deferred because:

1. **TSan already runs in CI** and catches data races
2. **Clang thread safety annotations** (compile-time) would provide higher ROI — consider for a future phase
3. The statistical equivalence tests require careful threshold tuning and are prone to flakiness

If thread-related bugs surface, revisit this sub-task.

---

## Verification Gate

After all sub-tasks complete:

- [ ] CI coverage job passes (no negative count errors)
- [ ] Codecov integration active: badge in README, PR annotations working
- [ ] `codecov.yml` enforces auto-ratchet + 80% patch coverage
- [ ] `--gtest_shuffle` runs in CI without failures
- [ ] ErrorStats, QualityStats, FragmentDuplicationStats each have >= 3 unit tests
- [ ] BamIngestionEngine and ReadSequenceStats each have >= 2 unit tests
- [ ] E2E error-path tests added (>= 3 new regression tests)
- [ ] `make build && make test` passes

---

## Files Created/Modified

| File | Action |
|------|--------|
| `reseq/CMakeLists.txt` | Add `-fprofile-update=atomic` to coverage options; add new test sources to `reseq_test` target |
| `reseq/test_main.cpp` | Add `#include` and `Register()` calls for new test classes |
| `.github/workflows/ci.yml` | Fix lcov, add Codecov upload, add test-isolation job |
| `codecov.yml` | New — auto-ratchet configuration |
| `README.md` | Add Codecov badge |
| `reseq/ErrorStatsTest.cpp` | Add direct fixture-based tests |
| `reseq/QualityStatsTest.cpp` | Add direct fixture-based tests |
| `reseq/FragmentDuplicationStatsTest.cpp` | Add direct fixture-based tests |
| `reseq/BamIngestionEngineTest.cpp` | New test file |
| `reseq/BamIngestionEngineTest.h` | New test header |
| `reseq/ReadSequenceStatsTest.cpp` | New test file |
| `reseq/ReadSequenceStatsTest.h` | New test header |
| `reseq/RegressionTest.cpp` | Add edge-case E2E tests |

---

## What This Phase Does NOT Do

- No hard coverage threshold (aspirational 70-75%, enforced only via auto-ratchet)
- No multi-thread invariant tests (deferred — TSan covers this)
- No CLI module unit tests (covered by existing 25 regression tests)
- No Clang thread safety annotations (future phase candidate)
- No test parallelization (tests remain serial due to existing `Register()` dependencies)
