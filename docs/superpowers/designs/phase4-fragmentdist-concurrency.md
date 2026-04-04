# FragmentDistributionStats Concurrency Modernization — Design Document

**Date:** 2026-04-03
**Phase:** 4a (Tranche 2, Task 8)
**Status:** Draft
**Prerequisite:** Phase 3 (PR #5, `c914db2`) and Tranche 1 (PR #6, `296c273`) merged to master.

---

## 1. Current Architecture

FragmentDistributionStats implements a lock-free producer-consumer bias calculation
pipeline. A fixed-size queue of 100 slots (`claimed_bias_bins_[100]`) is managed via
`atomic_flag` for slot ownership. Workers call `AddNewBiasCalculations` to claim slots
and populate parameters, then `ExecuteBiasCalculations` to process them. Publication
ordering is enforced by setting `current_bias_param_` last.

Temporary computation vectors (`bias_calc_vects_`) use a separate spin-wait acquisition
pattern: each thread scans the deque of flag-guarded vectors until finding a free one.

Two VLA `thread[]` arrays remain unconverted:
- `BiasSumThread` at `FragmentDistributionStats.cpp:2767` — bias sum computation
- `BiasNormalizationThread` at `FragmentDistributionStats.cpp:3723` — normalization

Both use atomic fetch-and-increment work stealing with manual join.

### Key Data Structures

```
FragmentDistributionStats
  // Queue mechanism (100-slot bounded queue)
  claimed_bias_bins_       : array<atomic_flag, 100>         — slot ownership
  current_bias_param_      : array<atomic<uintNumFits>, 100> — per-slot work index
  finished_bias_calcs_     : array<atomic<uintNumFits>, 100> — per-slot completion count
  current_bias_result_     : atomic<uintNumFits>             — converged result counter
  params_left_for_calculation_ : atomic<uintNumFits>         — termination counter
  params_fitted_           : atomic<uintNumFits>             — progress tracking

  // Temporary vector pools (one per thread)
  bias_calc_vects_  : deque<pair<atomic_flag, BiasCalculationVectors>>
  tmp_frag_count_   : deque<pair<atomic_flag, vector<uintFragCount>>>

  // Reference sequence progress
  num_handled_reference_sequence_bins_ : atomic<uintRefSeqBin>

  // Bias calculation parameters
  bias_calc_params_ : vector<BiasCalculationParams>
```

### Queue Mechanism

**AddNewBiasCalculations** (`FragmentDistributionStats.cpp:2957-2985`):
1. Scan `claimed_bias_bins_` via `test_and_set` to find free slot
2. If queue full (`queue_spot >= 100`): **break** (non-blocking guarantee)
3. If slot found: `compare_exchange_strong` on `num_handled_reference_sequence_bins_`
   to claim a reference sequence bin
4. If CAS succeeds: populate parameters via `UpdateBiasCalculationParams`
5. If CAS fails: release slot immediately

**UpdateBiasCalculationParams** (`FragmentDistributionStats.cpp:2169-2248`):
1. Fill `bias_calc_params_` entries for the queue slot
2. Reset `finished_bias_calcs_` to 0 (line 2239)
3. Set `current_bias_param_` to 0 **LAST** (line 2240) — publication signal

**ExecuteBiasCalculations** (`FragmentDistributionStats.cpp:2987-3036`):
1. For each of 100 queue bins: atomic fetch-and-increment `current_bias_param_[bin]`
2. Process claimed parameter: acquire `bias_calc_vects_` slot (spin-wait), compute, release
3. When `++finished_bias_calcs_[bin] == queue_bin_size`: release slot via `claimed_bias_bins_.clear()`

### Mutexes

| Mutex | Purpose | Code Site |
|-------|---------|-----------|
| `file_mutex_` (static) | Protects file output | `FragmentDistributionStats.h:64` |
| `print_mutex` (passed) | Progress logging | Throughout call stack |
| `result_mutex` (local) | Normalization result aggregation | `FragmentDistributionStats.cpp:3718` |

### Thread Model

**BiasSumThread** (`FragmentDistributionStats.cpp:2767-2775`): VLA `thread[]`,
atomic `current_param` work stealing, writes to shared `bias_sum` vector
(race-free because each parameter is processed by exactly one thread).

**BiasNormalizationThread** (`FragmentDistributionStats.cpp:3723-3731`): VLA `thread[]`,
atomic `current_param` work stealing, aggregation under `result_mutex`.

**Main bias pipeline** threads: Created by DataStats (already jthread from Tranche 1).
These threads call into `HandleReferenceSequencesUntil` which calls
`AddNewBiasCalculations` + `ExecuteBiasCalculations`.

---

## 2. Invariant Catalog

| # | Invariant | Current Mechanism | Code Site | Safety Status |
|---|-----------|-------------------|-----------|---------------|
| F1 | Queue full -> break, never block | `test_and_set` scan exits at `queue_spot >= 100` -> `break` | `FragmentDistributionStats.cpp:2962-2982` | Safe — critical |
| F2 | Publication ordering: parameters visible before work index | `current_bias_param_` set to 0 LAST | `FragmentDistributionStats.cpp:2239-2242` | Safe |
| F3 | Slot released only when ALL calculations in bin complete | `++finished_bias_calcs_ == queue_bin_size` guard | `FragmentDistributionStats.cpp:3025-3032` | Safe |
| F4 | Vector acquisition never deadlocks | Exactly `num_threads` vectors for `num_threads` workers; each holds at most one | `FragmentDistributionStats.cpp:3001-3003` | Safe |
| F5 | Work stealing is fair (no starvation) | Atomic fetch-and-increment on `current_bias_param_[bin]` | `FragmentDistributionStats.cpp:2991` | Safe |
| F6 | Termination is guaranteed | `params_left_for_calculation_` decremented to 0; `FinishThreads` loops until zero | `FragmentDistributionStats.cpp:2995, 3317` | Safe |
| F7 | Empty bins released without publication | `claimed_bias_bins_.clear()` immediately, skip `current_bias_param_` | `FragmentDistributionStats.cpp:2244-2246` | Safe |
| F8 | Result aggregation is race-free | `current_bias_result_` atomic gives unique slot; readers wait for thread join | `FragmentDistributionStats.cpp:2627, 3357` | Safe |
| F9 | BiasSumThread: no two threads write same parameter | Atomic `current_param` fetch-and-increment ensures unique assignment | `FragmentDistributionStats.cpp:2767-2775` | Safe |
| F10 | BiasNormalizationThread: result aggregation is mutex-guarded | `result_mutex` protects accumulation | `FragmentDistributionStats.cpp:3115, 3718` | Safe |

### Invariant Details

**F1 — Non-blocking progress guarantee (CRITICAL):** This is the most important invariant.
When `AddNewBiasCalculations` cannot find a free queue slot (all 100 claimed), the function
**breaks out of the loop** (`FragmentDistributionStats.cpp:2982`). It does NOT block, wait,
or spin. The caller (`HandleReferenceSequencesUntil`) then calls `ExecuteBiasCalculations`
which drains existing work. Eventually `FinishThreads` retries.

Why this matters: if ALL worker threads blocked on a full queue, no thread would be
available to execute the queued calculations that would free slots. This would be a
classic deadlock. The break-not-block pattern is a deliberate design that must be preserved
in any redesign.

**F2 — Publication ordering:** `UpdateBiasCalculationParams` (`FragmentDistributionStats.cpp:2169-2248`)
fills the queue slot's parameter entries, then as the very last step sets
`current_bias_param_[queue_spot] = 0` (line 2240). Workers in `ExecuteBiasCalculations`
use `current_bias_param_[bin]++` to claim work items. Setting it to 0 last ensures workers
see fully initialized parameters when they start processing.

Note: this relies on the `atomic` store having release semantics and the `atomic`
fetch-add having acquire semantics (both are defaults for `std::atomic`). No explicit
memory order annotations needed.

**F3 — Completion-guarded release:** In `ExecuteBiasCalculations`
(`FragmentDistributionStats.cpp:3025-3032`), after each calculation completes,
`++finished_bias_calcs_[bin]` is checked against `queue_bin_size`. Only when ALL
calculations in the bin are done does `claimed_bias_bins_[bin].clear()` run. Premature
release would let `AddNewBiasCalculations` overwrite in-use parameter data.

**F4 — Vector pool sizing:** `bias_calc_vects_` is resized to exactly `num_threads`
(line 1907). Each worker acquires at most one vector via spin-wait
(`FragmentDistributionStats.cpp:3001-3003`). Since there are `num_threads` vectors for
`num_threads` workers, and each worker releases its vector before trying to acquire
another, there is always a free vector. The spin-wait is bounded and brief.

**F6 — Guaranteed termination:** `params_left_for_calculation_` is initialized to the
total number of bias calculations needed (line 1898) and decremented in
`ExecuteBiasCalculations` (line 2995). `FinishThreads` (`FragmentDistributionStats.cpp:3311-3323`)
loops with 1-second sleeps until this counter reaches 0. The counter is guaranteed to
reach 0 because: (a) every parameter is either claimed by a worker or available for
claiming, (b) the non-blocking guarantee (F1) ensures workers always make progress,
(c) each parameter is decremented exactly once.

---

## 3. Target Architecture

The lock-free queue design is **already well-engineered**. Unlike CoverageStats and
Simulator where intrusive linked lists create ownership hazards, FragmentDistributionStats
uses fixed-size arrays of atomics — no pointer-based ownership, no dangling references.

The redesign focuses on three improvements:
1. Replace `atomic_flag` slot scan with `std::mutex`-guarded bounded channel
2. Replace spin-wait vector acquisition with thread-index-based pool
3. Replace VLA thread arrays with `std::jthread`

### Bounded Work Queue

Replace the three parallel arrays (`claimed_bias_bins_`, `current_bias_param_`,
`finished_bias_calcs_`) with a `BoundedWorkQueue` class that encapsulates slot
management.

```cpp
struct QueueSlot {
    std::atomic<uintNumFits> current_param{0};   // Work-stealing index
    std::atomic<uintNumFits> finished_count{0};   // Completion counter
    uintNumFits total_params{0};                  // Slot size (set on publish)
    bool published{false};                        // Whether slot has work
};

class BoundedWorkQueue {
    static constexpr size_t kMaxSlots = 100;      // kMaxBinsQueuedForBiasCalc
    std::array<QueueSlot, kMaxSlots> slots_;
    std::vector<size_t> free_indices_;             // Available slot indices
    mutable std::mutex mutex_;                     // Protects free_indices_ only

public:
    // NON-BLOCKING: returns slot index, or SIZE_MAX if full
    size_t try_acquire() {
        std::lock_guard lock(mutex_);
        if (free_indices_.empty()) return SIZE_MAX;
        size_t idx = free_indices_.back();
        free_indices_.pop_back();
        return idx;
    }

    // Called by producer after populating slot parameters
    // Sets current_param to 0 LAST (publication signal)
    void publish(size_t idx, uintNumFits total) {
        slots_[idx].finished_count = 0;
        slots_[idx].total_params = total;
        slots_[idx].published = true;
        slots_[idx].current_param = 0;  // LAST — publication ordering (F2)
    }

    // Release without publication (empty bin)
    void release_empty(size_t idx) {
        std::lock_guard lock(mutex_);
        slots_[idx].published = false;
        free_indices_.push_back(idx);
    }

    // Called when all calculations in slot complete
    void release(size_t idx) {
        std::lock_guard lock(mutex_);
        slots_[idx].published = false;
        free_indices_.push_back(idx);
    }

    QueueSlot& slot(size_t idx) { return slots_[idx]; }
    static constexpr size_t max_slots() { return kMaxSlots; }
};
```

**Key design decisions:**
- `try_acquire()` returns `SIZE_MAX` when full — preserves non-blocking guarantee (F1)
- Mutex protects only the `free_indices_` list, not the slot data
- `current_param` and `finished_count` remain `std::atomic` — workers access them
  without holding the queue mutex (work stealing is lock-free within a slot)
- Publication ordering preserved: `current_param = 0` is the last write in `publish()`

### Thread-Index-Based Pool

Replace `deque<pair<atomic_flag, BiasCalculationVectors>>` with direct indexing:

```cpp
std::vector<BiasCalculationVectors> bias_calc_vects_;   // One per thread
std::vector<std::vector<uintFragCount>> tmp_frag_count_; // One per thread
```

Each jthread lambda captures its thread index `i` (0..num_threads-1). Thread `i`
accesses `bias_calc_vects_[i]` exclusively. No atomic flags needed — ownership is
structural, not dynamic.

### Container Changes Summary

| Current | Target |
|---------|--------|
| `array<atomic_flag, 100> claimed_bias_bins_` | `BoundedWorkQueue::free_indices_` + `mutex_` |
| `array<atomic<uintNumFits>, 100> current_bias_param_` | `BoundedWorkQueue::slots_[].current_param` |
| `array<atomic<uintNumFits>, 100> finished_bias_calcs_` | `BoundedWorkQueue::slots_[].finished_count` |
| `deque<pair<atomic_flag, BiasCalculationVectors>>` | `vector<BiasCalculationVectors>` (thread-indexed) |
| `deque<pair<atomic_flag, vector<uintFragCount>>>` | `vector<vector<uintFragCount>>` (thread-indexed) |
| `thread threads[num_threads]` (line 2767) | `vector<jthread>` |
| `thread threads[num_threads]` (line 3723) | `vector<jthread>` |

### Thread Changes

Replace both VLA `thread[]` arrays with `vector<jthread>`:

```cpp
// BiasSumThread (FragmentDistributionStats.cpp:2767)
{
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (decltype(num_threads) i = 0; i < num_threads; ++i) {
        threads.emplace_back([this, &reference, &params, &current_param,
                              &finished_params, &bias_sum, &print_mutex](std::stop_token) {
            BiasSumThread(*this, reference, params, current_param,
                          finished_params, bias_sum, print_mutex);
        });
    }
    // jthread destructors join on scope exit
}

// BiasNormalizationThread (FragmentDistributionStats.cpp:3723) — same pattern
```

Workers accept `stop_token` but do NOT check it — natural completion semantics.

For `ExecuteBiasCalculations`, the thread index is passed to the lambda so the worker
can access `bias_calc_vects_[thread_idx]` directly:

```cpp
// In DataStats ReadThread, thread index is passed through ThreadData
void ExecuteBiasCalculations(const Reference& reference,
                            FragmentDuplicationStats& duplications,
                            std::mutex& print_mutex,
                            size_t thread_idx) {
    // Use bias_calc_vects_[thread_idx] directly — no spin-wait
    auto& calc_vect = bias_calc_vects_[thread_idx];
    // ...
}
```

### Invariant Preservation

**F1 (non-blocking):** `try_acquire()` returns `SIZE_MAX` when `free_indices_` empty.
Caller breaks. Same semantics, O(1) instead of O(100) scan.

**F2 (publication ordering):** `publish()` sets `current_param` to 0 as the last step.
Workers read `current_param` via atomic fetch-and-increment without holding the queue
mutex — same acquire/release semantics as the current code.

**F3 (slot release):** `release()` called only when `++finished_count == total_params`.
Pushes index back to `free_indices_` under mutex.

**F4 (vector pool):** Thread-index-based — each thread owns exactly one vector by
construction. No acquisition needed, no contention possible. Strictly better than
current spin-wait pattern.

**F5 (fair work stealing):** Atomic fetch-and-increment on `current_param` per slot —
unchanged.

**F6 (termination):** `params_left_for_calculation_` unchanged. `FinishThreads` loop
unchanged.

**F7 (empty bin release):** `try_acquire()` + immediate `release_empty()` without
`publish()`.

**F8 (result aggregation):** Unchanged — atomic `current_bias_result_` gives unique slots.

**F9/F10 (thread-specific aggregation):** Unchanged.

---

## 4. Migration Plan

Each step compiles and passes the full test suite.

1. **Add BoundedWorkQueue class.** New file `reseq/BoundedWorkQueue.h` (header-only,
   ~60 lines). Add to `reseq_lib` sources. Compile — no behavioral change.

2. **Add BoundedWorkQueue member to FragmentDistributionStats.** Initialize in `Prepare`.
   Keep old arrays alongside for dual-path validation.

3. **Convert AddNewBiasCalculations to use BoundedWorkQueue.** Replace `test_and_set`
   scan with `try_acquire()`. Replace `UpdateBiasCalculationParams` completion with
   `publish()`. Replace empty-bin release with `release_empty()`.

4. **Convert ExecuteBiasCalculations to use BoundedWorkQueue.** Access slots via
   `queue_.slot(bin)`. Replace `claimed_bias_bins_.clear()` with `queue_.release()`.

5. **Remove old atomic arrays.** Delete `claimed_bias_bins_`, `current_bias_param_`,
   `finished_bias_calcs_`. All access goes through `BoundedWorkQueue`.

6. **Convert bias_calc_vects_ to thread-indexed vector.** Pass thread index through
   `ThreadData` or lambda capture. Remove `atomic_flag` from vector pairs. Resize to
   `num_threads` in `Prepare`.

7. **Convert tmp_frag_count_ to thread-indexed vector.** Same pattern as step 6.

8. **Replace VLA thread arrays with vector<jthread>.** Convert both sites
   (lines 2767 and 3723). Pass thread index to lambda for BiasSumThread and
   BiasNormalizationThread.

---

## 5. Testing

### New Unit Tests

| Test | What It Validates |
|------|-------------------|
| `BoundedWorkQueueTest.TryAcquireWhenFull` | Fill all 100 slots, verify `try_acquire()` returns `SIZE_MAX` |
| `BoundedWorkQueueTest.AcquirePublishRelease` | Acquire slot, publish with parameters, verify workers see `current_param == 0`, complete all work, release, verify slot recycled |
| `BoundedWorkQueueTest.ConcurrentAcquireRelease` | N threads racing `try_acquire()` + `release()`, verify no double-acquire, no lost slots, all slots eventually returned |
| `BoundedWorkQueueTest.ReleaseEmpty` | Acquire + release_empty without publish, verify slot immediately reusable |
| `FragmentDistributionStats.NonBlockingProgressGuarantee` | Simulate full queue scenario with limited threads, verify all threads make forward progress (no deadlock within timeout) |
| `FragmentDistributionStats.ThreadIndexPoolExclusiveAccess` | Verify each thread accesses only its own `bias_calc_vects_[i]` — no cross-thread access |
| `FragmentDistributionStats.JthreadShutdown` | Verify clean shutdown: jthread vector goes out of scope, all threads join, no hanging |
| `FragmentDistributionStats.BiasPublicationOrdering` | Use test harness that reads slot state immediately when `current_param` transitions to 0 — verify all parameter fields are populated |
| `FragmentDistributionStats.TerminationGuarantee` | Run full bias pipeline, verify `params_left_for_calculation_` reaches 0 and `FinishThreads` exits |

### Existing Tests

All tests in `FragmentDistributionStatsTest` suite plus the 10 golden-file regression
tests must pass unchanged at every migration step.

---

## 6. Verification Criteria

- `make build && make test` — all 65 unit tests + 10 regression tests pass
- TSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=thread"` — no data race warnings
- ASan+UBSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` — no errors
- No performance regression on ecoli test data
- BoundedWorkQueue unit tests pass under TSan (concurrent acquire/release stress test)
- Spin-wait elimination verified: no `test_and_set` loops remain in production code
  (replaced by thread-index pool or mutex-guarded free list)
