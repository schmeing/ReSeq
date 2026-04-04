# CoverageStats Concurrency Modernization — Design Document

**Date:** 2026-04-03
**Phase:** 4a (Tranche 2, Task 6)
**Status:** Draft
**Prerequisite:** Phase 3 (PR #5, `c914db2`) and Tranche 1 (PR #6, `296c273`) merged to master.

---

## 1. Current Architecture

CoverageBlock forms a doubly-linked list: atomic `next_block_` for forward traversal,
raw `previous_block_` for reverse. `first_block_` and `last_block_` are atomic pointers
at the CoverageStats level. Blocks are created by `CreateBlock` (try_lock on `reuse_mutex_`
for pool reuse, fallback to `new`), processed in-place by worker threads, then detached
under `clean_up_mutex_` and returned to pool under `reuse_mutex_`.

Phase 3 already converted the pool to `vector<unique_ptr<CoverageBlock>>` and `reads_`
to `vector<unique_ptr<FullRecord>>`. Tranche 1 converted the DataStats threads that
use CoverageStats to `std::jthread`. CoverageStats itself creates no threads.

### Key Data Structures

```
CoverageStats
  first_block_  : atomic<CoverageBlock*>     — front of live list
  last_block_   : atomic<CoverageBlock*>     — end of live list
  reusable_blocks_ : vector<unique_ptr<CoverageBlock>>  — reuse pool

CoverageBlock
  sequence_id_  : uintRefSeqId
  start_pos_    : uintSeqLen
  previous_block_ : CoverageBlock*           — raw, reverse traversal
  next_block_   : atomic<CoverageBlock*>     — atomic, forward traversal
  coverage_     : vector<CoveragePosition>   — per-position atomic counters
  previous_coverage_ : vector<ProcessedCoveragePosition>
  unprocessed_fragments_ : atomic<uintFragCount>
  reads_        : vector<unique_ptr<FullRecord>>
  scheduled_for_processing_ : atomic_flag
  processed_    : atomic<bool>

CoveragePosition
  coverage_forward_ : array<atomic<uintCovCount>, 5>
  coverage_reverse_ : array<atomic<uintCovCount>, 5>
```

### Mutexes

| Mutex | Purpose | Code Site |
|-------|---------|-----------|
| `clean_up_mutex_` | Protects `first_block_`/`last_block_` pointer updates during detach | `CoverageStats.h:147` |
| `reuse_mutex_` | Protects `reusable_blocks_` pool access | `CoverageStats.h:148` |
| `variant_loading_mutex_` | Serializes reference variant loading | `CoverageStats.h:149` |

### Thread Model

CoverageStats does not create threads. DataStats creates worker threads (already jthread
from Tranche 1, `DataStats.cpp:1307-1314`) that call into CoverageStats via:
- `CreateBlock` — allocate or reuse a block
- `AddFragment` / `RemoveFragment` — track active reads per block
- `AddForward` / `AddReverse` — increment atomic coverage counters
- `CleanUp` — process finished blocks, detach, return to pool
- `FindBlock` — reverse traversal to locate block by position
- `Finalize` — called by last thread, processes all remaining blocks

---

## 2. Invariant Catalog

| # | Invariant | Current Mechanism | Code Site | Safety Status |
|---|-----------|-------------------|-----------|---------------|
| C1 | Forward traversal is lock-free during read processing | `atomic<CoverageBlock*> next_block_` | `CoverageStats.h:93` | Safe |
| C2 | Reverse traversal (FindBlock) only touches live blocks | `clean_up_mutex_` nullifies `previous_block_` before pool return | `CoverageStats.cpp:996` | Safe but fragile — raw pointer |
| C3 | Block creation never blocks on pool contention | `try_lock` on `reuse_mutex_`; fallback to `new` | `CoverageStats.cpp:449` | Safe |
| C4 | Block detach and pool return are two-phase | `clean_up_mutex_` for unlinking, then `reuse_mutex_` for pool | `CoverageStats.cpp:984,1023` | Safe |
| C5 | Coverage increments are lock-free | `atomic<uintCovCount>` per base per strand | `CoverageStats.h:42-44` | Safe |
| C6 | Only one thread processes a block's coverage data | `atomic_flag scheduled_for_processing_` | `CoverageStats.h:100` | Safe |
| C7 | Block not reclaimed until fully processed | `!unprocessed_fragments_ && processed_` guard | `CoverageStats.cpp:988-990` | Safe |
| C8 | Gap handling propagates systematic errors across non-contiguous blocks | `previous_coverage_` populated with dummies for gaps | `CoverageStats.cpp:559-582` | Safe |
| C9 | Variant loading is serialized | `variant_loading_mutex_` with try_lock | `CoverageStats.h:149` | Safe |

### Invariant Details

**C1 — Lock-free forward traversal:** Workers traverse blocks via `next_block_` (atomic load)
during read processing. This must remain lock-free — adding a mutex around traversal would
create a bottleneck in the hot path (coverage increment is per-read, per-position).

**C2 — Reverse traversal safety:** `FindBlock` (`CoverageStats.cpp:827-836`) walks backward
from `last_block_` via `previous_block_` (raw pointer). Safety depends on `CleanUp` nullifying
`previous_block_` of the new first block (`CoverageStats.cpp:996`) before returning detached
blocks to the pool. A raw pointer to a deleted block would be a use-after-free. Phase 3's
`unique_ptr` conversion of the pool prevents double-free but doesn't prevent dangling
`previous_block_` if the nullification order is wrong.

**C3 — Non-blocking creation:** `CreateBlock` (`CoverageStats.cpp:447-472`) uses
`try_to_lock` on `reuse_mutex_`. If the mutex is held (cleanup returning blocks to pool),
creation falls through to `new`. This prevents a deadlock where the creating thread needs
a block but the pool is locked by cleanup.

**C4 — Two-phase detach:** CleanUp has three phases:
1. Process ready blocks (no mutex — `scheduled_for_processing_` atomic_flag serializes)
2. Detach from live list under `clean_up_mutex_` — advance `first_block_`, nullify pointers
3. Count blocks and return to pool under `reuse_mutex_`

The two-mutex separation ensures no thread traverses a block being reclaimed.

**C8 — Gap handling:** When `next_block_` starts at a non-contiguous position,
`ProcessBlock` (`CoverageStats.cpp:559-582`) inserts dummy `ProcessedCoveragePosition`
entries into `previous_coverage_` of the next block. Three cases:
- Adjacent blocks (`start_pos_ + kBlockSize == next->start_pos_`): direct data sharing
- Small gap (gap < `reset_distance_`): systematic errors propagated with dummy padding
- Large gap (gap >= `reset_distance_`): coverage data not propagated, error regions reset

---

## 3. Target Architecture

Replace the intrusive doubly-linked list with a **`std::deque<std::unique_ptr<CoverageBlock>>`**
indexed by block ID, plus an explicit free list.

### Why Deque

Blocks are consumed from the front (cleanup) and appended at the back (creation).
`std::deque` supports O(1) operations at both ends without invalidating pointers/references
to other elements. A `std::vector` would invalidate all references on reallocation.
Pointer stability is critical because worker threads hold `CoverageBlock*` references
during read processing.

### Container Changes

| Current | Target | Notes |
|---------|--------|-------|
| `atomic<CoverageBlock*> first_block_` | `std::atomic<size_t> first_live_idx_` | Index into deque (atomic for cross-thread visibility) |
| `atomic<CoverageBlock*> last_block_` | `std::atomic<size_t> last_live_idx_` | Index into deque (atomic — publication signal for new blocks) |
| `CoverageBlock* previous_block_` | `size_t prev_block_idx_` | `SIZE_MAX` = null |
| `atomic<CoverageBlock*> next_block_` | `size_t next_block_idx_` | `SIZE_MAX` = null |
| `vector<unique_ptr<CoverageBlock>> reusable_blocks_` | `vector<size_t> free_indices_` | Indices of reusable slots |

### CoverageBlock Changes

```cpp
struct CoverageBlock {
    uintRefSeqId sequence_id_;
    uintSeqLen start_pos_;
    size_t prev_block_idx_;               // Replaces previous_block_
    size_t next_block_idx_;               // Replaces next_block_ (atomic no longer needed)
    std::vector<CoveragePosition> coverage_;
    std::vector<ProcessedCoveragePosition> previous_coverage_;
    std::atomic<uintFragCount> unprocessed_fragments_;
    std::vector<std::unique_ptr<FullRecord>> reads_;
    intVariantId first_variant_id_;
    std::atomic_flag scheduled_for_processing_;
    std::atomic<bool> processed_;
};
```

### CoverageStats Changes

```cpp
class CoverageStats {
    std::deque<std::unique_ptr<CoverageBlock>> blocks_;  // Owns all blocks
    std::vector<size_t> free_indices_;                    // Recycled block slots
    std::atomic<size_t> first_live_idx_;                    // Front of live range (atomic)
    std::atomic<size_t> last_live_idx_;                     // End of live range (atomic — publication)

    std::mutex clean_up_mutex_;         // Unchanged
    std::mutex reuse_mutex_;            // Now protects free_indices_
    std::mutex variant_loading_mutex_;  // Unchanged
};
```

### Invariant Preservation

**C1 (lock-free forward traversal):** `next_block_idx_` is a plain `size_t` set before
publication (`last_live_idx_` update). Workers look up `blocks_[next_block_idx_].get()`
to get the block pointer. The deque guarantees pointer stability — the `CoverageBlock*`
returned by `get()` remains valid as long as the `unique_ptr` exists in the deque.

No atomic needed on `next_block_idx_` because: the creator writes it in `CreateBlock`
before publishing the block (i.e., before updating `last_live_idx_`). The publication
of `last_live_idx_` must use `std::atomic<size_t>` with release semantics. Readers
that discover the new block via `last_live_idx_.load(acquire)` are guaranteed to see
the initialized `next_block_idx_`. This is the same happens-before pattern as the
current atomic `next_block_` pointer store — just decomposed into index + atomic
publication flag.

**C2 (reverse traversal safety):** `FindBlock` uses `prev_block_idx_` indexing into the
deque. Detached blocks have `prev_block_idx_ = SIZE_MAX`. No dangling pointer — the
deque element still exists until the index is pushed to `free_indices_`. The free list
is only consumed by `CreateBlock`, which reinitializes the block before use.

**C3 (non-blocking creation):** Keep try_lock pattern for `free_indices_`. If lock
contention, `emplace_back` a new `unique_ptr<CoverageBlock>` at end of deque instead.

**C4 (two-phase detach):** Same two-mutex pattern operating on indices:
1. Under `clean_up_mutex_`: advance `first_live_idx_`, set `prev_block_idx_ = SIZE_MAX`
2. Under `reuse_mutex_`: push freed indices to `free_indices_`

**C5-C7:** Unchanged — these operate on block internals (atomic counters, flags), not
the container structure.

**C8 (gap handling):** Unchanged — gap detection uses block fields (`start_pos_`,
`sequence_id_`), not container structure. `next_block_idx_` lookup replaces pointer
dereference.

**C9 (variant loading):** Unchanged.

### Mutex Changes

No new mutexes. `clean_up_mutex_` and `reuse_mutex_` serve the same roles. `reuse_mutex_`
now protects `free_indices_` instead of `reusable_blocks_`. `variant_loading_mutex_` unchanged.

---

## 4. Migration Plan

Each step compiles and passes the full test suite.

1. **Add deque and index infrastructure.** Add `blocks_` deque, `free_indices_` vector,
   and index fields (`block_id_`, `prev_block_idx_`, `next_block_idx_`) to CoverageBlock
   alongside existing pointer fields. Initialize in `Prepare`. No behavioral change.

2. **Convert CreateBlock to use deque.** Allocate blocks via deque `emplace_back` or
   free list reuse. Set both pointer and index fields. Keep pointer fields as aliases
   derived from `blocks_[idx].get()`.

3. **Convert CleanUp to detach by index.** Under `clean_up_mutex_`, advance
   `first_live_idx_` and set `prev_block_idx_ = SIZE_MAX`. Under `reuse_mutex_`, push
   freed indices. Keep pointer nullification for safety during transition.

4. **Convert FindBlock to index-based traversal.** Walk backward via `prev_block_idx_`
   and deque lookup instead of `previous_block_` pointer.

5. **Convert all remaining pointer-based access to index-based.** Update `AddFragment`,
   `RemoveFragment`, `EnsureSpace`, `ProcessBlock` gap handling, `Finalize`.

6. **Remove pointer fields.** Delete `previous_block_`, `next_block_`, `first_block_`,
   `last_block_`. Remove `reusable_blocks_`. All access is now index-based.

7. **Remove atomicity from next_block.** `next_block_idx_` is plain `size_t` — immutable
   once set, no concurrent writes.

---

## 5. Testing

### New Unit Tests

| Test | What It Validates |
|------|-------------------|
| `CoverageBlockDequeLifecycleTest` | Create blocks via deque, verify index assignment, free blocks, verify reuse from free list |
| `CoverageBlockFindByIndexTest` | Populate deque with blocks across multiple sequences, verify `FindBlock` index-based reverse traversal returns correct block |
| `CoverageBlockGapHandlingTest` | Create non-contiguous blocks, verify `previous_coverage_` dummy insertion for small gaps and reset for large gaps |
| `CoverageBlockCleanupTest` | Process blocks, verify two-phase detach (clean_up_mutex then reuse_mutex), verify freed indices recycled correctly |
| `CoverageBlockTryLockFallbackTest` | Hold `reuse_mutex_` from test thread, call `CreateBlock` from another thread, verify new block allocated (not blocked) |
| `CoverageBlockConcurrentIncrementTest` | Multiple threads increment coverage counters on same block, verify final counts are correct (atomic correctness) |

### Existing Tests

All tests in `CoverageStatsTest` suite plus the 10 golden-file regression tests must
pass unchanged at every migration step.

---

## 6. Verification Criteria

- `make build && make test` — all 65 unit tests + 10 regression tests pass
- TSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=thread"` — no data race warnings
- ASan+UBSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` — no errors
- No performance regression on ecoli test data (block creation/cleanup is not the hot path)
