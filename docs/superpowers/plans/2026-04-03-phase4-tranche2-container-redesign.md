# Phase 4 Tranche 2: Complex Subsystem Container Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace intrusive linked lists in CoverageStats, Simulator, and FragmentDistributionStats with indexed owning containers (`std::deque`), adopt `std::jthread` for remaining VLA thread arrays, and add unit tests for all modernized concurrency patterns.

**Architecture:** Three subsystems are modernized independently. CoverageStats and Simulator replace pointer-based linked lists with `std::deque<std::unique_ptr<>>` / `std::deque<>` indexed by `size_t`. FragmentDistributionStats replaces its `atomic_flag`-based queue with a `BoundedWorkQueue` class and converts spin-wait vector pools to thread-indexed access. Each subsystem follows a dual-representation migration: new index fields are added alongside old pointer fields, code is converted function-by-function, then old fields are removed.

**Tech Stack:** C++20 (`std::jthread`, `std::stop_token`, `std::atomic`, `std::scoped_lock`), GoogleTest, CMake, TSan, ASan+UBSan.

**Prerequisite:** Phase 3 (PR #5, `c914db2`) and Tranche 1 (PR #6, `296c273`) are on master. `make build && make test` passes (65 tests + 10 regression tests).

**Design documents:**
- `docs/superpowers/designs/phase4-coveragestats-concurrency.md`
- `docs/superpowers/designs/phase4-simulator-concurrency.md`
- `docs/superpowers/designs/phase4-fragmentdist-concurrency.md`

---

## Task 1: CoverageStats — Add Deque Infrastructure and Index Fields

**Files:**
- Modify: `reseq/CoverageStats.h:89-109` (CoverageBlock struct), `reseq/CoverageStats.h:280-282` (member declarations), `reseq/CoverageStats.h:395-397` (constructor)
- Modify: `reseq/CoverageStats.cpp:762-825` (Prepare function)

This task adds the new deque + free list alongside existing pointer fields. No behavioral change — old code paths remain active.

- [ ] **Step 1: Add index fields to CoverageBlock**

In `reseq/CoverageStats.h`, add three fields to the CoverageBlock struct (after `processed_`, before the constructor):

```cpp
size_t block_idx_;          // Index of this block in blocks_ deque
size_t prev_block_idx_;     // Index of previous block (SIZE_MAX = null)
size_t next_block_idx_;     // Index of next block (SIZE_MAX = null)
```

Update the CoverageBlock constructor to initialize these:

```cpp
CoverageBlock(uintRefSeqId seq_id, uintSeqLen start_pos, CoverageBlock* prev_block)
    : sequence_id_(seq_id), start_pos_(start_pos), previous_block_(prev_block), next_block_(nullptr),
      unprocessed_fragments_(0), first_variant_id_(0),
      block_idx_(SIZE_MAX), prev_block_idx_(SIZE_MAX), next_block_idx_(SIZE_MAX) {
    scheduled_for_processing_.clear();
    processed_ = false;
}
```

- [ ] **Step 2: Add deque and free list to CoverageStats**

In `reseq/CoverageStats.h`, add after `reusable_blocks_` (line 282):

```cpp
std::deque<std::unique_ptr<CoverageBlock>> blocks_;  // Indexed block storage
std::vector<size_t> free_indices_;                    // Recycled block slots
std::atomic<size_t> first_live_idx_{SIZE_MAX};        // Front of live range
std::atomic<size_t> last_live_idx_{SIZE_MAX};         // End of live range (publication signal)
```

Add `#include <deque>` to the include block if not already present.

- [ ] **Step 3: Initialize deque in Prepare**

In `reseq/CoverageStats.cpp`, in the `Prepare` function (around line 764, after `reusable_blocks_.reserve(100)`), add:

```cpp
blocks_.clear();
free_indices_.clear();
first_live_idx_ = SIZE_MAX;
last_live_idx_ = SIZE_MAX;
```

- [ ] **Step 4: Build and run full test suite**

```bash
make build && make test
```

Expected: all 65 unit tests + 10 regression tests pass. No behavioral change.

- [ ] **Step 5: Commit**

```bash
git add reseq/CoverageStats.h reseq/CoverageStats.cpp
git commit -m "refactor(4a): add deque infrastructure and index fields to CoverageStats

Add blocks_ deque, free_indices_ vector, first_live_idx_/last_live_idx_
atomic indices to CoverageStats. Add block_idx_, prev_block_idx_,
next_block_idx_ to CoverageBlock. No behavioral change — old pointer
fields remain active. Dual-representation migration step 1."
```

---

## Task 2: CoverageStats — Convert CreateBlock to Use Deque

**Files:**
- Modify: `reseq/CoverageStats.cpp:447-480` (CreateBlock)

Convert block allocation to use the deque + free list. Set both pointer AND index fields so old code paths still work.

- [ ] **Step 1: Rewrite CreateBlock to use deque**

Replace the CreateBlock function body (`reseq/CoverageStats.cpp:447-480`) with:

```cpp
CoverageBlock* CoverageStats::CreateBlock(uintRefSeqId seq_id, uintSeqLen start_pos) {
    CoverageBlock* new_block;
    size_t new_idx;

    if (std::unique_lock lock(reuse_mutex_, std::try_to_lock); lock.owns_lock()) {
        if (!free_indices_.empty()) {
            new_idx = free_indices_.back();
            free_indices_.pop_back();
            lock.unlock();

            new_block = blocks_[new_idx].get();
            new_block->sequence_id_ = seq_id;
            new_block->start_pos_ = start_pos;
            new_block->previous_block_ = last_block_;
            new_block->next_block_ = nullptr;
            new_block->coverage_.clear();
            new_block->previous_coverage_.clear();
            new_block->reads_.clear();
            new_block->scheduled_for_processing_.clear();
            new_block->processed_ = false;
            new_block->block_idx_ = new_idx;
            new_block->prev_block_idx_ = last_live_idx_.load(std::memory_order_relaxed);
            new_block->next_block_idx_ = SIZE_MAX;
        } else {
            lock.unlock();
            blocks_.emplace_back(std::make_unique<CoverageBlock>(seq_id, start_pos, last_block_));
            new_idx = blocks_.size() - 1;
            new_block = blocks_[new_idx].get();
            new_block->block_idx_ = new_idx;
            new_block->prev_block_idx_ = last_live_idx_.load(std::memory_order_relaxed);
            new_block->next_block_idx_ = SIZE_MAX;
            new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
        }
    } else {
        blocks_.emplace_back(std::make_unique<CoverageBlock>(seq_id, start_pos, last_block_));
        new_idx = blocks_.size() - 1;
        new_block = blocks_[new_idx].get();
        new_block->block_idx_ = new_idx;
        new_block->prev_block_idx_ = last_live_idx_.load(std::memory_order_relaxed);
        new_block->next_block_idx_ = SIZE_MAX;
        new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
    }

    new_block->reads_.reserve((*last_block_).reads_.capacity());
    new_block->coverage_.resize(kBlockSize);
    new_block->first_variant_id_ = (*last_block_).first_variant_id_;

    // Update old pointer chain (backward compat during migration)
    (*last_block_).next_block_ = new_block;

    // Update index chain
    if (last_live_idx_.load(std::memory_order_relaxed) != SIZE_MAX) {
        blocks_[last_live_idx_.load(std::memory_order_relaxed)]->next_block_idx_ = new_idx;
    }
    last_live_idx_.store(new_idx, std::memory_order_release);
    if (first_live_idx_.load(std::memory_order_relaxed) == SIZE_MAX) {
        first_live_idx_.store(new_idx, std::memory_order_release);
    }

    return new_block;
}
```

- [ ] **Step 2: Build and run full test suite**

```bash
make build && make test
```

Expected: all tests pass. Blocks now live in the deque but old pointer traversal still works.

- [ ] **Step 3: Commit**

```bash
git add reseq/CoverageStats.cpp
git commit -m "refactor(4a): convert CoverageStats::CreateBlock to use deque

Blocks are now allocated in the blocks_ deque or reused from free_indices_.
Both pointer fields (for backward compat) and index fields are maintained.
try_lock pattern preserved for non-blocking creation."
```

---

## Task 3: CoverageStats — Convert CleanUp to Index-Based Detach

**Files:**
- Modify: `reseq/CoverageStats.cpp:965-1032` (CleanUp)
- Modify: `reseq/CoverageStats.cpp:755-760` (RemoveBlock)

- [ ] **Step 1: Update RemoveBlock to push to free list**

Replace the RemoveBlock function (`reseq/CoverageStats.cpp:755-760`) with:

```cpp
CoverageBlock* CoverageStats::RemoveBlock(CoverageBlock* block) {
    CoverageBlock* next = block->next_block_;
    free_indices_.push_back(block->block_idx_);
    // Keep block alive in deque — will be reused via free_indices_
    return next;
}
```

Note: the old code moved blocks into `reusable_blocks_`. Now we push the index to `free_indices_` instead. The block stays in `blocks_` deque and is reused via `CreateBlock`.

- [ ] **Step 2: Update CleanUp phase 2 to maintain index chain**

In the CleanUp function (`reseq/CoverageStats.cpp`), in the phase 2 section under `clean_up_mutex_`, after `first_block_ = until_block;` add:

```cpp
first_live_idx_.store(until_block->block_idx_, std::memory_order_release);
```

And after `(*first_block_).previous_block_ = nullptr;` add:

```cpp
blocks_[first_live_idx_.load(std::memory_order_relaxed)]->prev_block_idx_ = SIZE_MAX;
```

- [ ] **Step 3: Build and run full test suite**

```bash
make build && make test
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/CoverageStats.cpp
git commit -m "refactor(4a): convert CoverageStats cleanup to use index-based free list

RemoveBlock now pushes block_idx_ to free_indices_ instead of moving to
reusable_blocks_. CleanUp maintains first_live_idx_ alongside first_block_
pointer. Dual-representation still active."
```

---

## Task 4: CoverageStats — Convert FindBlock to Index-Based Traversal

**Files:**
- Modify: `reseq/CoverageStats.cpp:827-837` (FindBlock)

- [ ] **Step 1: Rewrite FindBlock to use index traversal**

Replace the FindBlock function body with:

```cpp
CoverageStats::CoverageBlock* CoverageStats::FindBlock(uintRefSeqId ref_seq_id, uintSeqLen ref_pos) {
    size_t idx = last_live_idx_.load(std::memory_order_acquire);
    CoverageBlock* block = blocks_[idx].get();
    while (block->sequence_id_ > ref_seq_id) {
        idx = block->prev_block_idx_;
        block = blocks_[idx].get();
    }
    while (block->start_pos_ > ref_pos) {
        idx = block->prev_block_idx_;
        block = blocks_[idx].get();
    }
    return block;
}
```

- [ ] **Step 2: Build and run full test suite**

```bash
make build && make test
```

Expected: all tests pass. FindBlock now uses index-based traversal but returns the same `CoverageBlock*`.

- [ ] **Step 3: Commit**

```bash
git add reseq/CoverageStats.cpp
git commit -m "refactor(4a): convert CoverageStats::FindBlock to index-based traversal

FindBlock now walks backward via prev_block_idx_ and deque lookup instead
of previous_block_ raw pointer. Returns same CoverageBlock* for callers."
```

---

## Task 5: CoverageStats — Convert Remaining Access and Remove Old Fields

**Files:**
- Modify: `reseq/CoverageStats.h` (remove old pointer fields)
- Modify: `reseq/CoverageStats.cpp` (update Finalize, Prepare, all remaining sites)

This is the final CoverageStats migration step. Convert all remaining pointer-based access to index-based, then remove the old pointer fields.

- [ ] **Step 1: Audit all remaining pointer usage sites**

Search for all sites that use `previous_block_`, `next_block_` (as pointer), `first_block_`, `last_block_`, `reusable_blocks_`:

```bash
cd /home/bernt-popp/development/ReSeq && grep -n 'previous_block_\|->next_block_\|first_block_\|last_block_\|reusable_blocks_' reseq/CoverageStats.cpp reseq/CoverageStats.h | grep -v '//' | head -60
```

Convert each site to use index-based access (`prev_block_idx_`, `next_block_idx_`, `first_live_idx_`, `last_live_idx_`, `blocks_[]`, `free_indices_`). Key sites:

- `Prepare` initialization: replace `first_block_ = new CoverageBlock(...)` with deque emplace
- `AddFragment`/`RemoveFragment`: replace `last_block_` pointer reads with `blocks_[last_live_idx_]`
- `EnsureSpace`: replace `last_block_` and `next_block_` pointer navigation
- `ProcessBlock` gap handling (lines 559-582): replace `next_block_` pointer with `next_block_idx_`
- `Finalize`: replace pointer traversal with index traversal
- All sites reading `(*last_block_)` → `*blocks_[last_live_idx_.load()].get()`

- [ ] **Step 2: Remove old pointer fields from CoverageBlock**

In `reseq/CoverageStats.h`, remove from CoverageBlock:
- `CoverageBlock* previous_block_;`
- `std::atomic<CoverageBlock*> next_block_;`

Update the constructor to remove the `prev_block` parameter and `previous_block_`/`next_block_` initializers:

```cpp
CoverageBlock(uintRefSeqId seq_id, uintSeqLen start_pos, size_t prev_idx)
    : sequence_id_(seq_id), start_pos_(start_pos),
      unprocessed_fragments_(0), first_variant_id_(0),
      block_idx_(SIZE_MAX), prev_block_idx_(prev_idx), next_block_idx_(SIZE_MAX) {
    scheduled_for_processing_.clear();
    processed_ = false;
}
```

- [ ] **Step 3: Remove old pointer fields from CoverageStats**

In `reseq/CoverageStats.h`, remove:
- `std::atomic<CoverageBlock*> first_block_;`
- `std::atomic<CoverageBlock*> last_block_;`
- `std::vector<std::unique_ptr<CoverageBlock>> reusable_blocks_;`

Update the constructor in `CoverageStats.h` to remove `first_block_(nullptr), last_block_(nullptr)`.

- [ ] **Step 4: Build and fix all compilation errors**

```bash
make build 2>&1 | head -80
```

Fix any remaining references to removed fields. This may require several iterations.

- [ ] **Step 5: Run full test suite**

```bash
make test
```

Expected: all 65 unit tests + 10 regression tests pass.

- [ ] **Step 6: Run sanitizers**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 7: Commit**

```bash
git add reseq/CoverageStats.h reseq/CoverageStats.cpp
git commit -m "refactor(4a): complete CoverageStats container migration to deque

Remove previous_block_, next_block_ pointers, first_block_/last_block_
atomic pointers, and reusable_blocks_ vector. All access now through
blocks_ deque indexed by size_t. Free list replaces reuse pool.
Pointer stability guaranteed by std::deque."
```

---

## Task 6: CoverageStats — Add Unit Tests for New Container

**Files:**
- Modify: `reseq/CoverageStatsTest.cpp` (add new test cases)
- Modify: `reseq/CoverageStatsTest.h` (add test declarations if needed)

- [ ] **Step 1: Add CoverageBlockDequeLifecycleTest**

Add to `reseq/CoverageStatsTest.cpp`:

```cpp
TEST_F(CoverageStatsTest, DequeLifecycle) {
    // Verify blocks are created in deque and indices assigned correctly
    // Setup: call Prepare with test parameters, then create several blocks
    // Assert: block_idx_ matches position in deque, free_indices_ empty initially
    // Free a block, assert it appears in free_indices_
    // Create another block, assert it reuses the freed index
}
```

Implement the test body using the friend class access (CoverageStatsTest is already a friend of CoverageStats). Use the test reference data from `BasicTestClassWithReference`.

- [ ] **Step 2: Add CoverageBlockFindByIndexTest**

```cpp
TEST_F(CoverageStatsTest, FindByIndex) {
    // Create blocks across two reference sequences with different start positions
    // Call FindBlock with various (seq_id, pos) combinations
    // Assert: correct block returned for each query
    // Assert: edge cases — first block, last block, position at block boundary
}
```

- [ ] **Step 3: Add CoverageBlockCleanupTest**

```cpp
TEST_F(CoverageStatsTest, CleanupRecyclesIndices) {
    // Create blocks, mark some as processed (unprocessed_fragments_ = 0, processed_ = true)
    // Call CleanUp
    // Assert: processed blocks' indices appear in free_indices_
    // Assert: unprocessed blocks' indices do NOT appear in free_indices_
    // Assert: first_live_idx_ advanced past processed blocks
}
```

- [ ] **Step 4: Add CoverageBlockConcurrentIncrementTest**

```cpp
TEST_F(CoverageStatsTest, ConcurrentCoverageIncrement) {
    // Create a single block with coverage vector
    // Launch N threads, each incrementing coverage counters for different positions
    // Join all threads
    // Assert: final coverage counts equal expected totals (atomic correctness)
}
```

- [ ] **Step 5: Build and run tests**

```bash
make build && make test
```

Expected: all existing tests pass + new tests pass.

- [ ] **Step 6: Commit**

```bash
git add reseq/CoverageStatsTest.cpp reseq/CoverageStatsTest.h
git commit -m "test(4a): add unit tests for CoverageStats deque container

Tests cover: deque lifecycle (create/free/reuse), FindBlock index-based
traversal, cleanup index recycling, concurrent coverage increment."
```

---

## Task 7: Simulator — Add Deque Infrastructure and Index Fields

**Files:**
- Modify: `reseq/Simulator.h:116-145` (SimBlock, SimUnit structs), `reseq/Simulator.h:276-282` (member declarations)

- [ ] **Step 1: Add index fields to SimBlock**

In `reseq/Simulator.h`, add to SimBlock struct (after `first_methylation_id_`, before constructor):

```cpp
size_t block_idx_;          // Index in blocks_ deque
size_t next_block_idx_;     // Index of next block (SIZE_MAX = null)  — will replace atomic next_block_
size_t partner_block_idx_;  // Index of partner block (SIZE_MAX = null) — will replace partner_block_
```

Note: name these `next_block_idx_` and `partner_block_idx_` to avoid collision with existing `next_block_` and `partner_block_`. Update the constructor:

```cpp
SimBlock(uintRefSeqBin id, uintSeqLen start_pos, SimBlock* partner_block, uintSeed seed)
    : id_(id), start_pos_(start_pos), finished_(false), next_block_(nullptr), partner_block_(partner_block),
      seed_(seed), first_variant_id_(0), first_methylation_id_(0),
      block_idx_(SIZE_MAX), next_block_idx_(SIZE_MAX), partner_block_idx_(SIZE_MAX) {}
```

- [ ] **Step 2: Add index fields to SimUnit**

Add to SimUnit struct (after `next_unit_`, before constructor):

```cpp
size_t unit_idx_;            // Index in units_ deque
size_t first_block_idx_;     // Index of first block in chain — will replace first_block_
size_t last_block_idx_;      // Index of last block in chain — will replace last_block_
size_t next_unit_idx_;       // Index of next unit — will replace next_unit_
```

Note: again suffix with `_idx_` to distinguish from existing `first_block_` etc. Update constructor:

```cpp
SimUnit(uintRefSeqId ref_seq_id)
    : ref_seq_id_(ref_seq_id), first_block_(nullptr), last_block_(nullptr), next_unit_(nullptr),
      unit_idx_(SIZE_MAX), first_block_idx_(SIZE_MAX), last_block_idx_(SIZE_MAX), next_unit_idx_(SIZE_MAX) {}
```

- [ ] **Step 3: Add deque containers to Simulator**

In `reseq/Simulator.h`, add after `req_deletion_buffer_` (line 282):

```cpp
std::deque<SimBlock> blocks_;                 // Indexed block storage (value semantics)
std::deque<SimUnit> units_;                   // Indexed unit storage (value semantics)
std::vector<size_t> free_block_indices_;       // Recycled block slots
std::vector<size_t> free_unit_indices_;        // Recycled unit slots
size_t first_unit_idx_{SIZE_MAX};
size_t last_unit_idx_{SIZE_MAX};
size_t current_unit_idx_{SIZE_MAX};
size_t current_block_idx_{SIZE_MAX};
```

Add `#include <deque>` to the include block if not present.

- [ ] **Step 4: Build and run full test suite**

```bash
make build && make test
```

Expected: all tests pass. No behavioral change.

- [ ] **Step 5: Commit**

```bash
git add reseq/Simulator.h
git commit -m "refactor(4a): add deque infrastructure and index fields to Simulator

Add blocks_/units_ deques, free index vectors, and _idx_ fields to
SimBlock/SimUnit alongside existing pointer fields. No behavioral change."
```

---

## Task 8: Simulator — Convert CreateUnit and CreateBlock to Use Deques

**Files:**
- Modify: `reseq/Simulator.cpp:932-1041` (CreateUnit)
- Modify: `reseq/Simulator.cpp:1210-1339` (CreateBlock)

This is the largest single task. Both functions allocate blocks/units and must set both pointer and index fields during the dual-representation phase.

- [ ] **Step 1: Convert CreateUnit to emplace into deques**

In `reseq/Simulator.cpp`, in `CreateUnit` (line 932+), for each `new SimUnit(...)` call, replace with deque emplacement and set index fields. The unit is created at line 934:

```cpp
// Old: unit = new SimUnit(ref_id);
// New:
units_.emplace_back(ref_id);
size_t unit_idx = units_.size() - 1;
unit = &units_[unit_idx];
unit->unit_idx_ = unit_idx;
```

For each `new SimBlock(...)` in the reverse block creation loop (lines 938-945), replace with:

```cpp
// Old: auto block = new SimBlock(first_block_id, 0, nullptr, block_seed_gen_());
// New:
blocks_.emplace_back(first_block_id, 0, nullptr, block_seed_gen_());
size_t blk_idx = blocks_.size() - 1;
auto* block = &blocks_[blk_idx];
block->block_idx_ = blk_idx;
```

For each subsequent `new SimBlock` in the loop, do the same. After linking `next_block_` pointers (line 943), also set index fields:

```cpp
block->next_block_idx_ = old_block->block_idx_;      // reverse direction
old_block->partner_block_idx_ = block->block_idx_;    // forward-to-reverse reference
```

Set `first_reverse_block`'s `block_idx_` for use by CreateBlock's partner linkage.

- [ ] **Step 2: Convert CreateBlock to use deque + free list**

In the CreateBlock function (line 1210+), for the three cases (first block, next unit, same unit), replace `new SimBlock(...)` with deque emplacement or free list reuse:

```cpp
// Allocation pattern:
size_t blk_idx;
if (!free_block_indices_.empty()) {
    blk_idx = free_block_indices_.back();
    free_block_indices_.pop_back();
    blocks_[blk_idx] = SimBlock(id, start_pos, partner, seed);
} else {
    blocks_.emplace_back(id, start_pos, partner, seed);
    blk_idx = blocks_.size() - 1;
}
auto* block = &blocks_[blk_idx];
block->block_idx_ = blk_idx;
```

After each pointer link (`unit->last_block_->next_block_ = block`), also set the index:

```cpp
blocks_[unit->last_block_->block_idx_].next_block_idx_ = blk_idx;
```

For the partner linkage double-hop at line 1270, also set:

```cpp
block->partner_block_idx_ = blocks_[blocks_[unit->last_block_->block_idx_].partner_block_idx_].partner_block_idx_;
```

- [ ] **Step 3: Build and run full test suite**

```bash
make build && make test
```

Expected: all tests pass. Blocks now live in deques but old pointer traversal still works.

- [ ] **Step 4: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "refactor(4a): convert Simulator CreateUnit/CreateBlock to use deques

Blocks and units are emplaced into deques. Index fields maintained
alongside pointer fields. Free list used for block reuse in CreateBlock."
```

---

## Task 9: Simulator — Convert Cleanup and GetNextBlock to Index-Based

**Files:**
- Modify: `reseq/Simulator.cpp:1308-1339` (cleanup in CreateBlock)
- Modify: `reseq/Simulator.cpp:1341-1396` (GetNextBlock)
- Modify: `reseq/Simulator.cpp:2995-3018` (final cleanup)

- [ ] **Step 1: Convert cleanup in CreateBlock to mark-as-free**

In the cleanup loop inside CreateBlock (lines 1308-1330), replace `delete` calls with free list pushes:

```cpp
// Old:
// delete del_block->partner_block_;
// delete del_block;

// New:
free_block_indices_.push_back(del_block->partner_block_->block_idx_);
free_block_indices_.push_back(del_block->block_idx_);
```

For unit deletion:

```cpp
// Old: delete del_unit;
// New: free_unit_indices_.push_back(del_unit->unit_idx_);
```

Keep the pointer-based traversal intact during this step (we still have both representations).

- [ ] **Step 2: Convert GetNextBlock to maintain index fields**

In GetNextBlock (lines 1341-1396), after `current_block_` and `current_unit_` pointer advances (lines 1380-1393), also advance the index versions:

```cpp
current_block_idx_ = block->block_idx_;
current_unit_idx_ = unit->unit_idx_;

// When advancing:
if (current_block_->next_block_) {
    current_block_ = current_block_->next_block_;
    current_block_idx_ = current_block_->block_idx_;
} else {
    current_unit_ = current_unit_->next_unit_;
    current_unit_idx_ = current_unit_ ? current_unit_->unit_idx_ : SIZE_MAX;
    if (current_unit_) {
        current_block_ = current_unit_->first_block_;
        current_block_idx_ = current_block_ ? current_block_->block_idx_ : SIZE_MAX;
    }
}
```

- [ ] **Step 3: Convert final cleanup to mark-as-free**

In the final cleanup code (lines 2995-3018), replace `delete` calls with free list pushes, similar to step 1. At the very end, clear deques:

```cpp
// After the cleanup loop:
blocks_.clear();
units_.clear();
free_block_indices_.clear();
free_unit_indices_.clear();
first_unit_idx_ = SIZE_MAX;
last_unit_idx_ = SIZE_MAX;
```

Keep `first_unit_ = nullptr; last_unit_ = nullptr;` for pointer-based compat.

- [ ] **Step 4: Build and run full test suite**

```bash
make build && make test
```

- [ ] **Step 5: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "refactor(4a): convert Simulator cleanup and GetNextBlock to index-based

Cleanup marks blocks/units as free via index push instead of delete.
GetNextBlock maintains both pointer and index dispatch. Final cleanup
clears deques after mark-as-free loop."
```

---

## Task 10: Simulator — Remove Old Pointer Fields and Adopt jthread

**Files:**
- Modify: `reseq/Simulator.h` (remove pointer fields, add `<thread>` include)
- Modify: `reseq/Simulator.cpp` (convert all remaining pointer access, replace VLA threads)

- [ ] **Step 1: Convert all remaining pointer-based access to index-based**

Search for remaining pointer usage:

```bash
grep -n 'first_block_\|last_block_\|next_block_\|partner_block_\|first_unit_\|last_unit_\|current_unit_\|current_block_\|next_unit_' reseq/Simulator.cpp reseq/Simulator.h | grep -v '_idx_\|//' | head -80
```

Convert each site. Key patterns:
- `block->next_block_` → `blocks_[block->next_block_idx_]`
- `block->partner_block_` → `blocks_[block->partner_block_idx_]`
- `unit->first_block_` → `blocks_[unit->first_block_idx_]`
- `first_unit_` → `units_[first_unit_idx_]`

- [ ] **Step 2: Remove old pointer fields from SimBlock**

Remove from SimBlock in `Simulator.h`:
- `std::atomic<SimBlock*> next_block_;`
- `SimBlock* partner_block_;`

Update the constructor to remove `partner_block` parameter.

- [ ] **Step 3: Remove old pointer fields from SimUnit**

Remove from SimUnit:
- `SimBlock* first_block_;`
- `SimBlock* last_block_;`
- `SimUnit* next_unit_;`

- [ ] **Step 4: Remove old pointer fields from Simulator**

Remove:
- `SimUnit* first_unit_;`
- `SimUnit* last_unit_;`
- `SimUnit* current_unit_;`
- `SimBlock* current_block_;`

- [ ] **Step 5: Replace VLA thread arrays with vector<jthread>**

At `Simulator.cpp:2975`:

```cpp
// Old:
// thread threads[num_threads];
// for (auto i = num_threads; i--;) {
//     threads[i] = thread(SimulationThread, std::ref(*this), std::ref(ref), std::cref(stats), std::cref(estimates));
// }
// for (auto i = num_threads; i--;) {
//     threads[i].join();
// }

// New:
{
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (auto i = num_threads; i--;) {
        threads.emplace_back([this, &ref, &stats, &estimates](std::stop_token) {
            SimulationThread(*this, ref, stats, estimates);
        });
    }
    // jthread destructors join on scope exit
}
```

Same pattern at `Simulator.cpp:3129` for `ErrorModelOnlyThread`.

Add `#include <thread>` if needed (should already be present for `std::thread`).

- [ ] **Step 6: Build and fix compilation errors**

```bash
make build 2>&1 | head -80
```

Iterate until clean build.

- [ ] **Step 7: Run full test suite**

```bash
make test
```

- [ ] **Step 8: Run sanitizers**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 9: Commit**

```bash
git add reseq/Simulator.h reseq/Simulator.cpp
git commit -m "refactor(4a): complete Simulator container migration and adopt jthread

Remove all pointer fields from SimBlock/SimUnit/Simulator. All access
through deque indices. Replace VLA thread arrays with vector<jthread>
for both SimulationThread and ErrorModelOnlyThread. Workers run to
natural completion (stop_token accepted but not checked)."
```

---

## Task 11: Simulator — Add Unit Tests

**Files:**
- Modify: `reseq/SimulatorTest.cpp` (add new test cases)
- Modify: `reseq/SimulatorTest.h` (add test declarations if needed)

- [ ] **Step 1: Add SimulatorBlockLifecycleTest**

```cpp
TEST_F(SimulatorTest, BlockLifecycle) {
    // Use friend access to directly manipulate blocks_ deque
    // Create blocks, verify indices, push to free list, verify reuse
}
```

- [ ] **Step 2: Add SimulatorPartnerLinkageTest**

```cpp
TEST_F(SimulatorTest, PartnerLinkage) {
    // Create a unit, verify forward-reverse partner_block_idx_ symmetry
    // For each forward block i: blocks_[i].partner_block_idx_ points to reverse block
    // Reverse block's partner_block_idx_ points to next reverse block
}
```

- [ ] **Step 3: Add SimulatorGetNextBlockSequenceTest**

```cpp
TEST_F(SimulatorTest, GetNextBlockSequence) {
    // Call GetNextBlock repeatedly
    // Verify blocks are dispatched in order across units
    // Verify current_block_idx_ and current_unit_idx_ advance correctly
}
```

- [ ] **Step 4: Add SimulatorCleanupWhileSimulatingTest**

```cpp
TEST_F(SimulatorTest, CleanupWhileSimulating) {
    // Create blocks, mark some as finished
    // Trigger cleanup (via CreateBlock which runs cleanup internally)
    // Verify: finished blocks' indices in free_block_indices_
    // Verify: unfinished blocks NOT in free list
}
```

- [ ] **Step 5: Build and run tests**

```bash
make build && make test
```

- [ ] **Step 6: Commit**

```bash
git add reseq/SimulatorTest.cpp reseq/SimulatorTest.h
git commit -m "test(4a): add unit tests for Simulator deque container

Tests cover: block lifecycle (create/free/reuse), partner linkage
symmetry, GetNextBlock dispatch sequence, cleanup-while-simulating."
```

---

## Task 12: FragmentDistributionStats — Add BoundedWorkQueue Class

**Files:**
- Create: `reseq/BoundedWorkQueue.h`

- [ ] **Step 1: Write BoundedWorkQueue unit test first (TDD)**

Add to `reseq/FragmentDistributionStatsTest.cpp` (or a new test file — but using existing is simpler since it already links `reseq_lib`):

```cpp
#include "BoundedWorkQueue.h"

TEST(BoundedWorkQueueTest, TryAcquireWhenFull) {
    reseq::BoundedWorkQueue<10> queue;
    std::vector<size_t> acquired;
    for (int i = 0; i < 10; ++i) {
        size_t idx = queue.try_acquire();
        ASSERT_NE(idx, SIZE_MAX) << "Failed to acquire slot " << i;
        acquired.push_back(idx);
    }
    // Queue is now full
    EXPECT_EQ(queue.try_acquire(), SIZE_MAX);

    // Release one and try again
    queue.release(acquired[0]);
    size_t reacquired = queue.try_acquire();
    EXPECT_NE(reacquired, SIZE_MAX);
    EXPECT_EQ(reacquired, acquired[0]);
}

TEST(BoundedWorkQueueTest, AcquirePublishRelease) {
    reseq::BoundedWorkQueue<10> queue;
    size_t idx = queue.try_acquire();
    ASSERT_NE(idx, SIZE_MAX);

    queue.publish(idx, 5);
    auto& slot = queue.slot(idx);
    EXPECT_EQ(slot.current_param.load(), 0u);
    EXPECT_EQ(slot.total_params, 5u);

    // Simulate completion
    slot.finished_count = 5;
    queue.release(idx);

    // Slot should be reusable
    size_t idx2 = queue.try_acquire();
    EXPECT_EQ(idx2, idx);
}

TEST(BoundedWorkQueueTest, ReleaseEmpty) {
    reseq::BoundedWorkQueue<10> queue;
    size_t idx = queue.try_acquire();
    ASSERT_NE(idx, SIZE_MAX);

    queue.release_empty(idx);

    // Slot should be immediately reusable
    size_t idx2 = queue.try_acquire();
    EXPECT_EQ(idx2, idx);
}

TEST(BoundedWorkQueueTest, ConcurrentAcquireRelease) {
    reseq::BoundedWorkQueue<100> queue;
    constexpr int kThreads = 8;
    constexpr int kOpsPerThread = 1000;
    std::atomic<int> total_acquired{0};

    std::vector<std::jthread> threads;
    for (int t = 0; t < kThreads; ++t) {
        threads.emplace_back([&queue, &total_acquired](std::stop_token) {
            for (int i = 0; i < kOpsPerThread; ++i) {
                size_t idx = queue.try_acquire();
                if (idx != SIZE_MAX) {
                    ++total_acquired;
                    queue.release(idx);
                }
            }
        });
    }
    threads.clear(); // join all

    // All slots should be back in the free list
    int available = 0;
    for (int i = 0; i < 100; ++i) {
        if (queue.try_acquire() != SIZE_MAX) ++available;
    }
    EXPECT_EQ(available, 100);
}
```

- [ ] **Step 2: Run tests to verify they fail (class doesn't exist yet)**

```bash
make build 2>&1 | head -20
```

Expected: compilation failure — `BoundedWorkQueue.h` not found.

- [ ] **Step 3: Create BoundedWorkQueue.h**

Create `reseq/BoundedWorkQueue.h`:

```cpp
#ifndef BOUNDEDWORKQUEUE_H
#define BOUNDEDWORKQUEUE_H

#include <array>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <vector>

namespace reseq {

struct QueueSlot {
    std::atomic<uint32_t> current_param{0};
    std::atomic<uint32_t> finished_count{0};
    uint32_t total_params{0};
    bool published{false};
};

template <size_t MaxSlots>
class BoundedWorkQueue {
    std::array<QueueSlot, MaxSlots> slots_;
    std::vector<size_t> free_indices_;
    mutable std::mutex mutex_;

  public:
    BoundedWorkQueue() {
        free_indices_.reserve(MaxSlots);
        for (size_t i = MaxSlots; i--;) {
            free_indices_.push_back(i);
        }
    }

    // NON-BLOCKING: returns slot index, or SIZE_MAX if full
    size_t try_acquire() {
        std::lock_guard lock(mutex_);
        if (free_indices_.empty()) return SIZE_MAX;
        size_t idx = free_indices_.back();
        free_indices_.pop_back();
        return idx;
    }

    // Called by producer after populating slot parameters.
    // Sets current_param to 0 LAST (publication ordering).
    void publish(size_t idx, uint32_t total) {
        slots_[idx].finished_count = 0;
        slots_[idx].total_params = total;
        slots_[idx].published = true;
        slots_[idx].current_param = 0; // LAST — publication signal
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
    const QueueSlot& slot(size_t idx) const { return slots_[idx]; }
    static constexpr size_t max_slots() { return MaxSlots; }
};

} // namespace reseq

#endif // BOUNDEDWORKQUEUE_H
```

- [ ] **Step 4: Build and run tests**

```bash
make build && make test
```

Expected: all tests pass including new BoundedWorkQueue tests.

- [ ] **Step 5: Commit**

```bash
git add reseq/BoundedWorkQueue.h reseq/FragmentDistributionStatsTest.cpp
git commit -m "feat(4a): add BoundedWorkQueue with non-blocking try_acquire

Header-only bounded work queue for FragmentDistributionStats bias
calculation pipeline. Preserves non-blocking guarantee: try_acquire()
returns SIZE_MAX when full. Publication ordering: current_param set
last. Includes unit tests for full queue, publish/release, concurrent
acquire/release."
```

---

## Task 13: FragmentDistributionStats — Integrate BoundedWorkQueue

**Files:**
- Modify: `reseq/FragmentDistributionStats.h:391-393` (replace atomic arrays)
- Modify: `reseq/FragmentDistributionStats.cpp:2957-2985` (AddNewBiasCalculations)
- Modify: `reseq/FragmentDistributionStats.cpp:2987-3036` (ExecuteBiasCalculations)
- Modify: `reseq/FragmentDistributionStats.cpp:2169-2248` (UpdateBiasCalculationParams)

- [ ] **Step 1: Add BoundedWorkQueue member alongside old arrays**

In `reseq/FragmentDistributionStats.h`, add after the old atomic arrays (line 393):

```cpp
#include "BoundedWorkQueue.h"
// ...
BoundedWorkQueue<kMaxBinsQueuedForBiasCalc> bias_queue_;
```

Initialize in `Prepare` (after existing queue initialization code):

```cpp
// BoundedWorkQueue self-initializes in constructor — no additional setup needed
```

- [ ] **Step 2: Convert AddNewBiasCalculations to use BoundedWorkQueue**

Replace the `test_and_set` scan loop in `AddNewBiasCalculations` (`FragmentDistributionStats.cpp:2962-2982`):

```cpp
// Old:
// uint32_t queue_spot = 0;
// while (queue_spot < kMaxBinsQueuedForBiasCalc && claimed_bias_bins_.at(queue_spot).test_and_set()) {
//     ++queue_spot;
// }
// if (queue_spot < kMaxBinsQueuedForBiasCalc) { ... } else { break; }

// New:
size_t queue_spot = bias_queue_.try_acquire();
if (queue_spot != SIZE_MAX) {
    // CAS on ref_seq_bin as before
    if (num_handled_reference_sequence_bins_.compare_exchange_strong(ref_seq_bin, ref_seq_bin + 1)) {
        UpdateBiasCalculationParams(ref_seq_bin, queue_spot, thread.bias_calc_tmp_params_, print_mutex);
        ref_seq_bin = num_handled_reference_sequence_bins_;
    } else {
        bias_queue_.release_empty(queue_spot);
    }
} else {
    break; // Queue full — non-blocking guarantee preserved
}
```

- [ ] **Step 3: Convert UpdateBiasCalculationParams to use BoundedWorkQueue**

At the end of `UpdateBiasCalculationParams` (lines 2239-2246), replace:

```cpp
// Old:
// finished_bias_calcs_.at(queue_spot) = 0;
// current_bias_param_.at(queue_spot) = 0;

// New:
auto queue_bin_size = bias_calc_params_.size() / kMaxBinsQueuedForBiasCalc;
bias_queue_.publish(queue_spot, queue_bin_size);
```

For the empty-bin early return:

```cpp
// Old:
// params_left_for_calculation_ -= qbin_size;
// claimed_bias_bins_.at(queue_spot).clear();

// New:
params_left_for_calculation_ -= qbin_size;
bias_queue_.release_empty(queue_spot);
```

- [ ] **Step 4: Convert ExecuteBiasCalculations to use BoundedWorkQueue**

Replace slot access in `ExecuteBiasCalculations` (lines 2987-3036):

```cpp
// Work stealing per slot:
// Old: auto cur_par = queue_bin * queue_bin_size + current_bias_param_.at(queue_bin)++;
// New: auto cur_par = queue_bin * queue_bin_size + bias_queue_.slot(queue_bin).current_param++;

// Completion check:
// Old: if (++finished_bias_calcs_.at(queue_bin) == queue_bin_size) { claimed_bias_bins_.at(queue_bin).clear(); }
// New: if (++bias_queue_.slot(queue_bin).finished_count == queue_bin_size) { bias_queue_.release(queue_bin); }
```

- [ ] **Step 5: Remove old atomic arrays**

Remove from `FragmentDistributionStats.h`:
- `std::array<std::atomic_flag, kMaxBinsQueuedForBiasCalc> claimed_bias_bins_;`
- `std::array<std::atomic<uintNumFits>, kMaxBinsQueuedForBiasCalc> current_bias_param_;`
- `std::array<std::atomic<uintNumFits>, kMaxBinsQueuedForBiasCalc> finished_bias_calcs_;`

Remove their initialization in `Prepare`.

- [ ] **Step 6: Build and run full test suite**

```bash
make build && make test
```

- [ ] **Step 7: Commit**

```bash
git add reseq/FragmentDistributionStats.h reseq/FragmentDistributionStats.cpp
git commit -m "refactor(4a): integrate BoundedWorkQueue into FragmentDistributionStats

Replace claimed_bias_bins_/current_bias_param_/finished_bias_calcs_
atomic arrays with BoundedWorkQueue. Non-blocking guarantee preserved:
try_acquire() returns SIZE_MAX when full. Publication ordering preserved:
current_param set last via publish()."
```

---

## Task 14: FragmentDistributionStats — Convert Vector Pools and Adopt jthread

**Files:**
- Modify: `reseq/FragmentDistributionStats.h:401-404` (vector pool declarations)
- Modify: `reseq/FragmentDistributionStats.cpp:3001-3018` (spin-wait removal)
- Modify: `reseq/FragmentDistributionStats.cpp:2767-2775` (BiasSumThread jthread)
- Modify: `reseq/FragmentDistributionStats.cpp:3723-3731` (BiasNormalizationThread jthread)

- [ ] **Step 1: Convert bias_calc_vects_ to thread-indexed vector**

In `reseq/FragmentDistributionStats.h`, replace:

```cpp
// Old:
// std::deque<std::pair<std::atomic_flag, BiasCalculationVectors>> bias_calc_vects_;
// std::deque<std::pair<std::atomic_flag, std::vector<uintFragCount>>> tmp_frag_count_;

// New:
std::vector<BiasCalculationVectors> bias_calc_vects_;
std::vector<std::vector<uintFragCount>> tmp_frag_count_;
```

In `Prepare`, replace initialization:

```cpp
// Old:
// bias_calc_vects_.resize(num_threads);
// tmp_frag_count_.resize(num_threads);
// for (auto i = tmp_frag_count_.size(); i--;) {
//     bias_calc_vects_.at(i).first.clear();
//     tmp_frag_count_.at(i).first.clear();
// }

// New:
bias_calc_vects_.resize(num_threads);
tmp_frag_count_.resize(num_threads);
// No atomic flags to clear — thread-index access has no contention
```

- [ ] **Step 2: Pass thread index through to ExecuteBiasCalculations**

Add a `size_t thread_idx` parameter to `ExecuteBiasCalculations` and `HandleReferenceSequencesUntil`. Thread index is available from the DataStats `ReadThread` lambda (Tranche 1 already uses indexed jthread). Update call chain to pass it through.

In `ExecuteBiasCalculations`, replace the spin-wait:

```cpp
// Old:
// uintNumFits nvar = 0;
// while (bias_calc_vects_.at(nvar).first.test_and_set()) { ++nvar; }
// ... use bias_calc_vects_.at(nvar).second ...
// bias_calc_vects_.at(nvar).first.clear();

// New:
auto& calc_vect = bias_calc_vects_[thread_idx];
// ... use calc_vect directly ...
// No release needed — exclusive by construction
```

Same pattern for `tmp_frag_count_`.

- [ ] **Step 3: Replace VLA thread arrays with vector<jthread>**

At `FragmentDistributionStats.cpp:2767` (BiasSumThread):

```cpp
// Old:
// thread threads[num_threads];
// for (auto i = num_threads; i--;) {
//     threads[i] = thread(BiasSumThread, ...);
// }
// for (auto i = num_threads; i--;) { threads[i].join(); }

// New:
{
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (auto i = num_threads; i--;) {
        threads.emplace_back([this, &reference, &params, &current_param,
                              &finished_params, &bias_sum, &print_mutex](std::stop_token) {
            BiasSumThread(*this, reference, params, current_param,
                          finished_params, bias_sum, print_mutex);
        });
    }
    // jthread destructors join on scope exit
}
```

Same pattern at line 3723 for BiasNormalizationThread.

- [ ] **Step 4: Build and run full test suite**

```bash
make build && make test
```

- [ ] **Step 5: Run sanitizers**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 6: Commit**

```bash
git add reseq/FragmentDistributionStats.h reseq/FragmentDistributionStats.cpp
git commit -m "refactor(4a): convert FragmentDistributionStats to thread-indexed pools and jthread

Replace spin-wait bias_calc_vects_/tmp_frag_count_ with thread-indexed
vectors (no contention). Replace VLA thread arrays for BiasSumThread
and BiasNormalizationThread with vector<jthread>. Workers run to
natural completion (stop_token accepted but not checked)."
```

---

## Task 15: FragmentDistributionStats — Add Integration Tests

**Files:**
- Modify: `reseq/FragmentDistributionStatsTest.cpp`

- [ ] **Step 1: Add NonBlockingProgressGuarantee test**

```cpp
TEST_F(FragmentDistributionStatsTest, NonBlockingProgressGuarantee) {
    // Fill the BoundedWorkQueue completely (acquire all 100 slots)
    // From a worker thread, call AddNewBiasCalculations
    // Verify: it returns without blocking (does not hang)
    // Verify: ExecuteBiasCalculations drains work, freeing slots
    // Use a timeout to detect deadlock: std::future with wait_for
}
```

- [ ] **Step 2: Add TerminationGuarantee test**

```cpp
TEST_F(FragmentDistributionStatsTest, TerminationGuarantee) {
    // Run a small bias calculation pipeline end-to-end
    // Verify: params_left_for_calculation_ reaches 0
    // Verify: FinishThreads exits cleanly
}
```

- [ ] **Step 3: Build and run tests**

```bash
make build && make test
```

- [ ] **Step 4: Commit**

```bash
git add reseq/FragmentDistributionStatsTest.cpp
git commit -m "test(4a): add integration tests for FragmentDistributionStats concurrency

Tests cover: non-blocking progress guarantee (full queue doesn't deadlock),
termination guarantee (params_left_for_calculation_ reaches 0)."
```

---

## Task 16: Full Sanitizer Verification

**Files:** None modified — verification only.

- [ ] **Step 1: Run TSan on full test suite**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

Document any TSan warnings. If pre-existing races appear (from code outside scope), note them but do not fix in this task.

- [ ] **Step 2: Run ASan+UBSan on full test suite**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 3: Run normal build + test (verify non-sanitizer build still works)**

```bash
rm -rf build && make build && make test
```

All 65 unit tests + 10 regression tests + new tests must pass.

- [ ] **Step 4: Format check**

```bash
make format-check
```

If formatting drift, run `make format` and commit.

- [ ] **Step 5: Commit any formatting fixes**

```bash
# Only if Step 4 found issues:
make format
git add -A
git commit -m "style: fix clang-format violations in Phase 4 Tranche 2 files"
```

---

## Verification Protocol

**After every task:**

```bash
make build && make test
```

All 65 existing unit tests + 10 regression tests + new tests must pass. Do NOT use `--gtest_filter` as a reliable checkpoint — tests have inter-test dependencies (shared `Register()` state).

**After Tasks 5, 10, 14 (major milestones):**

Run ASan+UBSan:
```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

**After Task 16 (final):**

Run TSan, ASan+UBSan, and normal build. All must pass clean.
