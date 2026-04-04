# Simulator Concurrency Modernization — Design Document

**Date:** 2026-04-03
**Phase:** 4a (Tranche 2, Task 7)
**Status:** Draft
**Prerequisite:** Phase 3 (PR #5, `c914db2`) and Tranche 1 (PR #6, `296c273`) merged to master.

---

## 1. Current Architecture

Simulator uses two interlinked singly-linked lists: **SimUnit** chains (one per reference
sequence) and **SimBlock** chains (genomic intervals within a unit). Each forward block
has a `partner_block_` pointing to its reverse-strand counterpart. `next_block_` is
`atomic<SimBlock*>` for lock-free traversal; all other pointers are raw.

Two VLA thread arrays remain unconverted:
- `SimulationThread` at `Simulator.cpp:2975` — full simulation pipeline
- `ErrorModelOnlyThread` at `Simulator.cpp:3129` — error model application only

Block creation and cleanup are serialized by `block_creation_mutex_` in `GetNextBlock`.
Output uses a 3-mutex cascade (`output_mutex_` + `flush_mutex_[2]`) to guarantee
paired-end write ordering.

### Key Data Structures

```
Simulator
  first_unit_    : SimUnit*                  — front of unit chain
  last_unit_     : SimUnit*                  — end of unit chain
  current_unit_  : SimUnit*                  — next unit to dispatch
  current_block_ : SimBlock*                 — next block to dispatch
  deletion_buffer_     : uintRefSeqBin       — current buffer size
  req_deletion_buffer_ : atomic<uintRefSeqBin>  — requested buffer size
  written_records_     : atomic<uintFragCount>

SimUnit
  ref_seq_id_    : uintRefSeqId (const)
  first_block_   : SimBlock*                 — owning, front of block chain
  last_block_    : SimBlock*                 — non-owning, end of block chain
  next_unit_     : SimUnit*                  — owning, next unit in chain

SimBlock
  id_            : uintRefSeqBin (const)
  start_pos_     : uintSeqLen (const)
  finished_      : atomic<bool>              — simulation complete flag
  next_block_    : atomic<SimBlock*>         — owning, forward traversal
  partner_block_ : SimBlock*                 — non-owning, reverse strand cross-ref
  seed_          : uintSeed
  sys_errors_    : vector<pair<Dna5, uintPercent>>
  err_variants_  : vector<SysErrorVariant>
  first_variant_id_     : intVariantId
  first_methylation_id_ : intVariantId
```

### Ownership Model (Phase 3 Documentation)

- SimUnit owns its chain of SimBlocks (`first_block_` through `next_block_` links).
- Simulator owns its chain of SimUnits (`first_unit_` through `next_unit_` links).
- `partner_block_` is a non-owning cross-reference between forward/reverse block chains.
- Forward block owns its reverse partner: `delete del_block->partner_block_` before
  `delete del_block` (`Simulator.cpp:1327-1328`).

### Mutexes

| Mutex | Purpose | Code Site |
|-------|---------|-----------|
| `print_mutex_` | Console output | `Simulator.h:263` |
| `output_mutex_` | Output buffer + Flush entry | `Simulator.h:264` |
| `flush_mutex_[2]` | Paired-end write ordering cascade | `Simulator.h:265` |
| `block_creation_mutex_` | Serializes block creation and cleanup | `Simulator.h:267` |
| `var_read_mutex_` | Variant file reading | `Simulator.h:268` |
| `methylation_read_mutex_` | Methylation file reading | `Simulator.h:269` |

### Thread Model

Two thread creation sites, both VLA arrays with manual join:

**SimulationThread** (`Simulator.cpp:2522-2541`): Main loop calls `GetNextBlock()`,
simulates block, sets `finished_ = true`. After all blocks processed, one thread
runs adapter-only simulation (guarded by `adapter_only_simulated_` atomic_flag).

**ErrorModelOnlyThread** (`Simulator.cpp:2661-2709`): Reads batches of sequences,
applies error model. Different entry point, same output pipeline.

---

## 2. Invariant Catalog

| # | Invariant | Current Mechanism | Code Site | Safety Status |
|---|-----------|-------------------|-----------|---------------|
| S1 | Block creation and cleanup are mutually exclusive | `block_creation_mutex_` in `GetNextBlock` | `Simulator.cpp:1368` | Safe |
| S2 | Forward block owns its reverse partner | `delete partner_block_` before `delete del_block` | `Simulator.cpp:1327-1328` | Safe but manual |
| S3 | Lock-free forward traversal during simulation | `atomic<SimBlock*> next_block_` | `Simulator.h:122` | Safe |
| S4 | Block completion is signaled atomically | `atomic<bool> finished_` written by worker, read by cleanup | `Simulator.cpp:2533, 1311` | Safe |
| S5 | Paired-end segments written in synchronized order | Cascade: lock flush[0], write seg0, lock flush[1], unlock flush[0], write seg1, unlock flush[1] | `Simulator.cpp:169-177` | Safe |
| S6 | Output buffer swap is atomic with respect to appends | `output_mutex_` held during append and Flush entry | `Simulator.cpp:227, 167` | Safe |
| S7 | Variant/methylation loading is opportunistic | `try_lock` on `var_read_mutex_`/`methylation_read_mutex_` | `Simulator.cpp:1346, 1356` | Safe |
| S8 | Deletion buffer prevents premature cleanup | `req_deletion_buffer_` (atomic) requests extra blocks; `CheckDeletionBuffer` creates them under `block_creation_mutex_` | `Simulator.h:280-282` |  Safe |
| S9 | Adapter-only simulation runs exactly once | `atomic_flag adapter_only_simulated_` with `test_and_set` | `Simulator.cpp:2538` | Safe |
| S10 | Reverse blocks created before forward blocks | `CreateUnit` builds all reverse blocks first; forward derives partner via double-hop | `Simulator.cpp:932-946, 1270` | Safe |
| S11 | Current block/unit advance under mutex | `current_block_`/`current_unit_` updated under `block_creation_mutex_` | `Simulator.cpp:1380-1393` | Safe |

### Invariant Details

**S1 — Block creation/cleanup serialization:** `GetNextBlock` (`Simulator.cpp:1368`)
acquires `block_creation_mutex_` before calling `CreateBlock` and `CheckDeletionBuffer`.
Cleanup of finished blocks also happens inside `CreateBlock` (lines 1308-1330). This
serialization prevents concurrent block creation and deletion, which would corrupt the
linked list. The mutex also protects `current_block_`/`current_unit_` advancement.

**S2 — Forward owns reverse:** During cleanup (`Simulator.cpp:1327-1328`), the forward
block's `partner_block_` (reverse strand) is deleted first, then the forward block itself.
This is the ownership convention: forward blocks are the primary chain, reverse blocks
are satellites. `CreateUnit` builds reverse blocks first, then forward blocks, linking
each forward block's `partner_block_` to its corresponding reverse block.

**S5 — Flush cascade:** The 3-mutex cascade in `Flush` (`Simulator.cpp:154-191`):
1. `output_mutex_` locked by caller (`Output`, line 227)
2. Swap buffers under `output_mutex_`, unlock (`line 167`)
3. Lock `flush_mutex_[0]`, write segment 0 (`lines 169-171`)
4. Lock `flush_mutex_[1]` BEFORE unlocking `flush_mutex_[0]` (`line 172`)
5. Unlock `flush_mutex_[0]` (`line 174`)
6. Write segment 1, unlock `flush_mutex_[1]` (`lines 175-177`)

This hand-over-hand locking ensures paired-end reads appear in the same order in both
output files. The lock acquisition order is always 0->1, preventing deadlock.
This pattern is correct and should be preserved as-is.

**S8 — Deletion buffer:** When systematic errors include long deletions that span block
boundaries, `RequestBufferSize` (`Simulator.cpp:460-467`) atomically requests extra
buffer blocks. `CheckDeletionBuffer` (`Simulator.cpp:468-474`) creates them under
`block_creation_mutex_`. This prevents a block from being deleted while a deletion
from a previous block still references positions in it.

**S10 — Reverse-first creation:** `CreateUnit` (`Simulator.cpp:932-1041`) creates ALL
reverse blocks first, linking them via `next_block_` (in reverse direction) and
`partner_block_` (previous reverse block). Then `CreateBlock` creates forward blocks,
deriving the partner via double-hop: `last_block_->partner_block_->partner_block_`
(`Simulator.cpp:1270`). This double-hop works because:
- `last_block_->partner_block_` = reverse block corresponding to the last forward block
- That reverse block's `partner_block_` = the NEXT reverse block (set during
  `CreateUnit`, line 944) = the partner for the new forward block

---

## 3. Target Architecture

Replace intrusive linked lists with **indexed owning containers** using `std::deque`
for pointer stability. Replace VLA thread arrays with `std::jthread`.

### Why Deque

Same rationale as CoverageStats: front consumption (cleanup) + back appending with
pointer/reference stability. Worker threads receive `SimBlock&` references from
`GetNextBlock` and use them throughout simulation — those references must remain valid
even as new blocks are appended at the back.

### Container Changes

| Current | Target | Notes |
|---------|--------|-------|
| `SimBlock* first_block_` (in SimUnit) | `size_t first_block_idx_` | Index into blocks deque |
| `SimBlock* last_block_` (in SimUnit) | `size_t last_block_idx_` | Index into blocks deque |
| `atomic<SimBlock*> next_block_` | `size_t next_block_idx_` | `SIZE_MAX` = null, immutable once set |
| `SimBlock* partner_block_` | `size_t partner_block_idx_` | `SIZE_MAX` = null |
| `SimUnit* first_unit_` | `size_t first_unit_idx_` | Index into units deque |
| `SimUnit* last_unit_` | `size_t last_unit_idx_` | Index into units deque |
| `SimUnit* next_unit_` | `size_t next_unit_idx_` | `SIZE_MAX` = null |
| `SimUnit* current_unit_` | `size_t current_unit_idx_` | Dispatch pointer |
| `SimBlock* current_block_` | `size_t current_block_idx_` | Dispatch pointer |

### New Simulator Members

```cpp
class Simulator {
    std::deque<SimBlock> blocks_;              // Owns all blocks (value semantics)
    std::deque<SimUnit> units_;                // Owns all units (value semantics)
    std::vector<size_t> free_block_indices_;   // Recycled block slots
    std::vector<size_t> free_unit_indices_;    // Recycled unit slots
    size_t first_unit_idx_;
    size_t last_unit_idx_;
    size_t current_unit_idx_;
    size_t current_block_idx_;
};
```

**Value semantics for SimBlock/SimUnit:** Unlike CoverageBlock (which has atomic members
preventing move), SimBlock and SimUnit can be stored by value in the deque. The deque's
pointer stability guarantee means `&blocks_[i]` remains valid across push_back.
SimBlock's `finished_` (atomic<bool>) is not movable, so blocks are constructed in-place
via `emplace_back`.

### SimBlock Changes

```cpp
struct SimBlock {
    const uintRefSeqBin id_;
    const uintSeqLen start_pos_;
    std::atomic<bool> finished_;
    size_t next_block_idx_;          // Replaces atomic<SimBlock*> next_block_
    size_t partner_block_idx_;       // Replaces SimBlock* partner_block_
    uintSeed seed_;
    std::vector<std::pair<seqan::Dna5, uintPercent>> sys_errors_;
    std::vector<SysErrorVariant> err_variants_;
    intVariantId first_variant_id_;
    intVariantId first_methylation_id_;
};
```

### SimUnit Changes

```cpp
struct SimUnit {
    const uintRefSeqId ref_seq_id_;
    size_t first_block_idx_;         // Replaces SimBlock* first_block_
    size_t last_block_idx_;          // Replaces SimBlock* last_block_
    size_t next_unit_idx_;           // Replaces SimUnit* next_unit_
};
```

### Thread Changes

Replace both VLA `thread[]` arrays with `vector<jthread>`:

```cpp
// SimulationThread (Simulator.cpp:2975)
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

// ErrorModelOnlyThread (Simulator.cpp:3129) — same pattern
```

Workers accept `stop_token` but do NOT check it — they run to natural completion,
preserving identical semantics to the original code. Same pattern as Tranche 1's
ProbabilityEstimates and DataStats conversions.

### Invariant Preservation

**S1 (mutual exclusion):** `block_creation_mutex_` unchanged — still serializes
`GetNextBlock`. The mutex now protects index-based operations instead of pointer-based.

**S2 (forward owns reverse):** Explicit via indexed containers. Cleanup marks both the
forward block's index and its `partner_block_idx_` as free (push to `free_block_indices_`)
instead of `delete`. Both blocks remain in the deque until their slots are reused.
Cleanup frees forward and partner indices together in the correct order.

**S3 (lock-free traversal):** `next_block_idx_` is set before publication (before
`current_block_idx_` advances under `block_creation_mutex_`). The index is immutable
once set. Workers look up `blocks_[next_block_idx_]` — deque guarantees reference
stability. No atomic needed.

**S4 (completion signal):** `atomic<bool> finished_` stays — it's a block field, not
a container concern.

**S5 (Flush cascade):** **Unchanged.** The 3-mutex cascade is a correct hand-over-hand
locking pattern using standard `std::mutex`. It cannot deadlock (acquisition order is
always 0->1). Simplifying to a single `flush_mutex_` would serialize both segments
unnecessarily, halving write throughput. Keep as-is with added documentation comment.

**S6 (output buffer):** Unchanged — `output_mutex_` pattern stays.

**S7 (opportunistic loading):** Unchanged — try_lock pattern stays.

**S8 (deletion buffer):** `deletion_buffer_` / `req_deletion_buffer_` unchanged
conceptually. `CheckDeletionBuffer` creates blocks via deque emplacement instead of `new`.

**S9 (adapter-only):** Unchanged — atomic_flag stays.

**S10 (reverse-first creation):** `CreateUnit` still creates all reverse blocks first,
emplaced into `blocks_` deque. Partner linkage by stored index instead of double-hop
pointer chase: when creating a forward block, the partner index is computed directly
from the reverse block layout (reverse blocks are contiguous in the deque from
`CreateUnit`).

**S11 (current pointer advance):** Index-based: `current_block_idx_` and
`current_unit_idx_` advance under `block_creation_mutex_`.

### Flush Analysis

The cascade lock pattern in `Flush` (`Simulator.cpp:154-191`) is analyzed here for
completeness:

```
Thread A (Flush #1):                    Thread B (Flush #2):
  lock output_mutex_                      (blocked on output_mutex_)
  swap buffers
  unlock output_mutex_
  lock flush_mutex_[0]                    lock output_mutex_
  write segment 0                         swap buffers
  lock flush_mutex_[1]                    unlock output_mutex_
  unlock flush_mutex_[0]                  lock flush_mutex_[0]  ← blocks until A releases
  write segment 1                         ...
  unlock flush_mutex_[1]
```

The cascade ensures:
1. Segment 0 of flush N completes before segment 0 of flush N+1
2. Segment 1 of flush N completes before segment 1 of flush N+1
3. Within a single flush, segment 0 is written before segment 1

This guarantees paired-end reads appear in the same order in both output files.
No redesign needed.

---

## 4. Migration Plan

Each step compiles and passes the full test suite.

1. **Add deque containers and index fields.** Add `blocks_` deque, `units_` deque,
   `free_block_indices_`, `free_unit_indices_` to Simulator. Add index fields to
   SimBlock and SimUnit alongside existing pointer fields. Initialize in `Simulate`
   and `ApplyErrorModel`. No behavioral change.

2. **Convert CreateUnit to use deques.** Emplace reverse blocks into `blocks_`,
   emplace unit into `units_`. Set index fields. Keep pointer fields as aliases
   derived from `&blocks_[idx]`.

3. **Convert CreateBlock to use deque + free list.** New blocks emplaced at back or
   reused from free list. Both index and pointer fields set.

4. **Convert cleanup in CreateBlock.** Replace `delete partner_block_; delete del_block;`
   with push of both indices to `free_block_indices_`. Replace `delete del_unit;` with
   push to `free_unit_indices_`. Clear block state for reuse.

5. **Convert GetNextBlock to index-based dispatch.** `current_block_idx_` and
   `current_unit_idx_` advance. Return block reference via `blocks_[idx]`.

6. **Convert SimulationThread and ErrorModelOnlyThread.** Receive block/unit by
   reference (looked up from deque by index in GetNextBlock). Internal simulation
   code operates on `SimBlock&` — no changes needed inside simulation logic.

7. **Convert final cleanup** (`Simulator.cpp:2995-3018`). Walk via indices. Clear
   deques at end.

8. **Remove pointer fields.** Delete `next_block_` (atomic), `partner_block_`,
   `first_block_`, `last_block_`, `next_unit_`, `first_unit_`, `last_unit_`,
   `current_unit_`, `current_block_` pointer fields. All access is now index-based.

9. **Replace VLA thread arrays with vector<jthread>.** Convert both sites
   (lines 2975 and 3129) to `vector<jthread>` with natural completion semantics.

---

## 5. Testing

### New Unit Tests

| Test | What It Validates |
|------|-------------------|
| `SimulatorBlockLifecycleTest` | Create blocks via deque, verify index assignment, free blocks, verify reuse from free list |
| `SimulatorPartnerLinkageTest` | Create a unit with CreateUnit, verify forward-reverse partner indices are symmetric (forward's partner = reverse block, reverse's partner = next reverse = next forward's partner) |
| `SimulatorDeletionBufferTest` | Request buffer sizes via `RequestBufferSize`, verify `CheckDeletionBuffer` creates extra blocks in deque |
| `SimulatorFlushOrderingTest` | Concurrent Flush calls from multiple threads, verify paired-end output order preserved using mock output that captures write sequence |
| `SimulatorJthreadShutdownTest` | Verify clean shutdown: jthread vector goes out of scope, all threads join, no hanging, no use-after-free |
| `SimulatorGetNextBlockSequenceTest` | Call GetNextBlock N times, verify returned block indices advance correctly through units and blocks |
| `SimulatorCleanupWhileSimulatingTest` | Mark some blocks as finished, call CreateBlock (which triggers cleanup), verify only finished blocks are freed and in-flight blocks are untouched |

### Existing Tests

All tests in `SimulatorTest` suite plus the 10 golden-file regression tests must pass
unchanged at every migration step.

---

## 6. Verification Criteria

- `make build && make test` — all 65 unit tests + 10 regression tests pass
- TSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=thread"` — no data race warnings
- ASan+UBSan clean: `cmake -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined"` — no errors
- No performance regression on ecoli test data (simulation throughput unchanged)
- Flush ordering verified by new test under TSan
