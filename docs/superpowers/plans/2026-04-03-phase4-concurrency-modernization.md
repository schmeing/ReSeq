# Phase 4: Concurrency Modernization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize thread management in ProbabilityEstimates and DataStats to use `std::jthread`, make shared global state atomic, add TSan CI, and fix tooling gaps. Container redesigns for Simulator, CoverageStats, and FragmentDistributionStats are **deferred to Tranche 2** pending design documents that capture ownership, publication, and reclamation invariants.

**Architecture:** Phase 4 is structured in two tranches. **Tranche 1** (Tasks 1-5) addresses safe, mechanical wins: atomic globals, jthread adoption where threading is simple, TSan CI, tooling fixes. **Tranche 2** (Tasks 6-8) produces design documents for the complex subsystems before any code changes.

**Tech Stack:** C++20 (`std::jthread`, `std::stop_token`, `std::atomic`, `std::scoped_lock`), GoogleTest, CMake, TSan, clang-tidy.

**Prerequisite:** Phase 3 (PR #5, merged as `c914db2`) is on master.

**Baseline (verified 2026-04-03):** The repo has C++20, target-based CMake, `reseq_test` binary with CTest, GitHub Actions CI (GCC 13, Clang 17, ASan+UBSan, coverage, format-check). `make build && make test` passes (65 tests). Tests must run as a single invocation — shared `Register()` state requires running together, not individually (`reseq/CMakeLists.txt:60`).

**Scope note:** Manual thread arrays also exist in `reseq/DataStats.cpp:1304`, `reseq/FragmentDistributionStats.cpp:2764`, and `reseq/Simulator.cpp:2974`. DataStats is addressed in Task 3. FragmentDistributionStats and Simulator thread modernization is deferred to Tranche 2 because their threading is coupled to data structures with non-trivial ownership and publication semantics.

---

## Tranche 1: Safe Mechanical Wins

### Task 1: Make kVerbosityLevel Atomic (Phase 4c)

**Files:**
- Modify: `reseq/logging.hpp:23` (extern declaration)
- Modify: `reseq/main.cpp:19-23` (definition + include order)
- Modify: `reseq/test_main.cpp:6-10` (definition + include order)
- Modify: `reseq/BasicTestClass.hpp:24-29,95` (read/write sites)

Non-atomic `uint16_t` read by multiple threads is a data race per the C++ memory model.

**Include order constraint:** Both `main.cpp:20` and `test_main.cpp:7` define `kVerbosityLevel` BEFORE including `logging.hpp`. When changing the type to `std::atomic<uint16_t>`, those definition sites need `<atomic>` included before the definition — not just in `logging.hpp`.

- [ ] **Step 1: Add `<atomic>` to main.cpp and test_main.cpp before the definition**

In `reseq/main.cpp`, add `#include <atomic>` in the include block (before line 18). Then change line 20:
```cpp
// Old:
uint16_t kVerbosityLevel = 99;
// New:
std::atomic<uint16_t> kVerbosityLevel{99};
```

In `reseq/test_main.cpp`, add `#include <atomic>` before line 6. Then change line 7:
```cpp
// Old:
uint16_t kVerbosityLevel = 2;
// New:
std::atomic<uint16_t> kVerbosityLevel{2};
```

- [ ] **Step 2: Change the extern declaration in logging.hpp**

In `reseq/logging.hpp:23`, find:
```cpp
extern uint16_t kVerbosityLevel;
```
Replace with:
```cpp
extern std::atomic<uint16_t> kVerbosityLevel;
```
Ensure `#include <atomic>` is present in `logging.hpp`.

- [ ] **Step 3: Fix boost::program_options interaction in main.cpp**

At `reseq/main.cpp:413`, `value<uint16_t>(&reseq::kVerbosityLevel)` takes a `uint16_t*` which is not layout-compatible with `std::atomic<uint16_t>`. Fix by using a temporary:

Add a local `uint16_t verbosity_opt = 4;` before the options parsing block, use `&verbosity_opt` in the option definition. After parsing, assign: `kVerbosityLevel.store(verbosity_opt);`.

- [ ] **Step 4: Verify all read sites compile**

`std::atomic<uint16_t>` supports implicit conversion to `uint16_t`, so the many `if (0 < kVerbosityLevel)` comparisons in `main.cpp` and the logging macros in `logging.hpp` (`kVerbosityLevel` appears in `printErr`, `printWarn`, `printInfo`, `printSucc`, `printDebug` macro guards) compile without changes.

In `reseq/BasicTestClass.hpp`, fix `std::min` at line 28:
```cpp
// Old:
inline void RestoreTestVerbosity() { kVerbosityLevel = std::min((uint16_t)kTestVerbosity, real_verbosity_); }
// New:
inline void RestoreTestVerbosity() { kVerbosityLevel.store(std::min(kTestVerbosity, real_verbosity_)); }
```

The constructor at line 95 (`real_verbosity_(kVerbosityLevel)`) works via implicit conversion — `kVerbosityLevel.load()` returns `uint16_t`.

- [ ] **Step 5: Build and run full test suite**

```bash
make build && make test
```
Expected: all 65 tests pass.

- [ ] **Step 6: Commit**

```bash
git add reseq/logging.hpp reseq/main.cpp reseq/test_main.cpp reseq/BasicTestClass.hpp
git commit -m "refactor(4c): make kVerbosityLevel std::atomic<uint16_t>

Non-atomic uint16_t read by multiple threads is a data race per C++
memory model. std::atomic makes the contract explicit at zero cost.
Include order preserved: <atomic> added before definition sites in
main.cpp and test_main.cpp which define kVerbosityLevel before
including logging.hpp."
```

---

### Task 2: ProbabilityEstimates — jthread Adoption (Phase 4b)

**Files:**
- Modify: `reseq/ProbabilityEstimates.h` (thread declaration, IPFThread signature)
- Modify: `reseq/ProbabilityEstimates.cpp:988,1195-1205` (IPFThread, thread creation/join)

ProbabilityEstimates has a simple threading model: atomic fetch-and-increment work distribution. Each thread pulls the next parameter index via `current_param_++` and processes it. Threads terminate when `cur_par >= params.size()` or `error_during_fitting_` is set.

**Critical: jthread destructor semantics.** `std::jthread::~jthread()` calls `request_stop()` then `join()`. If we store jthreads in a vector and let the vector destructor run (or call `.clear()`), each jthread will be asked to stop before join. This is semantically different from the current code where threads run to natural completion and are then joined. We must NOT use destructor-driven stop — we must explicitly join after all work is done.

- [ ] **Step 1: Replace VLA thread array with vector<jthread>, explicit scope for join**

In `reseq/ProbabilityEstimates.cpp` (~lines 1195-1205), replace:
```cpp
if (num_threads > params.size()) {
    num_threads = params.size();
}
thread t[num_threads];
for (auto i = num_threads; i--;) {
    t[i] = std::thread(IPFThread, std::ref(*this), std::cref(stats), std::ref(params), max_iterations,
                       precision_aim / 100);
}
for (auto i = num_threads; i--;) {
    t[i].join();
}
```
with:
```cpp
if (num_threads > params.size()) {
    num_threads = params.size();
}
{
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (auto i = num_threads; i--;) {
        threads.emplace_back([this, &stats, &params, max_iterations, precision_aim](std::stop_token) {
            // stop_token accepted but NOT checked — threads run to natural completion
            decltype(params.size()) cur_par(current_param_++);
            for (; cur_par < params.size() && !error_during_fitting_;
                 cur_par = current_param_++) {
                IterativeProportionalFitting(stats, params.at(cur_par).selected_data, /* ... remaining args ... */);
            }
        });
    }
    // Block scope exit: jthread destructors call request_stop() then join().
    // Workers ignore stop_token — they complete all remaining work before joining.
    // This preserves identical completion semantics to the original std::thread code.
}
```

**Why stop_token is accepted but ignored:** The jthread destructor calls `request_stop()` before `join()`. If workers checked `stop_requested()`, they would terminate early when the vector goes out of scope. By ignoring the token, workers process ALL remaining parameters exactly as the original code does. The `request_stop()` call is harmless — it sets a flag nobody reads.

- [ ] **Step 2: Remove the static IPFThread function**

The lambda in Step 1 replaces the static `IPFThread` function. Remove the static function declaration from `ProbabilityEstimates.h` and the definition from `ProbabilityEstimates.cpp` (~lines 981-989).

Verify the lambda captures all necessary arguments by checking what `IPFThread` passed to `IterativeProportionalFitting`. Copy the exact call signature from the existing `IPFThread` body.

- [ ] **Step 3: Build and run full test suite**

```bash
make build && make test
```
Expected: all 65 tests pass.

- [ ] **Step 4: Run under ASan+UBSan**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 5: Commit**

```bash
git add reseq/ProbabilityEstimates.h reseq/ProbabilityEstimates.cpp
git commit -m "refactor(4b): convert ProbabilityEstimates to std::jthread

Replace VLA std::thread array + manual join with std::vector<std::jthread>.
Workers accept stop_token but do NOT check it — they run to natural
completion, preserving identical semantics to the original code.
jthread destructors handle join automatically when the vector goes
out of scope."
```

---

### Task 3: DataStats — jthread Adoption (Phase 4b)

**Files:**
- Modify: `reseq/DataStats.h` (thread-related members if needed)
- Modify: `reseq/DataStats.cpp:1304-1313` (thread creation/join)

DataStats has a similar pattern to ProbabilityEstimates: VLA thread array, manual join.

- [ ] **Step 1: Examine the ReadThread pattern**

Read `reseq/DataStats.cpp` around line 1304 to understand:
- What `ReadThread` does
- How `running_threads_` and `finish_threads_` coordinate
- Whether there's a stop condition that maps to `stop_token`

- [ ] **Step 2: Replace VLA with vector<jthread>, preserve completion semantics**

Same pattern as Task 2 — accept `stop_token` but don't check it. Workers run to natural completion:
```cpp
{
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    running_threads_ = num_threads;
    finish_threads_ = false;
    for (auto i = num_threads; i--;) {
        threads.emplace_back([this, &bam](std::stop_token) {
            ReadThread(*this, bam);
        });
    }
    // jthread destructors join on scope exit
}
```

- [ ] **Step 3: Build and run full test suite**

```bash
make build && make test
```

- [ ] **Step 4: Commit**

```bash
git add reseq/DataStats.h reseq/DataStats.cpp
git commit -m "refactor(4b): convert DataStats to std::jthread

Replace VLA std::thread array + manual join with std::vector<std::jthread>.
Workers run to natural completion (stop_token accepted but not checked)."
```

---

### Task 4: Add TSan CI Job

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Add ThreadSanitizer CI job**

Add after the `sanitizers` job in `.github/workflows/ci.yml`:
```yaml
  thread-sanitizer:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y gcc-13 g++-13 \
            libboost-all-dev libbz2-dev zlib1g-dev

      - name: Download test data
        run: ./test/download_test_data.sh

      - name: Configure with TSan
        env:
          CC: gcc-13
          CXX: g++-13
        run: |
          cmake -S . -B build \
            -DCMAKE_BUILD_TYPE=Debug \
            -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
            -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"

      - name: Build
        run: cmake --build build -j$(nproc)

      - name: Test with TSan
        run: cd build && ctest --output-on-failure -j1
```

**Note:** TSan and ASan cannot be combined — they are separate CI jobs. TSan may report pre-existing data races. If so, set `continue-on-error: true` initially and document the races for later fixing.

- [ ] **Step 2: Test locally with TSan**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

Document any TSan warnings — they indicate pre-existing races that Tranche 2 should fix.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add ThreadSanitizer (TSan) CI job

TSan runs separately from ASan as they cannot be combined. Verifies
no data races remain after concurrency modernization."
```

---

### Task 5: Fix Tooling Gaps

**Files:**
- Modify: `.clang-format:4` (C++14 → C++20)
- Modify: `.github/workflows/ci.yml` (add clang-tidy job)

**Out of scope for this task:** The GTest PUBLIC link in `reseq/CMakeLists.txt:31` is NOT a cosmetic leak — production headers (`reseq/Vect.hpp:10`, `reseq/DataStats.h:246`, `reseq/FragmentDistributionStats.h:524`, `reseq/FragmentDuplicationStats.h:41`) use `FRIEND_TEST` which requires `gtest/gtest.h`. Removing GTest from the PUBLIC interface requires first removing `FRIEND_TEST` from all production headers and replacing with an alternative test access pattern (e.g., dedicated test accessors, or `#ifdef TESTING`). This is a separate task with its own design considerations. Similarly, GoogleTest is unconditionally fetched at `CMakeLists.txt:27-36` before the `RESEQ_BUILD_TESTS` conditional — gating the fetch also requires removing the public header dependency first.

- [ ] **Step 1: Fix .clang-format C++ standard**

In `.clang-format:4`:
```yaml
# Old:
Standard: c++14
# New:
Standard: c++20
```

- [ ] **Step 2: Check if standard change causes reformatting**

```bash
make format-check
```

If formatting drift is found, run `make format` and include the formatting changes in the commit.

- [ ] **Step 3: Add clang-tidy CI job**

clang-tidy works best with a compile command database for actual translation units. Run on `.cpp` files only — header diagnostics flow through the compile database via `HeaderFilterRegex` in `.clang-tidy:23`.

Add to `.github/workflows/ci.yml`:
```yaml
  lint:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y clang-tidy-17 clang-17 clang++-17 \
            libboost-all-dev libbz2-dev zlib1g-dev

      - name: Generate compile_commands.json
        env:
          CC: clang-17
          CXX: clang++-17
        run: cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

      - name: Run clang-tidy on production sources
        run: |
          git ls-files 'reseq/*.cpp' \
            | grep -v Test | grep -v test_main \
            | xargs clang-tidy-17 -p build
```

**Note:** Start without `--warnings-as-errors` (informational). Header diagnostics are controlled by `HeaderFilterRegex: '.*/reseq/[^/]*\.(h|hpp)$'` in `.clang-tidy:23` and will appear in output when a `.cpp` file includes them.

- [ ] **Step 4: Run clang-tidy locally to assess current state**

```bash
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
git ls-files 'reseq/*.cpp' | grep -v Test | grep -v test_main | head -3 | \
  xargs clang-tidy -p build 2>&1 | tail -20
```

Review output. Do NOT fix warnings in this task — just establish the CI job.

- [ ] **Step 5: Build and test**

```bash
make build && make test
```

- [ ] **Step 6: Commit**

```bash
git add .clang-format .github/workflows/ci.yml
git commit -m "build: update clang-format to C++20 and add clang-tidy CI

Update .clang-format Standard from c++14 to c++20 to match the build.
Add informational clang-tidy CI job running on production .cpp files
with compile_commands.json. Header diagnostics flow via HeaderFilterRegex."
```

---

## Tranche 2: Complex Subsystem Design Documents

These tasks produce **design documents, not code**. Each document must capture the invariants the review identified as missing.

### Task 6: CoverageStats — Design Document

**Files:**
- Create: `docs/superpowers/designs/phase4-coveragestats-concurrency.md`

The design document must address these specific invariants:

- [ ] **Step 1: Document linked-list neighbor semantics**

Read and document every site that uses `previous_block_` and `next_block_`. Key invariants to capture:
- **Gap handling** (`CoverageStats.cpp:558-583`): When `next_block_` starts at a non-contiguous position (gap between blocks), dummy entries are inserted into `previous_coverage_`. The code checks `start_pos_ + kBlockSize == next_block_->start_pos_` to distinguish contiguous from gapped blocks.
- **Reverse traversal** (`CoverageStats.cpp:827-837`): `FindBlock` walks backward via `previous_block_` matching `sequence_id_` and `start_pos_`. This requires `previous_block_` to point to the semantically previous live block, not just any block with a lower index.
- **Inline traversal** (`CoverageStats.h:498+`): Helper functions like `NextBlockWithInSysErrorResetDistance` traverse via `next_block_`.

- [ ] **Step 2: Document reuse and reclamation semantics**

Key sites:
- **Pool acquisition** (`CoverageStats.cpp:447-476`): `CreateBlock` uses `try_lock` on `reuse_mutex_` — if lock fails, allocates new. If lock succeeds and pool non-empty, reuses block and resets state.
- **Detach and reclaim** (`CoverageStats.cpp:965-1025`): After processing, blocks are unlinked from the live list (`first_block_` advances) under `clean_up_mutex_`, then added to reuse pool under `reuse_mutex_`. The two-phase unlock ensures no thread traverses a block being reclaimed.
- **Finalization** (`CoverageStats.cpp:1051-1056`): All reusable blocks deleted at end.

- [ ] **Step 3: Evaluate container options preserving these invariants**

Options:
1. **Keep linked list, add jthread to callers**: Minimal change. CoverageStats doesn't create threads itself — DataStats does. After Task 3, the thread creation is already modernized. Document why container redesign is deferred.
2. **Intrusive list with shared ownership**: Replace raw pointers with `shared_ptr`/`weak_ptr`. Pool stores `shared_ptr`. Neighbor links use `shared_ptr`. Prevents dangling.
3. **Keep raw linked list, guard with reader-writer lock**: Add `std::shared_mutex` around traversal vs modification.

- [ ] **Step 4: Write and commit design document**

```bash
git add docs/superpowers/designs/phase4-coveragestats-concurrency.md
git commit -m "docs(4a): CoverageStats concurrency design — neighbor and reclamation invariants"
```

---

### Task 7: Simulator — Design Document

**Files:**
- Create: `docs/superpowers/designs/phase4-simulator-concurrency.md`

- [ ] **Step 1: Document ownership and lifecycle invariants**

Key invariants the previous plan missed:
- **Reverse blocks built first** (`Simulator.cpp:932-946`): `CreateUnit` builds ALL reverse blocks first, linking them via `next_block_` (in reverse direction) and `partner_block_` (points to previous reverse block). Then forward blocks are created and linked to reverse blocks.
- **partner_block_ carries deletion semantics** (`Simulator.cpp:1327-1328`): When cleaning up, `delete del_block->partner_block_` is called before `delete del_block`. The forward block owns its reverse partner.
- **Hot-path partner traversal** (`Simulator.cpp:645+`): During simulation, `partner_block_` is accessed to get the paired block for the other strand. This is read-only during simulation but the block could be deleted by cleanup running concurrently.

- [ ] **Step 2: Document the Flush multi-mutex cascade**

The Flush pattern (`Simulator.cpp:152-187`):
1. `output_mutex_` locked by caller (Output)
2. Swap buffers under `output_mutex_`, then unlock
3. Lock `flush_mutex_[0]`, write segment 0
4. Lock `flush_mutex_[1]`, unlock `flush_mutex_[0]`
5. Write segment 1, unlock `flush_mutex_[1]`

This ensures both template segments are written in paired-end order. Any redesign must preserve this ordering guarantee.

- [ ] **Step 3: Document the block creation/cleanup race window**

In `CreateBlock` (`Simulator.cpp:1204-1331`):
- Block creation and cleanup both run under `block_creation_mutex_`
- Cleanup traverses from `first_unit_->first_block_` forward, deleting finished blocks
- `GetNextBlock` also runs under `block_creation_mutex_`
- The `deletion_buffer_` / `req_deletion_buffer_` mechanism ensures blocks aren't deleted too early when long deletions span block boundaries

- [ ] **Step 4: Evaluate container options**

Options:
1. **Keep linked lists, adopt jthread**: Replace `thread t[num_threads]` arrays (`Simulator.cpp:2974`, `Simulator.cpp:3128`) with `vector<jthread>`. Keep block/unit structures. Simplest, lowest risk.
2. **Vector with explicit neighbor indices**: Each block stores `next_index_`, `partner_index_` instead of pointers. Cleanup marks blocks as dead rather than deleting. Deferred reclamation. Requires careful analysis of all pointer dereference sites.
3. **Simplify Flush**: Could the 3-mutex cascade be replaced with a single `flush_mutex_`? Profile to check if the concurrent write optimization matters.

- [ ] **Step 5: Write and commit design document**

```bash
git add docs/superpowers/designs/phase4-simulator-concurrency.md
git commit -m "docs(4a): Simulator concurrency design — ownership, Flush, and cleanup invariants"
```

---

### Task 8: FragmentDistributionStats — Design Document

**Files:**
- Create: `docs/superpowers/designs/phase4-fragmentdist-concurrency.md`

- [ ] **Step 1: Document the non-blocking progress guarantee**

The current design ensures forward progress — no worker can be stuck waiting for a resource only it can release:

- **AddNewBiasCalculations** (`FragmentDistributionStats.cpp:2957-2984`): Workers `test_and_set()` on queue slots. If ALL slots claimed, worker **breaks** (not blocks) and continues draining.
- **ExecuteBiasCalculations** (`FragmentDistributionStats.cpp:2987-3035`): Workers process items from claimed bins and release them when done.

A naive `condvar::wait()` on a full queue deadlocks: if all workers block, no one completes work to free slots.

- [ ] **Step 2: Document publication and reuse rules**

Key ordering constraints:
- **Publication** (`FragmentDistributionStats.cpp:2238-2242`): `current_bias_param_` is set to 0 LAST — this allows other threads to start using the completed queue bin. Setting it earlier would let threads read incomplete data.
- **Slot release** (`FragmentDistributionStats.cpp:3024-3032`): `claimed_bias_bins_.at(queue_bin).clear()` only after ALL calculations in that bin are finished (`finished_bias_calcs_ == queue_bin_size`). Premature release would let `AddNewBiasCalculations` overwrite in-use data.
- **Resource acquisition** (`FragmentDistributionStats.cpp:3002-3004`): `bias_calc_vects_` spin-wait to acquire a temporary vector for bias calculation. These are separate from queue bins.

- [ ] **Step 3: Evaluate alternatives preserving progress**

Options:
1. **Keep atomic_flag spin-waits, adopt jthread**: Minimal change. Replace `std::thread` arrays (`FragmentDistributionStats.cpp:2764`, `FragmentDistributionStats.cpp:3724`) with `vector<jthread>`. Spin-waits are cheap — contention is brief and bounded by `kMaxBinsQueuedForBiasCalc`. Document the deliberate design.
2. **try_lock with fallback**: Replace `atomic_flag::test_and_set` with `std::mutex::try_lock`. If lock fails, skip (same non-blocking behavior). More idiomatic but functionally identical.
3. **Bounded channel with try_send/try_recv**: A bounded MPMC channel where `try_send` returns false when full. Worker breaks on false. Preserves non-blocking semantics.

- [ ] **Step 4: Write and commit design document**

```bash
git add docs/superpowers/designs/phase4-fragmentdist-concurrency.md
git commit -m "docs(4a): FragmentDistributionStats concurrency design — progress and publication invariants"
```

---

## Tranche 2 Implementation (Tasks 9+)

After design documents are reviewed and approved, create follow-up implementation tasks. These will be planned in a separate document after designs are finalized.

---

## Verification Protocol

**Test execution:** Always run the full test suite as a single invocation:
```bash
make test
# or equivalently:
cd build && ctest --output-on-failure -j1
```
Do NOT use `--gtest_filter` as a reliable checkpoint — tests have inter-test dependencies (shared `Register()` state) that require running together (`reseq/CMakeLists.txt:60`).

**Regression tests:** The 10 golden-file regression tests verify byte-identical output. They run as part of the full suite.

**Sanitizer verification:** After each task:
```bash
# ASan+UBSan
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1

# TSan (after Task 4 is merged)
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```
