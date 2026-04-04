# Phase 5b: DataStats Decomposition — Design Specification

**Date:** 2026-04-03
**Status:** Draft (rev 2 — incorporates code review findings)
**Scope:** Decompose the `DataStats` god object (1,844 lines) into focused, single-responsibility components
**Prerequisite:** Phase 5a (main.cpp decomposition) — merged
**Input:** Re-evaluation of original Phase 5 design against current codebase (post-Phases 3, 4, 6)

---

## Goal

Break `DataStats` into three components:
1. **ReadSequenceStats** — read-level histogram statistics produced solely by `EvalBaseLevelStats`
2. **BamIngestionEngine** — multi-threaded two-pass BAM reading pipeline with pair-matching, record evaluation, and validation
3. **DataStats** (thinned) — aggregator that owns sub-stats, exposes getters, handles serialization and lifecycle

The public API of `DataStats` does not change. External callers (CLI commands, Simulator, ProbabilityEstimates) continue to use it identically. Internal test code that accesses private members via friend will need minor path changes (see Testing section).

---

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| RecordClassifier extraction | No — keep predicates in BamIngestionEngine | Small inline predicates only called from within the BAM loop; separate class adds indirection without value |
| ReadSequenceStats scope | Only histograms where `EvalBaseLevelStats` is the sole writer | Reference-derived histograms (`tmp_sequence_content_reference_`, `tmp_gc_read_content_reference_`, `tmp_gc_read_content_mapped_`) are written by `EvalReferenceStatistics` and stay on DataStats |
| QualityStats coupling | `EvalBaseLevelStats` continues to call `qualities_.AddRawBase()` and `qualities_.AddRawHomoqualimer()` by receiving `QualityStats&` as a parameter | ReadSequenceStats does not own QualityStats; the BAM engine passes it in |
| BamIngestionEngine ownership model | Receives sub-stats by reference, does not own them | Avoids moving ownership; keeps sub-stat lifecycle under DataStats |
| BamFileMetadata struct | No — keep scalars in DataStats | Only 6 scalar fields; a struct adds a layer without meaningful encapsulation benefit |
| GetReadPosOnReference | Stays on DataStats (non-static, depends on `reference_`) | Cannot move to BamIngestionEngine's public interface since it requires instance state |
| Public API stability | No changes to DataStats public interface | All external callers work unchanged |
| Friend access | BamIngestionEngine and ReadSequenceStats are friends of DataStats; DataStatsTest gets friend access to ReadSequenceStats | Necessary during transition; can be narrowed in Phase 6 |
| Helper types | `RecordHasher`, `RecordEqual` move to BamIngestionEngine | They are only used by `first_read_records_` which moves to the engine |

---

## Component 1: ReadSequenceStats

### Responsibility

Accumulates read-level sequence statistics during the BAM scan via `EvalBaseLevelStats` and exposes them as histograms afterward. This class is a data sink — it does not orchestrate or coordinate other sub-stats.

### Write Ownership Analysis

`EvalBaseLevelStats` (DataStats.cpp:154-223) writes to the following `tmp_*` members. Only members where it is the **sole writer** move to ReadSequenceStats:

| Member | Writer | Moves? |
|--------|--------|--------|
| `tmp_sequence_content_` | `EvalBaseLevelStats` only | Yes |
| `tmp_homopolymer_distribution_` | `EvalBaseLevelStats` only | Yes |
| `tmp_gc_read_content_` | `EvalBaseLevelStats` only | Yes |
| `tmp_n_content_` | `EvalBaseLevelStats` only | Yes |
| `tmp_proper_pair_mapping_quality_` | `EvalRecord` (not EvalBaseLevelStats) | Yes (sole writer is EvalRecord, but it's a simple `++` in a single location) |
| `tmp_improper_pair_mapping_quality_` | `EvalRecord` | Yes (same pattern) |
| `tmp_single_read_mapping_quality_` | `EvalRecord` | Yes (same pattern) |
| `tmp_sequence_content_reference_` | `EvalReferenceStatistics` | **No — stays on DataStats** |
| `tmp_gc_read_content_reference_` | `EvalReferenceStatistics` | **No — stays on DataStats** |
| `tmp_gc_read_content_mapped_` | `EvalReferenceStatistics` | **No — stays on DataStats** |

`EvalBaseLevelStats` also mutates external state:
- Calls `qualities_.AddRawBase()`, `qualities_.AddRawHomoqualimer()`, `qualities_.AddRawRead()` — these are QualityStats mutations
- Sets `full_record->sequence_quality_` from computed `seq_qual_stats.mean_`

These external mutations are handled by passing `QualityStats&` and `FullRecord*` as parameters to `EvalBaseLevelStats`, not by owning them.

### Data Members (moved from DataStats)

**Atomic accumulators (thread-safe, filled during multi-threaded BAM pass):**
- `tmp_proper_pair_mapping_quality_` — `std::vector<utilities::VectorAtomic<uintFragCount>>`
- `tmp_improper_pair_mapping_quality_` — same type
- `tmp_single_read_mapping_quality_` — same type
- `tmp_gc_read_content_` — `std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2>`
- `tmp_n_content_` — same type
- `tmp_sequence_content_` — `std::array<std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5>, 2>`
- `tmp_homopolymer_distribution_` — `std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5>`

**Final histograms (serialized):**
- `proper_pair_mapping_quality_` — `Vect<uintFragCount>`
- `improper_pair_mapping_quality_` — same type
- `single_read_mapping_quality_` — same type
- `gc_read_content_` — `std::array<Vect<uintFragCount>, 2>`
- `n_content_` — same type
- `sequence_content_` — `std::array<std::array<Vect<uintNucCount>, 5>, 2>`
- `homopolymer_distribution_` — `std::array<Vect<uintNucCount>, 5>`

**Members that stay on DataStats** (written by `EvalReferenceStatistics`):
- `tmp_sequence_content_reference_`, `sequence_content_reference_`
- `tmp_gc_read_content_reference_`, `gc_read_content_reference_`
- `tmp_gc_read_content_mapped_`, `gc_read_content_mapped_`

### Methods

- `PrepareAccumulators(max_read_len, max_qual)` — allocates all `tmp_*` vectors to correct sizes
- `EvalBaseLevelStats(full_record, template_segment, strand, tile_id, paired_seq_qual, qualities, phred_offset, max_quality)` — moved from DataStats. Receives `QualityStats&` as parameter for the `AddRawBase`/`AddRawHomoqualimer`/`AddRawRead` calls. Writes to own `tmp_*` members and sets `full_record->sequence_quality_`.
- `IncrementMappingQuality(quality_type, value)` — thin wrapper for the mapping quality `++` operations currently inline in `EvalRecord`. Called by the engine, not by `EvalBaseLevelStats`.
- `Finalize()` — copies atomic accumulators into final `Vect` structures
- `Shrink()` — trims histogram ranges

### Const Getters

Same signatures as current DataStats getters:
- `ProperPairMappingQuality()`, `ImproperPairMappingQuality()`, `SingleReadMappingQuality()`
- `GCReadContent(template_segment)`, `NContent(template_segment)`
- `SequenceContent(template_segment, nucleotide)`
- `HomopolymerDistribution(nucleotide)`

### Serialization

Boost `serialize()` template covering all final histogram members. Called as part of DataStats's own serialization (DataStats serializes `read_sequence_stats_` as a member).

### Dependencies

- `Vect<T>` (histogram container)
- `utilities::VectorAtomic<T>` (thread-safe accumulation)
- `CoverageStats::FullRecord` (writes `sequence_quality_` field)
- `QualityStats&` (passed as parameter, called for base/read-level quality accumulation)
- `seqan::BamAlignmentRecord` (read from)

### Friend Declarations

```cpp
friend class DataStatsTest;  // tests access histograms directly
friend class boost::serialization::access;
```

---

## Component 2: BamIngestionEngine

### Responsibility

Executes the two-pass BAM reading pipeline: pre-run (pass 1 for sizing), multi-threaded evaluation (pass 2), and finalization. Populates a DataStats object and its sub-stats by reference.

### Data Members (moved from DataStats)

**Threading primitives:**
- `read_mutex_`, `print_mutex_` — `std::mutex`
- `running_threads_` — `std::atomic<uintNumThreads>`
- `finish_threads_mutex_` — `std::mutex`
- `finish_threads_cv_` — `std::condition_variable`
- `finish_threads_` — `bool`

**Pair-matching state and helper types:**
- `RecordHasher` struct (moved from DataStats) — hash function for FullRecord pointers
- `RecordEqual` struct (moved from DataStats) — equality comparison for FullRecord pointers
- `first_read_records_` — `std::unordered_set<CoverageStats::FullRecord*, RecordHasher, RecordEqual>`
- `reading_success_` — `std::atomic<bool>`
- `read_records_` — `uintFragCount`
- `kBatchSize` — `const uintFragCount` (10000)

**Pre-run temporaries:**
- `reads_per_frag_len_bin_` — `std::vector<uintFragCount>`
- `lowq_reads_per_frag_len_bin_` — `std::vector<uintFragCount>`

**Read-disposition counters:**
- `reads_in_unmapped_pairs_without_adapters_` — `std::atomic<uintFragCount>`
- `reads_in_unmapped_pairs_with_adapters_` — same type
- `reads_with_low_quality_with_adapters_` — same type
- `reads_with_low_quality_without_adapters_` — same type
- `reads_on_too_short_fragments_` — same type
- `reads_in_excluded_regions_` — same type
- `reads_used_` — same type

### Methods (moved from DataStats)

**Pipeline orchestration:**
- `Run(bam_file, adapter_file, adapter_matrix, variant_file, max_ref_seq_bin_size, num_threads, calculate_bias, DataStats& target)` — replaces `DataStats::ReadBam` as the implementation
- `PreRun(bam, bam_file, header, ...)` — first pass: read lengths, quality range, adapter evidence, tile enumeration
- `ReadRecords(bam, target, ...)` — locked batch reader; pairs reads, registers with CoverageStats blocks
- `ReadThread(engine, bam, thread_id)` — static worker entry point
- `PrepareReadIn(target, ...)` — allocates all temporary arrays on DataStats and sub-stats
- `FinishReadIn(target)` — copies atomics to finals, finalizes sub-stats
- `Shrink(target)` — trims all histogram arrays
- `Calculate(target, num_threads)` — triggers bias calculation and computes corrected coverage

**Per-record evaluation:**
- `EvalRecord(...)` — dispatches one read pair through filters and statistics. Calls `read_sequence_stats.EvalBaseLevelStats(...)` for read-level stats, `EvalReferenceStatistics(...)` for reference-level stats, and directly increments mapping quality via `read_sequence_stats.IncrementMappingQuality(...)`.
- `EvalReferenceStatistics(...)` — CIGAR walk. Writes to `tmp_sequence_content_reference_`, `tmp_gc_read_content_reference_`, `tmp_gc_read_content_mapped_` (which remain on DataStats, accessed via friend). Also calls `errors_.AddBasePlotting()`, `errors_.AddInDel()`, and coverage block methods.
- `CheckForAdapters(...)` — detects adapter contamination in unmapped/low-quality pairs
- `IsSecondRead(...)` — uses hash set to match mates

**Validation predicates (all inline):**
- `IsValidRecord(record)` — paired, consistent lengths, has sequence
- `QualitySufficient(record, min_quality)` — mapQ threshold (takes threshold as parameter instead of reading from member)
- `PotentiallyValidGeneral(...)`, `PotentiallyValidFirst(...)`, `PotentiallyValidSecond(...)`, `PotentiallyValid(...)` — multi-condition filters
- `InProperDirection(...)` — FR orientation and insert-length check
- `GetReadLengthOnReference(record)` (x2 static) — CIGAR-corrected length

**Cleanup:**
- `OrderOfBamFileCorrect(...)` — validates position-sort order
- `SignsOfPairsWithNamesNotIdentical(...)` — detects missing mates

### Interface

```cpp
class BamIngestionEngine {
public:
    // Hash/equality functors for FullRecord pair-matching
    struct RecordHasher {
        size_t operator()(CoverageStats::FullRecord* const& record) const;
    };
    struct RecordEqual {
        bool operator()(CoverageStats::FullRecord* const& lhs,
                       CoverageStats::FullRecord* const& rhs) const;
    };

    bool Run(const char* bam_file, const char* adapter_file,
             const char* adapter_matrix, const std::string& variant_file,
             uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
             bool calculate_bias, DataStats& target);

    // Static geometry utilities
    static uintReadLen GetReadLengthOnReference(
        const seqan::BamAlignmentRecord& record);
    static uintReadLen GetReadLengthOnReference(
        const seqan::BamAlignmentRecord& record, uintReadLen& max_indel);
};
```

All other methods are private. The engine is constructed on the stack inside `DataStats::ReadBam`, runs, and is destroyed.

**Note on GetReadPosOnReference:** This method is non-static — it depends on `reference_` for the soft-clip correction. It stays on DataStats. No external callers use it (confirmed by grep). BamIngestionEngine accesses it via the DataStats target reference during `Run`.

### Dependencies

- `DataStats` (friend, writes to target's members and sub-stats)
- `ReadSequenceStats` (calls `EvalBaseLevelStats`, `PrepareAccumulators`, `Finalize`, `Shrink`, `IncrementMappingQuality`)
- All sub-stats classes via DataStats: `AdapterStats`, `CoverageStats`, `ErrorStats`, `FragmentDistributionStats`, `FragmentDuplicationStats`, `QualityStats`, `TileStats`
- `Reference` (via DataStats's reference pointer)
- SeqAn BAM I/O

---

## Component 3: DataStats (thinned)

### What Remains

**Config:**
- `reference_` — `Reference*`
- `maximum_insert_length_` — `const uintSeqLen`
- `minimum_mapping_quality_` — `const uintQual`

**Owned sub-stats (unchanged):**
- `adapters_`, `coverage_`, `duplicates_`, `errors_`, `fragment_distribution_`, `qualities_`, `tiles_`

**New owned component:**
- `read_sequence_stats_` — `ReadSequenceStats`

**Metadata scalars (unchanged):**
- `creation_time_`, `phred_quality_offset_`, `minimum_quality_`, `maximum_quality_`
- `minimum_read_length_on_reference_`, `maximum_read_length_on_reference_`
- `corrected_coverage_`, `total_number_reads_`

**Read-length histograms (unchanged):**
- `read_lengths_`, `read_lengths_by_fragment_length_`, `non_mapped_read_lengths_by_fragment_length_`

**Reference-derived histograms (stay here, written by EvalReferenceStatistics via friend):**
- `tmp_sequence_content_reference_`, `sequence_content_reference_`
- `tmp_gc_read_content_reference_`, `gc_read_content_reference_`
- `tmp_gc_read_content_mapped_`, `gc_read_content_mapped_`

**Non-static method that stays:**
- `GetReadPosOnReference(start, end, record)` — depends on `reference_`

### Public API (unchanged for external callers)

All existing getters remain with identical signatures. Getters for ReadSequenceStats data forward:

```cpp
const Vect<uintFragCount>& ProperPairMappingQuality() const {
    return read_sequence_stats_.ProperPairMappingQuality();
}
// ... same pattern for all moved histogram getters
```

Getters for reference-derived histograms remain direct (data stays on DataStats):
```cpp
const Vect<uintNucCount>& SequenceContentReference(...) const {
    return sequence_content_reference_...;  // unchanged
}
```

`ReadBam(...)` becomes:
```cpp
bool DataStats::ReadBam(...) {
    BamIngestionEngine engine;
    return engine.Run(bam_file, adapter_file, adapter_matrix,
                      variant_file, max_ref_seq_bin_size, num_threads,
                      calculate_bias, *this);
}
```

`Load`/`Save` unchanged — `serialize()` includes `read_sequence_stats_` as a member.

### Friend Declarations

```cpp
friend class BamIngestionEngine;      // writes to sub-stats and ref-derived histograms during Run()
friend class ReadSequenceStats;       // not needed (ReadSequenceStats doesn't access DataStats privates)
// Existing friends remain: DataStatsTest, ProbabilityEstimatesTest, SimulatorTest
```

---

## Testing Impact

The spec previously claimed "all tests pass unchanged." This is true for **external callers** but not for **DataStatsTest**, which accesses private members directly via friend.

### Members accessed by DataStatsTest that move to ReadSequenceStats

Tests in `DataStatsTest.cpp` directly access (via `test_->` with friend access):
- `sequence_content_` (lines 60-73)
- `proper_pair_mapping_quality_` (line 120)
- `improper_pair_mapping_quality_`, `single_read_mapping_quality_`
- `gc_read_content_`, `n_content_`, `homopolymer_distribution_`

### Migration path

These test accesses change from `test_->sequence_content_` to `test_->read_sequence_stats_.sequence_content_`. This is mechanical — DataStatsTest already has friend access to DataStats, and `read_sequence_stats_` is a private member of DataStats, so the test can reach through. ReadSequenceStats also declares `friend class DataStatsTest` so the test can access its private histogram members.

No test logic changes. No new assertions. No new test files. Just member access path updates.

### Members accessed by DataStatsTest that stay on DataStats

- `sequence_content_reference_`, `gc_read_content_reference_`, `gc_read_content_mapped_` — these stay, so test code for them is unchanged.
- All metadata scalars (`minimum_quality_`, etc.) — unchanged.

---

## File Layout

| File | Status | Content |
|------|--------|---------|
| `reseq/ReadSequenceStats.h` | New | Class declaration, inline getters, friend declarations |
| `reseq/ReadSequenceStats.cpp` | New | PrepareAccumulators, EvalBaseLevelStats, Finalize, Shrink, IncrementMappingQuality |
| `reseq/BamIngestionEngine.h` | New | Class declaration with RecordHasher/RecordEqual, public Run + static utilities |
| `reseq/BamIngestionEngine.cpp` | New | All pipeline methods, threading, evaluation, validation predicates |
| `reseq/DataStats.h` | Modified | Remove extracted members/methods, add read_sequence_stats_ member, add friend BamIngestionEngine, keep ref-derived histograms |
| `reseq/DataStats.cpp` | Modified | Remove extracted method implementations, ReadBam becomes thin wrapper, keep EvalReferenceStatistics-related finalization |
| `reseq/DataStatsTest.cpp` | Modified | Update member access paths from `test_->X` to `test_->read_sequence_stats_.X` for moved histograms |
| `reseq/CMakeLists.txt` | Modified | Add ReadSequenceStats.cpp, BamIngestionEngine.cpp to reseq_lib |

---

## Extraction Order

Each step must build and pass all tests before proceeding:

1. **Extract ReadSequenceStats** — move histogram data + EvalBaseLevelStats + getters + serialize. DataStats forwards getters. Update DataStatsTest member access paths. Reference-derived histograms stay on DataStats.
2. **Extract BamIngestionEngine** — move pipeline methods + threading + predicates + counters + RecordHasher/RecordEqual. DataStats::ReadBam becomes thin wrapper.
3. **Clean up DataStats** — remove dead includes, verify all forwarding getters, tidy friend declarations.
4. **Format + full test suite** — `make format`, `make build`, `make test`, `pre-commit run --all-files`.

---

## Verification Gate

- All existing unit tests pass (DataStatsTest with updated access paths, ProbabilityEstimatesTest, SimulatorTest)
- All regression tests pass unchanged (golden file comparisons)
- `make format-check` clean
- `pre-commit run --all-files` clean
- No public API changes — grep for `DataStats::` usage across codebase confirms no external caller changes needed

---

## What This Does NOT Do

- Does not add new unit test files (deferred to Phase 7)
- Does not remove FRIEND_TEST or narrow friend access (deferred to Phase 6)
- Does not change any public API signatures
- Does not touch other god objects (Phases 5c-f)
- Does not refactor sub-stats classes (they are already single-responsibility)
- Does not move reference-derived histograms — their writer (`EvalReferenceStatistics`) is deeply coupled to coverage blocks and error stats
