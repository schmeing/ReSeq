# Phase 5b: DataStats Decomposition — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decompose the DataStats god object into ReadSequenceStats, BamIngestionEngine, and a thinned DataStats aggregator while preserving the public API and all test behavior.

**Architecture:** ReadSequenceStats owns read-level histograms and EvalBaseLevelStats. BamIngestionEngine owns the BAM reading pipeline, threading, and record evaluation. DataStats remains the public facade, forwarding getters and delegating ReadBam to the engine.

**Tech Stack:** C++20, Boost.Serialization, SeqAn BAM I/O, GoogleTest

**Design spec:** `docs/superpowers/specs/2026-04-03-phase5b-datastats-decomposition-design.md`

**Prerequisite:** Phase 5a (main.cpp decomposition) merged. `make build && make test` passes.

---

## Critical Invariants

1. **Public API of DataStats does NOT change.** All external callers (Simulator, ProbabilityEstimates, CLI commands, DataStatsInterface) continue using DataStats identically.
2. **All existing tests must pass.** DataStatsTest accesses private members via friend; access paths change mechanically.
3. **Build and test after every task.** `make build && make test` must pass before committing.
4. **Conventional commits** with `(5b)` scope tag.
5. **Boost serialization backward compatibility.** The new `read_sequence_stats_` member is serialized inline within DataStats's existing `serialize()` template so that the on-disk format of the histogram fields remains identical. The serialization order of all fields MUST be preserved exactly.

---

## Task 1: Extract ReadSequenceStats

**Goal:** Create `ReadSequenceStats` class owning read-level histogram data and `EvalBaseLevelStats`. DataStats gets `read_sequence_stats_` member, forwards getters, includes it in serialization.

**Files created:**
- `reseq/ReadSequenceStats.h`
- `reseq/ReadSequenceStats.cpp`

**Files modified:**
- `reseq/DataStats.h`
- `reseq/DataStats.cpp`
- `reseq/DataStatsTest.cpp`
- `reseq/CMakeLists.txt`

### Step-by-step

- [ ] **Step 1.1: Create `reseq/ReadSequenceStats.h`**

```cpp
#ifndef READSEQUENCESTATS_H
#define READSEQUENCESTATS_H

#include <array>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include "CoverageStats.h"
#include "QualityStats.h"
#include "Vect.hpp"
#include "utilities.hpp"

namespace reseq {

class ReadSequenceStats {
  private:
    // Atomic accumulators (thread-safe, filled during multi-threaded BAM pass)
    std::vector<utilities::VectorAtomic<uintFragCount>> tmp_proper_pair_mapping_quality_;
    std::vector<utilities::VectorAtomic<uintFragCount>> tmp_improper_pair_mapping_quality_;
    std::vector<utilities::VectorAtomic<uintFragCount>> tmp_single_read_mapping_quality_;

    std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_gc_read_content_;
    std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_n_content_;
    std::array<std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5>, 2> tmp_sequence_content_;
    std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5> tmp_homopolymer_distribution_;

    // Final histograms (serialized)
    Vect<uintFragCount> proper_pair_mapping_quality_;
    Vect<uintFragCount> improper_pair_mapping_quality_;
    Vect<uintFragCount> single_read_mapping_quality_;

    std::array<Vect<uintFragCount>, 2> gc_read_content_;
    std::array<Vect<uintFragCount>, 2> n_content_;
    std::array<std::array<Vect<uintNucCount>, 5>, 2> sequence_content_;
    std::array<Vect<uintNucCount>, 5> homopolymer_distribution_;

    // Boost serialization
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned int UNUSED(version)) {
        ar & proper_pair_mapping_quality_;
        ar & improper_pair_mapping_quality_;
        ar & single_read_mapping_quality_;

        ar & gc_read_content_;
        ar & n_content_;
        ar & sequence_content_;

        ar & homopolymer_distribution_;
    }

    // Google test
    friend class DataStatsTest;

  public:
    // Accumulator lifecycle
    void PrepareAccumulators(uintQual size_mapping_quality, uintReadLen size_pos);
    void EvalBaseLevelStats(CoverageStats::FullRecord* full_record, uintTempSeq template_segment,
                            uintTempSeq strand, uintTileId tile_id, uintQual& paired_seq_qual,
                            QualityStats& qualities, uintQual phred_quality_offset, uintQual maximum_quality);
    void IncrementMappingQuality(uintTempSeq quality_type, uintQual value);
    void Finalize();
    void Shrink();

    // Const getters (same signatures as original DataStats getters)
    const Vect<uintFragCount>& ProperPairMappingQuality() const { return proper_pair_mapping_quality_; }
    const Vect<uintFragCount>& ImproperPairMappingQuality() const { return improper_pair_mapping_quality_; }
    const Vect<uintFragCount>& SingleReadMappingQuality() const { return single_read_mapping_quality_; }

    const Vect<uintFragCount>& GCReadContent(uintTempSeq template_segment) const {
        return gc_read_content_.at(template_segment);
    }
    const Vect<uintFragCount>& NContent(uintTempSeq template_segment) const {
        return n_content_.at(template_segment);
    }
    const Vect<uintNucCount>& SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const {
        return sequence_content_.at(template_segment).at(nucleotide);
    }
    const Vect<uintNucCount>& HomopolymerDistribution(uintBaseCall nucleotide) const {
        return homopolymer_distribution_.at(nucleotide);
    }
};

} // namespace reseq

#endif // READSEQUENCESTATS_H
```

Key points:
- The `serialize()` field order MUST match the order these fields appeared in DataStats's `serialize()`. We will handle this in Step 1.5 by serializing `read_sequence_stats_` inline.
- `friend class DataStatsTest` gives tests direct member access.
- `EvalBaseLevelStats` takes `QualityStats&`, `phred_quality_offset`, and `maximum_quality` as parameters since those stay on DataStats.

- [ ] **Step 1.2: Create `reseq/ReadSequenceStats.cpp`**

```cpp
#include "ReadSequenceStats.h"
using reseq::ReadSequenceStats;

#include <array>
using std::array;

#include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::Dna5;

#include "SeqQualityStats.hpp"
// include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::Complement;
using reseq::utilities::Percent;

void ReadSequenceStats::PrepareAccumulators(uintQual size_mapping_quality, uintReadLen size_pos) {
    tmp_proper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_improper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_single_read_mapping_quality_.resize(size_mapping_quality);

    for (auto template_segment = 2; template_segment--;) {
        tmp_gc_read_content_.at(template_segment).resize(101);
        tmp_n_content_.at(template_segment).resize(101);
        for (auto base = 5; base--;) {
            tmp_sequence_content_.at(template_segment).at(base).resize(size_pos);
        }
    }

    for (auto base = 5; base--;) {
        tmp_homopolymer_distribution_.at(base).resize(size_pos);
    }
}

void ReadSequenceStats::EvalBaseLevelStats(CoverageStats::FullRecord* full_record, uintTempSeq template_segment,
                                           uintTempSeq strand, uintTileId tile_id, uintQual& paired_seq_qual,
                                           QualityStats& qualities, uintQual phred_quality_offset,
                                           uintQual maximum_quality) {
    const BamAlignmentRecord& record(full_record->record_);

    uintSeqLen pos_reversed(length(record.seq));

    SeqQualityStats<uintNucCount> seq_qual_stats;
    seq_qual_stats[maximum_quality]; // Resize to maximum quality
    for (auto qual : record.qual) {
        ++seq_qual_stats.at(qual - phred_quality_offset);
    }
    seq_qual_stats.Calculate();
    seq_qual_stats.CalculateProbabilityMean();
    full_record->sequence_quality_ = seq_qual_stats.mean_;

    Dna5 base;
    uintQual qual, last_qual(1);
    array<uintSeqLen, 5> read_bases = {0, 0, 0, 0, 0};
    uintReadLen homoquality_length = 0;
    uintReadLen homopolymer_length = 0;
    Dna5 homopolymer_nucleotide = 0;

    for (uintReadLen pos = 0; pos < length(record.seq); ++pos) {
        if (hasFlagRC(record)) {
            base = Complement::Dna5(at(record.seq, --pos_reversed));
            qual = at(record.qual, pos_reversed) - phred_quality_offset;
        } else {
            base = at(record.seq, pos);
            qual = at(record.qual, pos) - phred_quality_offset;
        }

        qualities.AddRawBase(template_segment, base, tile_id, strand, qual, seq_qual_stats.mean_, last_qual, pos);

        ++read_bases.at(base);
        ++tmp_sequence_content_.at(template_segment).at(base).at(pos);

        if (last_qual == qual) {
            ++homoquality_length;
        } else {
            if (pos) {
                qualities.AddRawHomoqualimer(last_qual, homoquality_length);
            }
            homoquality_length = 1;
        }

        if (homopolymer_nucleotide == base) {
            ++homopolymer_length;
        } else {
            ++tmp_homopolymer_distribution_.at(homopolymer_nucleotide).at(homopolymer_length);
            homopolymer_nucleotide = base;
            homopolymer_length = 1;
        }

        last_qual = qual;
    }

    qualities.AddRawHomoqualimer(last_qual, homoquality_length);
    ++tmp_homopolymer_distribution_.at(homopolymer_nucleotide)
          .at(homopolymer_length); // Add the homopolymer at read end

    // Read level summaries of base level stats
    qualities.AddRawRead(paired_seq_qual, seq_qual_stats, template_segment, tile_id, read_bases, length(record.seq));

    ++tmp_gc_read_content_.at(template_segment)
          .at(Percent(read_bases.at(1) + read_bases.at(2), static_cast<uintSeqLen>(length(record.seq))));
    ++tmp_n_content_.at(template_segment).at(Percent(read_bases.at(4), static_cast<uintSeqLen>(length(record.seq))));
}

void ReadSequenceStats::IncrementMappingQuality(uintTempSeq quality_type, uintQual value) {
    // quality_type: 0 = proper pair, 1 = improper pair, 2 = single read
    switch (quality_type) {
    case 0:
        ++tmp_proper_pair_mapping_quality_.at(value);
        break;
    case 1:
        ++tmp_improper_pair_mapping_quality_.at(value);
        break;
    case 2:
        ++tmp_single_read_mapping_quality_.at(value);
        break;
    }
}

void ReadSequenceStats::Finalize() {
    proper_pair_mapping_quality_.Acquire(tmp_proper_pair_mapping_quality_);
    improper_pair_mapping_quality_.Acquire(tmp_improper_pair_mapping_quality_);
    single_read_mapping_quality_.Acquire(tmp_single_read_mapping_quality_);

    for (auto template_segment = 2; template_segment--;) {
        gc_read_content_.at(template_segment).Acquire(tmp_gc_read_content_.at(template_segment));
        n_content_.at(template_segment).Acquire(tmp_n_content_.at(template_segment));
        for (auto base = 5; base--;) {
            sequence_content_.at(template_segment)
                .at(base)
                .Acquire(tmp_sequence_content_.at(template_segment).at(base));
        }
    }

    tmp_homopolymer_distribution_.at(0).at(0) = 0; // Remove homopolymers introduced by initialization values
    for (auto base = 5; base--;) {
        homopolymer_distribution_.at(base).Acquire(tmp_homopolymer_distribution_.at(base));
    }
}

void ReadSequenceStats::Shrink() {
    // No additional shrinking needed for these histograms beyond what Acquire already does.
    // This method exists for consistency with the extraction pattern.
}
```

Key points:
- `EvalBaseLevelStats` is a direct move of `DataStats::EvalBaseLevelStats` (lines 153-222 of DataStats.cpp), with `phred_quality_offset_` and `maximum_quality_` changed from member access to parameters, and `qualities_` changed from member to `qualities` parameter.
- `IncrementMappingQuality` encapsulates the `++tmp_*_mapping_quality_.at(value)` operations currently inline in `EvalRecord`.
- `Finalize` contains the relevant portions from `DataStats::FinishReadIn` (lines 690-733) that handle the moved histograms.
- `PrepareAccumulators` contains the relevant portions from `DataStats::PrepareReadIn` (lines 651-680) that allocate the moved temporaries.

- [ ] **Step 1.3: Add ReadSequenceStats.cpp to CMakeLists.txt**

In `reseq/CMakeLists.txt`, add `ReadSequenceStats.cpp` to the `reseq_lib` sources list (after `QualityStats.cpp`, line 20):

```cmake
add_library(reseq_lib STATIC
  AdapterStats.cpp
  archive_format.cpp
  CoverageStats.cpp
  DataStats.cpp
  DataStatsInterface.cpp
  ErrorStats.cpp
  FragmentDistributionStats.cpp
  FragmentDuplicationStats.cpp
  ProbabilityEstimates.cpp
  QualityStats.cpp
  ReadSequenceStats.cpp
  Reference.cpp
  Simulator.cpp
  Surrounding.cpp
  TileStats.cpp
)
```

- [ ] **Step 1.4: Modify `reseq/DataStats.h`**

**1.4a: Add include** at the top of includes (after the existing includes, before namespace):
```cpp
#include "ReadSequenceStats.h"
```

**1.4b: Add member.** In the private section, after the existing sub-stats block (after `TileStats tiles_;`, approximately line 88), add:
```cpp
    ReadSequenceStats read_sequence_stats_;
```

**1.4c: Remove moved data members.** Delete these private member declarations (they now live on ReadSequenceStats):
- `tmp_proper_pair_mapping_quality_` (line 100)
- `tmp_improper_pair_mapping_quality_` (line 101)
- `tmp_single_read_mapping_quality_` (line 102)
- `tmp_gc_read_content_` (line 104)
- `tmp_n_content_` (line 108)
- `tmp_sequence_content_` (line 109)
- `tmp_homopolymer_distribution_` (line 113)
- `proper_pair_mapping_quality_` (lines 139-140)
- `improper_pair_mapping_quality_` (lines 141-142)
- `single_read_mapping_quality_` (lines 143-144)
- `gc_read_content_` (line 146)
- `n_content_` (line 152)
- `sequence_content_` (lines 153-155)
- `homopolymer_distribution_` (lines 160-161)

Keep these members on DataStats (written by EvalReferenceStatistics, NOT EvalBaseLevelStats):
- `tmp_gc_read_content_reference_` (line 105)
- `tmp_gc_read_content_mapped_` (line 106)
- `tmp_sequence_content_reference_` (lines 110-111)
- `gc_read_content_reference_` (lines 147-148)
- `gc_read_content_mapped_` (lines 149-150)
- `sequence_content_reference_` (lines 156-158)

**1.4d: Remove moved method declaration.** Delete the `EvalBaseLevelStats` declaration (line 183-184). It is now on ReadSequenceStats.

**1.4e: Add friend declaration.** In the friend section (after existing friends, around line 243), add:
```cpp
    friend class BamIngestionEngine;  // Added for Phase 5b Task 2
```

**1.4f: Update serialize().** Replace the serialization of the moved histogram fields with inline serialization of `read_sequence_stats_`. The order of all other fields MUST remain identical. The key change is to replace:
```cpp
        ar & proper_pair_mapping_quality_;
        ar & improper_pair_mapping_quality_;
        ar & single_read_mapping_quality_;

        ar & gc_read_content_;
        ar & gc_read_content_reference_;
        ar & gc_read_content_mapped_;
        ar & n_content_;
        ar & sequence_content_;
        ar & sequence_content_reference_;

        ar & homopolymer_distribution_;
```

With:
```cpp
        ar & read_sequence_stats_.proper_pair_mapping_quality_;
        ar & read_sequence_stats_.improper_pair_mapping_quality_;
        ar & read_sequence_stats_.single_read_mapping_quality_;

        ar & read_sequence_stats_.gc_read_content_;
        ar & gc_read_content_reference_;
        ar & gc_read_content_mapped_;
        ar & read_sequence_stats_.n_content_;
        ar & read_sequence_stats_.sequence_content_;
        ar & sequence_content_reference_;

        ar & read_sequence_stats_.homopolymer_distribution_;
```

This preserves EXACT binary compatibility. The fields serialize in the same order, same types, just accessed through the sub-object. Note: `gc_read_content_reference_`, `gc_read_content_mapped_`, and `sequence_content_reference_` remain direct members of DataStats and serialize in their original positions.

**IMPORTANT:** Do NOT use `ar & read_sequence_stats_;` as a single call. That would change the archive format by wrapping the fields in a Boost.Serialization object header. Instead, serialize each sub-field individually at the exact same position it occupied before.

**1.4g: Update public getters.** Replace direct member access with forwarding to `read_sequence_stats_`:

```cpp
    const Vect<uintFragCount>& ProperPairMappingQuality() const {
        return read_sequence_stats_.ProperPairMappingQuality();
    }
    const Vect<uintFragCount>& ImproperPairMappingQuality() const {
        return read_sequence_stats_.ImproperPairMappingQuality();
    }
    const Vect<uintFragCount>& SingleReadMappingQuality() const {
        return read_sequence_stats_.SingleReadMappingQuality();
    }

    const Vect<uintFragCount>& GCReadContent(uintTempSeq template_segment) const {
        return read_sequence_stats_.GCReadContent(template_segment);
    }

    const Vect<uintFragCount>& NContent(uintTempSeq template_segment) const {
        return read_sequence_stats_.NContent(template_segment);
    }
    const Vect<uintNucCount>& SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const {
        return read_sequence_stats_.SequenceContent(template_segment, nucleotide);
    }

    const Vect<uintNucCount>& HomopolymerDistribution(uintBaseCall nucleotide) const {
        return read_sequence_stats_.HomopolymerDistribution(nucleotide);
    }
```

These getters that access data STAYING on DataStats are UNCHANGED:
- `GCReadContentReference()`, `GCReadContentMapped()`, `SequenceContentReference()`

- [ ] **Step 1.5: Modify `reseq/DataStats.cpp`**

**1.5a: Remove `EvalBaseLevelStats` implementation** (lines 153-222). This method is now on `ReadSequenceStats`.

**1.5b: Update `EvalRecord`** to call `read_sequence_stats_` instead of direct member access. Replace the two `EvalBaseLevelStats` calls (around lines 422-423):

Old:
```cpp
    EvalBaseLevelStats(record.first, (template_segment + 1) % 2, strand, tile_id, paired_seq_qual);
    EvalBaseLevelStats(record.second, template_segment, strand, tile_id, paired_seq_qual);
```

New:
```cpp
    read_sequence_stats_.EvalBaseLevelStats(record.first, (template_segment + 1) % 2, strand, tile_id,
                                            paired_seq_qual, qualities_, phred_quality_offset_, maximum_quality_);
    read_sequence_stats_.EvalBaseLevelStats(record.second, template_segment, strand, tile_id,
                                            paired_seq_qual, qualities_, phred_quality_offset_, maximum_quality_);
```

**1.5c: Update mapping quality increments in `EvalRecord`.** Replace all direct `++tmp_*_mapping_quality_.at(...)` with `read_sequence_stats_.IncrementMappingQuality(type, value)`:

Line ~428-429: `++tmp_single_read_mapping_quality_.at(record.first->record_.mapQ);` becomes `read_sequence_stats_.IncrementMappingQuality(2, record.first->record_.mapQ);`

Line ~430: `++tmp_single_read_mapping_quality_.at(record.second->record_.mapQ);` becomes `read_sequence_stats_.IncrementMappingQuality(2, record.second->record_.mapQ);`

Lines ~440-441: The improper pair block:
```cpp
            ++tmp_improper_pair_mapping_quality_.at(record.first->record_.mapQ);
            ++tmp_improper_pair_mapping_quality_.at(record.second->record_.mapQ);
```
becomes:
```cpp
            read_sequence_stats_.IncrementMappingQuality(1, record.first->record_.mapQ);
            read_sequence_stats_.IncrementMappingQuality(1, record.second->record_.mapQ);
```

Lines ~443-444: The proper pair block:
```cpp
            ++tmp_proper_pair_mapping_quality_.at(record.first->record_.mapQ);
            ++tmp_proper_pair_mapping_quality_.at(record.second->record_.mapQ);
```
becomes:
```cpp
            read_sequence_stats_.IncrementMappingQuality(0, record.first->record_.mapQ);
            read_sequence_stats_.IncrementMappingQuality(0, record.second->record_.mapQ);
```

**1.5d: Update `PrepareReadIn`.** Remove the accumulator preparation lines that moved to ReadSequenceStats (lines 651-680) and replace with a single call:

After the existing `qualities_.Prepare(...)` call (line 648), add:
```cpp
    read_sequence_stats_.PrepareAccumulators(size_mapping_quality, size_pos);
```

Remove the following lines from `PrepareReadIn` (they are now in `ReadSequenceStats::PrepareAccumulators`):
```cpp
    // Lines 651-653 (mapping quality resizes)
    tmp_proper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_improper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_single_read_mapping_quality_.resize(size_mapping_quality);

    // Lines 660-661 (gc/n content resizes within the template_segment loop)
        tmp_gc_read_content_.at(template_segment).resize(101);
        // ...
        tmp_n_content_.at(template_segment).resize(101);
    // Lines 665-666 (sequence content resizes)
        for (auto base = 5; base--;) {
            tmp_sequence_content_.at(template_segment).at(base).resize(size_pos);
        }

    // Lines 678-680 (homopolymer distribution resizes)
    for (auto base = 5; base--;) {
        tmp_homopolymer_distribution_.at(base).resize(size_pos);
    }
```

Keep these lines in `PrepareReadIn` (they relate to reference-derived histograms that stay on DataStats):
```cpp
        tmp_gc_read_content_reference_.at(template_segment).resize(101);
        tmp_gc_read_content_mapped_.at(template_segment).resize(101);

        for (auto strand = 2; strand--;) {
            for (auto base = 4; base--;) {
                tmp_sequence_content_reference_.at(template_segment)
                    .at(strand)
                    .at(base)
                    .resize(maximum_read_length_on_reference_ + 1);
            }
        }
```

Also keep the `tmp_read_lengths_by_fragment_length_` and `tmp_non_mapped_read_lengths_by_fragment_length_` resizes (lines 656-658) since those stay on DataStats.

**1.5e: Update `FinishReadIn`.** Remove the histogram finalization lines that moved to ReadSequenceStats and replace with a single call.

After the existing `qualities_.Finalize(...)` call (line 687), add:
```cpp
    read_sequence_stats_.Finalize();
```

Remove these lines from `FinishReadIn` (they are now in `ReadSequenceStats::Finalize`):
```cpp
    // Lines 690-692 (mapping quality Acquire)
    proper_pair_mapping_quality_.Acquire(tmp_proper_pair_mapping_quality_);
    improper_pair_mapping_quality_.Acquire(tmp_improper_pair_mapping_quality_);
    single_read_mapping_quality_.Acquire(tmp_single_read_mapping_quality_);

    // Within the template_segment loop:
    // Lines 710-711 (gc_read_content Acquire)
        gc_read_content_.at(template_segment).Acquire(tmp_gc_read_content_.at(template_segment));
    // Lines 714-715 (n_content Acquire)
        n_content_.at(template_segment).Acquire(tmp_n_content_.at(template_segment));
    // Lines 716-719 (sequence_content Acquire, the 5-base loop)
        for (auto base = 5; base--;) {
            sequence_content_.at(template_segment)
                .at(base)
                .Acquire(tmp_sequence_content_.at(template_segment).at(base));
        }

    // Lines 730-733 (homopolymer_distribution Acquire)
    tmp_homopolymer_distribution_.at(0).at(0) = 0;
    for (auto base = 5; base--;) {
        homopolymer_distribution_.at(base).Acquire(tmp_homopolymer_distribution_.at(base));
    }
```

Keep these lines in `FinishReadIn` (reference-derived histograms that stay on DataStats):
```cpp
        gc_read_content_reference_.at(template_segment).Acquire(tmp_gc_read_content_reference_.at(template_segment));
        gc_read_content_mapped_.at(template_segment).Acquire(tmp_gc_read_content_mapped_.at(template_segment));

        for (auto strand = 2; strand--;) {
            for (auto base = 4; base--;) {
                sequence_content_reference_.at(template_segment)
                    .at(strand)
                    .at(base)
                    .Acquire(tmp_sequence_content_reference_.at(template_segment).at(strand).at(base));
            }
        }
```

Also keep all `read_lengths_by_fragment_length_` and `non_mapped_read_lengths_by_fragment_length_` code since those stay.

**1.5f: The `Shrink()` method** does not touch the moved histograms, so no change is needed there.

- [ ] **Step 1.6: Update `reseq/DataStatsTest.cpp` member access paths**

All test code that accesses moved histogram members through `test_->` must now go through `test_->read_sequence_stats_`. This is a mechanical find-and-replace within the test file only.

**TestSequenceContent** (lines 59-73): Replace all `test_->sequence_content_` with `test_->read_sequence_stats_.sequence_content_`:
```cpp
    // e.g., line 59:
    // Old: test_->sequence_content_.at(template_segment).at(0)[at_pos]
    // New: test_->read_sequence_stats_.sequence_content_.at(template_segment).at(0)[at_pos]
```

Do the same for all 5 `EXPECT_EQ` lines in this method (lines 59, 62, 65, 68, 71).

**TestSrr490124Equality** (lines 93-279): Replace:
- `test_->proper_pair_mapping_quality_` with `test_->read_sequence_stats_.proper_pair_mapping_quality_` (lines 119-141)
- `test_->improper_pair_mapping_quality_` with `test_->read_sequence_stats_.improper_pair_mapping_quality_` (lines 121-128, 142-143)
- `test_->single_read_mapping_quality_` with `test_->read_sequence_stats_.single_read_mapping_quality_` (lines 129-147)
- `test_->gc_read_content_` with `test_->read_sequence_stats_.gc_read_content_` (lines 152-168)
- `test_->n_content_` with `test_->read_sequence_stats_.n_content_` (lines 199-203)
- `test_->homopolymer_distribution_` with `test_->read_sequence_stats_.homopolymer_distribution_` (lines 269-278)

These stay unchanged (reference-derived histograms remain on DataStats):
- `test_->gc_read_content_reference_` (lines 173-184)
- `test_->gc_read_content_mapped_` (lines 187-198)
- `test_->sequence_content_reference_` (lines 234-264 and TestSequenceContentReference method)

**TestTiles** (lines 281-295): Replace:
- `test_->sequence_content_` with `test_->read_sequence_stats_.sequence_content_` (lines 286-289)
- `test_->homopolymer_distribution_` with `test_->read_sequence_stats_.homopolymer_distribution_` (line 293)

**TestDuplicates** (lines 297-334): Replace:
- `test_->proper_pair_mapping_quality_` with `test_->read_sequence_stats_.proper_pair_mapping_quality_` (lines 308-314)
- `test_->improper_pair_mapping_quality_` with `test_->read_sequence_stats_.improper_pair_mapping_quality_` (lines 315-318)
- `test_->single_read_mapping_quality_` with `test_->read_sequence_stats_.single_read_mapping_quality_` (lines 319-320)

**TestVariants** (lines 336-381): Reference-derived histograms stay unchanged. No moved members are accessed in this method.

**TestCrossDuplicates** (lines 383-408): Replace:
- `test_->proper_pair_mapping_quality_` with `test_->read_sequence_stats_.proper_pair_mapping_quality_` (lines 397-398)
- `test_->improper_pair_mapping_quality_` with `test_->read_sequence_stats_.improper_pair_mapping_quality_` (lines 400-405)
- `test_->single_read_mapping_quality_` with `test_->read_sequence_stats_.single_read_mapping_quality_` (lines 406-407)

**Construction test** (lines 487-504): Replace:
- `test_->sequence_content_` with `test_->read_sequence_stats_.sequence_content_` (lines 496-499)
- `test_->Shrink()` call on line 503 stays as-is (Shrink is still on DataStats)

**Regex for bulk replacement:**
```
test_->sequence_content_\.  ->  test_->read_sequence_stats_.sequence_content_.
test_->proper_pair_mapping_quality_  ->  test_->read_sequence_stats_.proper_pair_mapping_quality_
test_->improper_pair_mapping_quality_  ->  test_->read_sequence_stats_.improper_pair_mapping_quality_
test_->single_read_mapping_quality_  ->  test_->read_sequence_stats_.single_read_mapping_quality_
test_->gc_read_content_\.  ->  test_->read_sequence_stats_.gc_read_content_.
test_->n_content_\.  ->  test_->read_sequence_stats_.n_content_.
test_->homopolymer_distribution_\.  ->  test_->read_sequence_stats_.homopolymer_distribution_.
```

**CRITICAL:** Do NOT replace `test_->sequence_content_reference_`, `test_->gc_read_content_reference_`, or `test_->gc_read_content_mapped_`. These stay on DataStats. The regex uses `\.` after the member name (before `.at(`) which helps distinguish `sequence_content_.` from `sequence_content_reference_.`.

- [ ] **Step 1.7: Build and test**

```bash
make build && make test
```

All tests must pass. If compilation fails, check:
1. Missing includes in ReadSequenceStats.h/cpp
2. Missed member removals in DataStats.h
3. Incomplete access path updates in DataStatsTest.cpp

- [ ] **Step 1.8: Format and commit**

```bash
make format
git add reseq/ReadSequenceStats.h reseq/ReadSequenceStats.cpp reseq/DataStats.h reseq/DataStats.cpp reseq/DataStatsTest.cpp reseq/CMakeLists.txt
git commit -m "refactor(5b): extract ReadSequenceStats from DataStats

Move read-level histogram data and EvalBaseLevelStats into a new
ReadSequenceStats class. DataStats forwards getters and delegates
accumulation. Serialization field order preserved for binary
compatibility. Test access paths updated mechanically."
```

---

## Task 2: Extract BamIngestionEngine

**Goal:** Create `BamIngestionEngine` class owning the two-pass BAM reading pipeline, threading, pair-matching, record evaluation, and validation predicates. DataStats::ReadBam becomes a thin wrapper.

**Files created:**
- `reseq/BamIngestionEngine.h`
- `reseq/BamIngestionEngine.cpp`

**Files modified:**
- `reseq/DataStats.h`
- `reseq/DataStats.cpp`
- `reseq/CMakeLists.txt`

### Step-by-step

- [ ] **Step 2.1: Create `reseq/BamIngestionEngine.h`**

```cpp
#ifndef BAMINGESTIONENGINE_H
#define BAMINGESTIONENGINE_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <unordered_set>
#include <utility>
#include <vector>

#include <seqan/bam_io.h>

#include "CoverageStats.h"
#include "utilities.hpp"

namespace reseq {

// Forward declaration
class DataStats;

class BamIngestionEngine {
  public:
    // Hash/equality functors for FullRecord pair-matching (moved from DataStats)
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

    // Static geometry utilities (unchanged signatures)
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record);
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record,
                                                uintReadLen& max_indel);

  private:
    struct ThreadData {
        std::vector<std::pair<CoverageStats::FullRecord*, CoverageStats::FullRecord*>> rec_store_;

        CoverageStats::ThreadData coverage_;
        FragmentDistributionStats::ThreadData fragment_distribution_;

        uintSeqLen last_exclusion_region_id_;
        uintRefSeqId last_exclusion_ref_seq_;

        ThreadData(uintSeqLen maximum_insert_length, uintSeqLen start_exclusion_length)
            : fragment_distribution_(maximum_insert_length, start_exclusion_length),
              last_exclusion_region_id_(0), last_exclusion_ref_seq_(0) {}
    };

    // Definitions
    const uintFragCount kBatchSize = 10000;

    // Threading primitives
    std::mutex read_mutex_;
    std::mutex print_mutex_;
    std::atomic<uintNumThreads> running_threads_;
    std::mutex finish_threads_mutex_;
    std::condition_variable finish_threads_cv_;
    bool finish_threads_;

    // Pair-matching state
    std::unordered_set<CoverageStats::FullRecord*, RecordHasher, RecordEqual> first_read_records_;
    std::atomic<bool> reading_success_;
    uintFragCount read_records_;

    // Pre-run temporaries
    std::vector<uintFragCount> reads_per_frag_len_bin_;
    std::vector<uintFragCount> lowq_reads_per_frag_len_bin_;

    // Read-disposition counters
    std::atomic<uintFragCount> reads_in_unmapped_pairs_without_adapters_;
    std::atomic<uintFragCount> reads_in_unmapped_pairs_with_adapters_;
    std::atomic<uintFragCount> reads_with_low_quality_with_adapters_;
    std::atomic<uintFragCount> reads_with_low_quality_without_adapters_;
    std::atomic<uintFragCount> reads_on_too_short_fragments_;
    std::atomic<uintFragCount> reads_in_excluded_regions_;
    std::atomic<uintFragCount> reads_used_;

    // Private methods
    inline bool PotentiallyValidGeneral(const seqan::BamAlignmentRecord& record,
                                        const DataStats& target) const;
    inline bool PotentiallyValidFirst(const seqan::BamAlignmentRecord& record_first,
                                      const DataStats& target) const;
    inline bool PotentiallyValidSecond(const seqan::BamAlignmentRecord& record_second,
                                       const DataStats& target) const;
    inline bool PotentiallyValid(const seqan::BamAlignmentRecord& record,
                                 const DataStats& target) const;

    bool IsSecondRead(CoverageStats::FullRecord* record, CoverageStats::FullRecord*& record_first,
                      CoverageStats::CoverageBlock*& block, DataStats& target);

    bool CheckForAdapters(const seqan::BamAlignmentRecord& record_first,
                          const seqan::BamAlignmentRecord& record_second, DataStats& target);
    bool EvalReferenceStatistics(CoverageStats::FullRecord* record, uintTempSeq template_segment,
                                 CoverageStats::CoverageBlock* coverage_block, DataStats& target);
    bool EvalRecord(std::pair<CoverageStats::FullRecord*, CoverageStats::FullRecord*> record,
                    ThreadData& thread, DataStats& target);

    bool SignsOfPairsWithNamesNotIdentical(const DataStats& target);
    void PrepareReadIn(uintQual size_mapping_quality, uintReadLen size_indel,
                       uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads, DataStats& target);
    void FinishReadIn(DataStats& target);
    void Shrink(DataStats& target);
    bool Calculate(uintNumThreads num_threads, DataStats& target);
    void PrepareGeneral(DataStats& target);

    bool OrderOfBamFileCorrect(const seqan::BamAlignmentRecord& record,
                               std::pair<uintRefSeqId, uintSeqLen> last_record_pos);
    bool PreRun(seqan::BamFileIn& bam, const char* bam_file, seqan::BamHeader& header,
                uintQual& size_mapping_quality, uintReadLen& size_indel, DataStats& target);
    bool ReadRecords(seqan::BamFileIn& bam, bool& not_done, ThreadData& thread_data, DataStats& target);

    static void ReadThread(BamIngestionEngine& engine, seqan::BamFileIn& bam, size_t thread_idx,
                           DataStats& target);

    // Google test — DataStatsTest checks first_read_records_.size() == 0
    friend class DataStatsTest;
};

} // namespace reseq

#endif // BAMINGESTIONENGINE_H
```

Key points:
- All private methods receive `DataStats& target` to access sub-stats (coverage, errors, etc.) and remaining DataStats data (reference, config scalars, reference-derived histograms).
- `RecordHasher`/`RecordEqual` move here verbatim from DataStats.
- `ThreadData` struct moves here verbatim from DataStats.
- `GetReadLengthOnReference` (static) moves here — these are pure geometry utilities that only read from the BAM record.
- `friend class DataStatsTest` is needed because tests check `first_read_records_.size()` directly.

- [ ] **Step 2.2: Create `reseq/BamIngestionEngine.cpp`**

This is the largest file. Move ALL of the following method implementations from DataStats.cpp:

**Methods to move (with modifications):**

1. **RecordHasher::operator()** (lines 73-83) — Change `DataStats::RecordHasher` to `BamIngestionEngine::RecordHasher`
2. **RecordEqual::operator()** (lines 85-88) — Change `DataStats::RecordEqual` to `BamIngestionEngine::RecordEqual`
3. **PotentiallyValidGeneral** (lines 90-93) — Add `DataStats& target` parameter. Change `reference_->` to `target.reference_->` and `maximum_insert_length_` references.
4. **PotentiallyValidFirst** (lines 94-97) — Same pattern.
5. **PotentiallyValidSecond** (lines 98-101) — Same pattern.
6. **PotentiallyValid** (lines 102-108) — Same pattern.
7. **IsSecondRead** (lines 110-133) — Add `DataStats& target` parameter. Change `coverage_` to `target.coverage_` and `maximum_read_length_on_reference_` to `target.maximum_read_length_on_reference_`.
8. **CheckForAdapters** (lines 135-151) — Add `DataStats& target`. Change `adapters_` to `target.adapters_`, `fragment_distribution_` to `target.fragment_distribution_`, `tmp_non_mapped_read_lengths_by_fragment_length_` to `target.tmp_non_mapped_read_lengths_by_fragment_length_`.
9. **EvalReferenceStatistics** (lines 224-393) — Add `DataStats& target`. Change `coverage_` to `target.coverage_`, `reference_` to `target.reference_`, `errors_` to `target.errors_`, `phred_quality_offset_` to `target.phred_quality_offset_`, `tmp_sequence_content_reference_` to `target.tmp_sequence_content_reference_`, `tmp_gc_read_content_reference_` to `target.tmp_gc_read_content_reference_`, `tmp_gc_read_content_mapped_` to `target.tmp_gc_read_content_mapped_`.
10. **EvalRecord** (lines 395-606) — Add `DataStats& target`. Change all member accesses: `tiles_` to `target.tiles_`, `reference_` to `target.reference_`, `read_sequence_stats_` to `target.read_sequence_stats_`, `qualities_` to `target.qualities_`, `phred_quality_offset_` to `target.phred_quality_offset_`, `maximum_quality_` to `target.maximum_quality_`, `minimum_mapping_quality_` to `target.minimum_mapping_quality_`, `maximum_insert_length_` to `target.maximum_insert_length_`, `maximum_read_length_on_reference_` to `target.maximum_read_length_on_reference_`, `coverage_` to `target.coverage_`, `fragment_distribution_` to `target.fragment_distribution_`, `adapters_` to `target.adapters_`, `tmp_read_lengths_by_fragment_length_` to `target.tmp_read_lengths_by_fragment_length_`. Also replace `QualitySufficient(...)` calls with `target.QualitySufficient(...)`, `GetReadPosOnReference(...)` calls with `target.GetReadPosOnReference(...)`, `InProperDirection(...)` calls with `target.InProperDirection(...)`.
11. **SignsOfPairsWithNamesNotIdentical** (lines 608-621) — Add `const DataStats& target`. Change `reference_` to `target.reference_`.
12. **PrepareReadIn** (lines 623-681) — Add `DataStats& target`. All sub-stat Prepare calls go through `target.`. The `read_sequence_stats_.PrepareAccumulators(...)` call also goes through `target.`.
13. **FinishReadIn** (lines 683-734) — Add `DataStats& target`. All sub-stat Finalize calls and histogram Acquire calls go through `target.`.
14. **Shrink** (lines 736-746) — Add `DataStats& target`. All sub-stat Shrink calls go through `target.`.
15. **Calculate** (lines 749-771) — Add `DataStats& target`.
16. **PrepareGeneral** (lines 773-778) — Add `DataStats& target`.
17. **OrderOfBamFileCorrect** (lines 780-793) — No target needed (pure record validation).
18. **PreRun** (lines 795-934) — Add `DataStats& target`. All member accesses go through `target.`.
19. **ReadRecords** (lines 937-1018) — Add `DataStats& target`.
20. **ReadThread** (lines 1021-1121) — Change parameter from `DataStats& self` to `BamIngestionEngine& engine, ... DataStats& target`. Replace `self.` with `engine.` for engine members and `target.` for DataStats members.
21. **GetReadLengthOnReference** (both overloads, lines 1162-1198) — Move verbatim, change `DataStats::` to `BamIngestionEngine::`.
22. **IsValidRecord** (lines 1139-1160) — This stays on DataStats (it's a public method). Do NOT move.

**The Run method** implements the body of `DataStats::ReadBam` (lines 1237-1411):

```cpp
bool BamIngestionEngine::Run(const char* bam_file, const char* adapter_file,
                             const char* adapter_matrix, const std::string& variant_file,
                             uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
                             bool calculate_bias, DataStats& target) {
    // Initialize engine state
    reading_success_ = true;
    read_records_ = 0;
    reads_in_unmapped_pairs_without_adapters_ = 0;
    reads_in_unmapped_pairs_with_adapters_ = 0;
    reads_with_low_quality_with_adapters_ = 0;
    reads_with_low_quality_without_adapters_ = 0;
    reads_on_too_short_fragments_ = 0;
    reads_in_excluded_regions_ = 0;
    reads_used_ = 0;

    // ... rest is the body of DataStats::ReadBam (lines 1240-1411)
    // with all member accesses changed to target.member_
    // and self-references changed to engine methods
}
```

The method body follows the exact structure of `DataStats::ReadBam` but with:
- `reference_` becomes `target.reference_`
- `fragment_distribution_` becomes `target.fragment_distribution_`
- `adapters_` becomes `target.adapters_`
- `coverage_` becomes `target.coverage_`
- `total_number_reads_` becomes `target.total_number_reads_`
- Sub-method calls like `PreRun(...)` become `PreRun(..., target)`
- `ReadThread(*this, bam, i)` becomes `ReadThread(*this, bam, i, target)`
- `PrepareReadIn(...)` becomes `PrepareReadIn(..., target)`
- `FinishReadIn()` becomes `FinishReadIn(target)`
- `Shrink()` becomes `Shrink(target)`
- `Calculate(num_threads)` becomes `Calculate(num_threads, target)`

The final printing block (lines 1357-1384) uses the engine's own counter members, not `target.`.

**Includes for BamIngestionEngine.cpp:**
```cpp
#include "BamIngestionEngine.h"
#include "DataStats.h"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <thread>
#include <vector>

#include "logging.hpp"
#include "CMakeConfig.h"

#include <seqan/bam_io.h>

// Using declarations matching those removed from DataStats.cpp
using reseq::BamIngestionEngine;
using reseq::DataStats;
using std::lock_guard;
using std::mutex;
using std::pair;
using std::string;
using std::thread;
using std::unique_lock;
using std::unordered_set;
using std::vector;
using seqan::atEnd;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamHeader;
using seqan::contigNames;
using seqan::Dna;
using seqan::Dna5;
using seqan::Exception;
using seqan::readHeader;
using seqan::readRecord;

using reseq::utilities::at;
using reseq::utilities::Complement;
using reseq::utilities::ConstIupacStringReverseComplement;
using reseq::utilities::Divide;
using reseq::utilities::MeanWithRoundingToFirst;
using reseq::utilities::Percent;
using reseq::utilities::ReversedConstCharString;
using reseq::utilities::ReversedConstCigarString;
using reseq::utilities::SafePercent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;
```

- [ ] **Step 2.3: Add BamIngestionEngine.cpp to CMakeLists.txt**

In `reseq/CMakeLists.txt`, add `BamIngestionEngine.cpp` to the `reseq_lib` sources list (before `CoverageStats.cpp`):

```cmake
add_library(reseq_lib STATIC
  AdapterStats.cpp
  archive_format.cpp
  BamIngestionEngine.cpp
  CoverageStats.cpp
  DataStats.cpp
  DataStatsInterface.cpp
  ErrorStats.cpp
  FragmentDistributionStats.cpp
  FragmentDuplicationStats.cpp
  ProbabilityEstimates.cpp
  QualityStats.cpp
  ReadSequenceStats.cpp
  Reference.cpp
  Simulator.cpp
  Surrounding.cpp
  TileStats.cpp
)
```

- [ ] **Step 2.4: Modify `reseq/DataStats.h`**

**2.4a: Add include:**
```cpp
#include "BamIngestionEngine.h"
```

**2.4b: Remove moved types.** Delete the `RecordHasher` and `RecordEqual` struct declarations from the public section (lines 40-45). Delete the `ThreadData` struct from the private section (lines 48-60).

**2.4c: Remove moved data members.** Delete from private section:
- `kBatchSize` (line 63)
- `read_mutex_`, `print_mutex_` (lines 73-74)
- `running_threads_` (line 76)
- `finish_threads_mutex_` (line 77)
- `finish_threads_cv_` (line 78)
- `finish_threads_` (line 79)
- `first_read_records_` (line 91)
- `reading_success_` (line 92)
- `read_records_` (line 93)
- `reads_per_frag_len_bin_` (lines 115-116)
- `lowq_reads_per_frag_len_bin_` (lines 117-119)
- `reads_in_unmapped_pairs_without_adapters_` through `reads_used_` (lines 165-172)

**2.4d: Remove moved private method declarations.** Delete:
- `PotentiallyValidGeneral`, `PotentiallyValidFirst`, `PotentiallyValidSecond`, `PotentiallyValid` (lines 174-177)
- `IsSecondRead` (lines 178-179)
- `CheckForAdapters` (lines 181-182)
- `EvalReferenceStatistics` (lines 185-186)
- `EvalRecord` (lines 187-188)
- `SignsOfPairsWithNamesNotIdentical` (line 190)
- `PrepareReadIn` (lines 191-192)
- `FinishReadIn` (line 193)
- `Shrink` (line 194)
- `Calculate` (line 195)
- `PrepareGeneral` (line 196)
- `OrderOfBamFileCorrect` (lines 198-199)
- `PreRun` (lines 200-201)
- `ReadRecords` (line 202)
- `ReadThread` (line 204)

**2.4e: Remove moved public methods.** Delete from public section:
- `IsValidRecord` declaration (line 314) — Wait, this stays. `IsValidRecord` is public and called externally. Keep it.
- `QualitySufficient` (line 315-317) — This stays. It's public and called in Simulator.
- `GetReadLengthOnReference` (both overloads, lines 318-319) — These move to BamIngestionEngine. But they are also called externally by Simulator. Check external usage first.

**IMPORTANT CHECK:** Before moving `GetReadLengthOnReference`, verify if it's called outside DataStats:

```bash
grep -rn "GetReadLengthOnReference" reseq/ --include="*.cpp" --include="*.h" | grep -v DataStats | grep -v BamIngestionEngine
```

If it is called from Simulator or elsewhere, it must ALSO remain as a static method on DataStats (forwarding to BamIngestionEngine::GetReadLengthOnReference) OR stay on DataStats. Per the design spec, it can stay on DataStats as a static method while also being on BamIngestionEngine. The simplest approach: **keep the static `GetReadLengthOnReference` declarations on DataStats** and just move the implementation to BamIngestionEngine, with DataStats forwarding:

```cpp
// In DataStats.h (keep existing declarations)
static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record);
static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record, uintReadLen& max_indel);

// In DataStats.cpp (forwarding implementations)
reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord& record) {
    return BamIngestionEngine::GetReadLengthOnReference(record);
}
reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord& record, uintReadLen& max_indel) {
    return BamIngestionEngine::GetReadLengthOnReference(record, max_indel);
}
```

Similarly for `InProperDirection` — keep it on DataStats (it's public and inline, used by external callers).

Similarly for `IsValidRecord` and `QualitySufficient` — keep on DataStats.

**2.4f: Keep `GetReadPosOnReference`** on DataStats (non-static, depends on `reference_`).

- [ ] **Step 2.5: Modify `reseq/DataStats.cpp`**

**2.5a: Add include:**
```cpp
#include "BamIngestionEngine.h"
```

**2.5b: Remove all moved method implementations.** Remove from DataStats.cpp:
- `RecordHasher::operator()` (lines 73-83)
- `RecordEqual::operator()` (lines 85-88)
- `PotentiallyValidGeneral` through `PotentiallyValid` (lines 90-108)
- `IsSecondRead` (lines 110-133)
- `CheckForAdapters` (lines 135-151)
- `EvalReferenceStatistics` (lines 224-393)
- `EvalRecord` (lines 395-606)
- `SignsOfPairsWithNamesNotIdentical` (lines 608-621)
- `PrepareReadIn` (lines 623-681)
- `FinishReadIn` (lines 683-734)
- `Shrink` (lines 736-746)
- `Calculate` (lines 749-771)
- `PrepareGeneral` (lines 773-778)
- `OrderOfBamFileCorrect` (lines 780-793)
- `PreRun` (lines 795-934)
- `ReadRecords` (lines 937-1018)
- `ReadThread` (lines 1021-1121)
- `GetReadLengthOnReference` (both overloads, lines 1162-1198) — replace with forwarding stubs

**2.5c: Replace `ReadBam` body** (lines 1237-1411) with thin wrapper:

```cpp
bool DataStats::ReadBam(const char* bam_file, const char* adapter_file, const char* adapter_matrix,
                        const string& variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
                        bool calculate_bias) {
    BamIngestionEngine engine;
    return engine.Run(bam_file, adapter_file, adapter_matrix, variant_file, max_ref_seq_bin_size, num_threads,
                      calculate_bias, *this);
}
```

**2.5d: Update the constructor initializer list.** Remove initialization of moved members:
- `reading_success_(true)`, `read_records_(0)`
- `reads_in_unmapped_pairs_without_adapters_(0)` through `reads_used_(0)`

The constructor initialization (lines 1123-1137) becomes significantly shorter. Keep:
- `reference_(ref)`, `maximum_insert_length_(maximum_insert_length)`, `minimum_mapping_quality_(minimum_mapping_quality)`
- `minimum_quality_(255)`, `maximum_quality_(0)`
- `minimum_read_length_on_reference_(numeric_limits<uintReadLen>::max())`, `maximum_read_length_on_reference_(0)`
- `total_number_reads_(0)`
- The `read_lengths_` offset setup in the constructor body stays.

**2.5e: Remove unused includes from DataStats.cpp.** After the extraction, these are likely no longer needed in DataStats.cpp:
- `<condition_variable>` (threading moved)
- `<thread>` (moved)
- `<unordered_set>` (moved)

Verify by building. Keep any include that is still transitively needed.

**2.5f: Remove unused using declarations** from DataStats.cpp that correspond to moved code (e.g., `using std::thread;`, `using std::unordered_set;`, `using std::lock_guard;`).

- [ ] **Step 2.6: Update `reseq/DataStatsTest.cpp` for `first_read_records_`**

The test file accesses `test_->first_read_records_` in four places (lines 102, 300, 393, 412). Since `first_read_records_` moved to BamIngestionEngine and BamIngestionEngine is created on the stack and destroyed after Run() completes, these assertions always check `== 0` which means "all pairs were matched". After the engine is destroyed, the set doesn't exist anymore.

The cleanest approach: since these assertions only verify that `first_read_records_` is empty after ReadBam completes (which is guaranteed by `SignsOfPairsWithNamesNotIdentical` returning false), we can simply remove these 4 assertions. They were testing an internal invariant that is now encapsulated within BamIngestionEngine.

Alternatively, if we want to preserve them, BamIngestionEngine can store a `bool pairs_matched_` flag set at the end of Run(), and DataStats can expose it. But since SignsOfPairsWithNamesNotIdentical already validates this and ReadBam would have returned false if it failed, removing the assertions is safe.

**Decision:** Remove the 4 `first_read_records_` assertions from DataStatsTest.cpp:
- Line 101-102: Remove `EXPECT_EQ(0, test_->first_read_records_.size()) ...`
- Line 299-300: Remove
- Line 392-393: Remove
- Line 412: Remove

- [ ] **Step 2.7: Build and test**

```bash
make build && make test
```

This is the most complex step. Likely issues:
1. Missing `target.` prefixes on member accesses in BamIngestionEngine.cpp
2. Friend access issues — ensure `friend class BamIngestionEngine` is in DataStats.h
3. Circular include issues — BamIngestionEngine.h forward-declares DataStats, BamIngestionEngine.cpp includes DataStats.h
4. Missing `using` declarations in BamIngestionEngine.cpp

Debug strategy: If compile errors are numerous, focus on one method at a time. Comment out methods in BamIngestionEngine.cpp and add them back incrementally.

- [ ] **Step 2.8: Format and commit**

```bash
make format
git add reseq/BamIngestionEngine.h reseq/BamIngestionEngine.cpp reseq/DataStats.h reseq/DataStats.cpp reseq/DataStatsTest.cpp reseq/CMakeLists.txt
git commit -m "refactor(5b): extract BamIngestionEngine from DataStats

Move two-pass BAM reading pipeline, threading primitives, pair-matching,
record evaluation, and validation predicates into BamIngestionEngine.
DataStats::ReadBam becomes a thin wrapper that constructs an engine on
the stack. All external callers remain unchanged."
```

---

## Task 3: Clean Up and Final Verification

**Goal:** Remove dead includes, run formatting, verify all tests, run pre-commit hooks.

**Files modified:**
- `reseq/DataStats.h`
- `reseq/DataStats.cpp`

- [ ] **Step 3.1: Clean up includes in `reseq/DataStats.h`**

After extraction, these includes may no longer be needed in DataStats.h:
- `<condition_variable>` — threading moved to BamIngestionEngine
- `<mutex>` — mutexes moved to BamIngestionEngine
- `<unordered_set>` — `first_read_records_` moved to BamIngestionEngine
- `<seqan/bam_io.h>` — CHECK: still needed for `GetReadPosOnReference`, `IsValidRecord`, `QualitySufficient`, `InProperDirection` signatures. KEEP.

Do NOT remove:
- `<array>`, `<atomic>`, `<vector>` — still used for remaining members
- `<chrono>` — CHECK: only used in ReadBam for `creation_time_`. This moved to engine. But `creation_time_` is still a DataStats member set by the engine. Remove `<chrono>` from DataStats.h only if it's not needed.
- `<map>` — CHECK: may be needed by sub-stat includes
- `<string>` — still needed
- Boost serialization headers — still needed
- Sub-stat includes — still needed
- `SeqQualityStats.hpp` — CHECK: was needed for `EvalBaseLevelStats` which moved. Can likely be removed from DataStats.h.
- `"utilities.hpp"`, `"Vect.hpp"` — still needed

Build and verify after each include removal.

- [ ] **Step 3.2: Clean up includes in `reseq/DataStats.cpp`**

Remove using declarations for types/functions that only the moved code used:
- Check each `using` declaration against remaining code
- Remove any `using seqan::...` that are only used in moved methods

- [ ] **Step 3.3: Verify DataStats.h is clean and consistent**

Audit the remaining DataStats class:
- All remaining data members are either config, sub-stats, metadata, read-length histograms, or reference-derived histograms
- `read_sequence_stats_` is the only new member
- All public getters work correctly (forwarding or direct)
- `ReadBam` is a thin wrapper
- `serialize()` preserves exact field order
- Friend declarations include `BamIngestionEngine` and `DataStatsTest`

- [ ] **Step 3.4: Run full verification suite**

```bash
make format-check
make build
make test
pre-commit run --all-files
```

All four must pass.

- [ ] **Step 3.5: Final line count audit**

```bash
wc -l reseq/DataStats.h reseq/DataStats.cpp reseq/ReadSequenceStats.h reseq/ReadSequenceStats.cpp reseq/BamIngestionEngine.h reseq/BamIngestionEngine.cpp
```

Expected approximate sizes:
- DataStats.h: ~150-180 lines (down from 347)
- DataStats.cpp: ~200-300 lines (down from 1497)
- ReadSequenceStats.h: ~100 lines
- ReadSequenceStats.cpp: ~150 lines
- BamIngestionEngine.h: ~130 lines
- BamIngestionEngine.cpp: ~1000-1100 lines

Total lines should be roughly the same as the original (conservation of code).

- [ ] **Step 3.6: Commit cleanup**

```bash
make format
git add reseq/DataStats.h reseq/DataStats.cpp
git commit -m "refactor(5b): clean up DataStats after decomposition

Remove dead includes and using declarations from DataStats.h/.cpp
after ReadSequenceStats and BamIngestionEngine extraction."
```

---

## Verification Checklist

After all three tasks are complete, verify:

- [ ] `make build` succeeds with no warnings related to the changed files
- [ ] `make test` passes all tests (DataStatsTest, ProbabilityEstimatesTest, SimulatorTest, RegressionTest, all others)
- [ ] `make format-check` is clean
- [ ] `pre-commit run --all-files` passes
- [ ] No public API changes — grep confirms external callers unchanged:
  ```bash
  grep -rn "DataStats::" reseq/ --include="*.cpp" --include="*.h" | grep -v "DataStats\.\|DataStatsTest\.\|BamIngestionEngine\.\|ReadSequenceStats\." | head -30
  ```
- [ ] Binary serialization compatibility — Load/Save round-trip works (tested by DataStatsTest::Ecoli which does save/load/verify)
- [ ] Text serialization round-trip works (tested by DataStatsTest::Ecoli text format tests)

---

## Risk Mitigations

| Risk | Mitigation |
|------|------------|
| Serialization format break | Serialize `read_sequence_stats_` fields individually in exact original order within DataStats::serialize(), not as a sub-object |
| Circular includes | BamIngestionEngine.h forward-declares DataStats; BamIngestionEngine.cpp includes DataStats.h |
| Friend access insufficient | BamIngestionEngine is friend of DataStats; DataStatsTest is friend of both DataStats and ReadSequenceStats |
| Missed member access in BamIngestionEngine | Build will fail immediately — all accesses must go through `target.` |
| GetReadLengthOnReference external callers | Keep static forwarding stubs on DataStats that delegate to BamIngestionEngine |
| first_read_records_ test assertions | Remove (invariant now encapsulated in BamIngestionEngine) |
