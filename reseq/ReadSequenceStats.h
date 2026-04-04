#ifndef READSEQUENCESTATS_H
#define READSEQUENCESTATS_H

#include <array>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include "container_types.hpp"
#include "CoverageStats.h"
#include "QualityStats.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq {

class ReadSequenceStats {
  private:
    // Atomic accumulators (thread-safe, filled during multi-threaded BAM pass)
    AtomicVec<uintFragCount> tmp_proper_pair_mapping_quality_;
    AtomicVec<uintFragCount> tmp_improper_pair_mapping_quality_;
    AtomicVec<uintFragCount> tmp_single_read_mapping_quality_;

    PerSegment<AtomicVec<uintFragCount>> tmp_gc_read_content_;
    PerSegment<AtomicVec<uintFragCount>> tmp_n_content_;
    PerSegmentPerBaseN<AtomicVec<uintNucCount>> tmp_sequence_content_;
    PerBaseN<AtomicVec<uintNucCount>> tmp_homopolymer_distribution_;

    // Final histograms (serialized)
    Vect<uintFragCount> proper_pair_mapping_quality_;
    Vect<uintFragCount> improper_pair_mapping_quality_;
    Vect<uintFragCount> single_read_mapping_quality_;

    PerSegment<Vect<uintFragCount>> gc_read_content_;
    PerSegment<Vect<uintFragCount>> n_content_;
    PerSegmentPerBaseN<Vect<uintNucCount>> sequence_content_;
    PerBaseN<Vect<uintNucCount>> homopolymer_distribution_;

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

    // DataStats needs direct access for inline serialization (binary compatibility)
    friend class DataStats;

    // Google test
    friend class DataStatsTest;

  public:
    // Accumulator lifecycle
    void PrepareAccumulators(uintQual size_mapping_quality, uintReadLen size_pos);
    void EvalBaseLevelStats(CoverageStats::FullRecord* full_record, uintTempSeq template_segment, uintTempSeq strand,
                            uintTileId tile_id, uintQual& paired_seq_qual, QualityStats& qualities,
                            uintQual phred_quality_offset, uintQual maximum_quality);
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
    const Vect<uintFragCount>& NContent(uintTempSeq template_segment) const { return n_content_.at(template_segment); }
    const Vect<uintNucCount>& SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const {
        return sequence_content_.at(template_segment).at(nucleotide);
    }
    const Vect<uintNucCount>& HomopolymerDistribution(uintBaseCall nucleotide) const {
        return homopolymer_distribution_.at(nucleotide);
    }
};

} // namespace reseq

#endif // READSEQUENCESTATS_H
