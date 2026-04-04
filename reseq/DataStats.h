#ifndef DATASTATS_H
#define DATASTATS_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <seqan/bam_io.h>

#include "AdapterStats.h"
#include "container_types.hpp"
#include "CoverageStats.h"
#include "ErrorStats.h"
#include "FragmentDistributionStats.h"
#include "FragmentDuplicationStats.h"
#include "QualityStats.h"
#include "ReadSequenceStats.h"
#include "Reference.h"
#include "TileStats.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq {

class BamIngestionEngine; // Forward declaration (full header in DataStats.cpp)

class DataStats {
  private:
    // User parameter
    Reference* reference_;

    const uintSeqLen maximum_insert_length_;
    const uintQual minimum_mapping_quality_; // Minimum mapping quality required for both reads to use read pair for
                                             // insert_lengths_, fragment_duplication_number_, ...

    // Subclasses
    AdapterStats adapters_;
    CoverageStats coverage_;
    FragmentDuplicationStats duplicates_;
    ErrorStats errors_;
    FragmentDistributionStats fragment_distribution_;
    QualityStats qualities_;
    TileStats tiles_;
    ReadSequenceStats read_sequence_stats_;

    // Temporary variables
    PerSegment<std::vector<AtomicVec<uintFragCount>>> tmp_read_lengths_by_fragment_length_;
    PerSegment<std::vector<AtomicVec<uintFragCount>>> tmp_non_mapped_read_lengths_by_fragment_length_;

    PerSegment<AtomicVec<uintFragCount>> tmp_gc_read_content_reference_;
    PerSegment<AtomicVec<uintFragCount>> tmp_gc_read_content_mapped_;

    PerSegmentPerStrand<PerBase<AtomicVec<uintNucCount>>> tmp_sequence_content_reference_;

    // Collected variables for simulation
    uint64_t creation_time_; // Store time when bam file was completelly read, can be used to check whether the stats
                             // file was updated
    PerSegment<Vect<uintFragCount>> read_lengths_; // read_lengths_[first/second][length] = #reads
    PerSegment<Vect<Vect<uintFragCount>>>
        read_lengths_by_fragment_length_; // read_lengths_by_fragment_length_[first/second][fragment_length][read_length]
                                          // = #reads
    PerSegment<Vect<Vect<uintFragCount>>>
        non_mapped_read_lengths_by_fragment_length_; // non_mapped_read_lengths_by_fragment_length_[first/second][fragment_length][read_length]
                                                     // = #reads
    uintQual phred_quality_offset_;
    uintQual minimum_quality_;
    uintQual maximum_quality_;
    uintReadLen minimum_read_length_on_reference_;
    uintReadLen maximum_read_length_on_reference_;
    double corrected_coverage_;

    // Collected variables for plotting
    PerSegment<Vect<uintFragCount>>
        gc_read_content_reference_; // gc_read_content_reference_[first/second][gcContentReference(%)] = #reads
    PerSegment<Vect<uintFragCount>>
        gc_read_content_mapped_; // gc_read_content_mapped_[first/second][gcContent(%)] = #reads

    PerSegmentPerStrand<PerBase<Vect<uintNucCount>>>
        sequence_content_reference_; // sequence_content_reference_[first/second][forward/reverse][A/C/G/T][readPosition]
                                     // = #(reads with given reference content at given reference position)

    // Collected variables for output
    uintFragCount total_number_reads_;

    // Private functions
    void Shrink();
    void PrepareGeneral();

    // boost serialization
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned int UNUSED(version)) {
        ar & adapters_;
        ar & coverage_;
        ar & errors_;
        ar & duplicates_;
        ar & fragment_distribution_;
        ar & qualities_;
        ar & tiles_;

        ar & creation_time_;
        ar & read_lengths_;
        ar & read_lengths_by_fragment_length_;
        ar & non_mapped_read_lengths_by_fragment_length_;
        ar & phred_quality_offset_;
        ar & minimum_quality_;
        ar & maximum_quality_;
        ar & minimum_read_length_on_reference_;
        ar & maximum_read_length_on_reference_;
        ar & corrected_coverage_;

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
    }

    // Google test
    friend class BamIngestionEngine;
    friend class DataStatsTest;
    friend class ProbabilityEstimatesTest;
    friend class SimulatorTest;

  public:
    DataStats(Reference* ref, uintSeqLen maximum_insert_length = 2000, uintQual minimum_mapping_quality = 10);

    // Getter functions
    const AdapterStats& Adapters() const { return adapters_; }
    const CoverageStats& Coverage() const { return coverage_; }
    const ErrorStats& Errors() const { return errors_; }
    const FragmentDuplicationStats& Duplicates() const { return duplicates_; }
    const FragmentDistributionStats& FragmentDistribution() const { return fragment_distribution_; }
    FragmentDistributionStats& FragmentDistribution() { return fragment_distribution_; }
    const QualityStats& Qualities() const { return qualities_; }
    const TileStats& Tiles() const { return tiles_; }

    uintQual PhredQualityOffset() const { return phred_quality_offset_; }
    uintFragCount TotalNumberReads() const { return total_number_reads_; }

    uint64_t CreationTime() const { return creation_time_; }

    const Vect<uintFragCount>& ReadLengths(uintTempSeq template_segment) const {
        return read_lengths_.at(template_segment);
    }
    uintReadLen MaxReadLenOnReference() const { return maximum_read_length_on_reference_; }
    const Vect<Vect<uintFragCount>>& ReadLengthsByFragmentLength(uintTempSeq template_segment) const {
        return read_lengths_by_fragment_length_.at(template_segment);
    }
    const Vect<Vect<uintFragCount>>& NonMappedReadLengthsByFragmentLength(uintTempSeq template_segment) const {
        return non_mapped_read_lengths_by_fragment_length_.at(template_segment);
    }

    double CorrectedCoverage() const { return corrected_coverage_; }
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
    const Vect<uintFragCount>& GCReadContentReference(uintTempSeq template_segment) const {
        return gc_read_content_reference_.at(template_segment);
    }
    const Vect<uintFragCount>& GCReadContentMapped(uintTempSeq template_segment) const {
        return gc_read_content_mapped_.at(template_segment);
    }

    const Vect<uintFragCount>& NContent(uintTempSeq template_segment) const {
        return read_sequence_stats_.NContent(template_segment);
    }
    const Vect<uintNucCount>& SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const {
        return read_sequence_stats_.SequenceContent(template_segment, nucleotide);
    }
    const Vect<uintNucCount>& SequenceContentReference(uintTempSeq template_segment, uintTempSeq strand,
                                                       uintBaseCall nucleotide) const {
        return sequence_content_reference_.at(template_segment).at(strand).at(nucleotide);
    }

    const Vect<uintNucCount>& HomopolymerDistribution(uintBaseCall nucleotide) const {
        return read_sequence_stats_.HomopolymerDistribution(nucleotide);
    }

    bool HasReference() const { return reference_; }

    // Setter functions
    void IgnoreTiles() { tiles_.IgnoreTiles(); }
    void ClearReference() { reference_ = nullptr; }
    void SetReference(Reference* ref) { reference_ = ref; }
    void SetUniformBias() { fragment_distribution_.SetUniformBias(); }

    // Public functions
    [[nodiscard]] bool IsValidRecord(const seqan::BamAlignmentRecord& record);
    inline bool QualitySufficient(const seqan::BamAlignmentRecord& record) const {
        return minimum_mapping_quality_ <= record.mapQ;
    }
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record);
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record, uintReadLen& max_indel);
    void GetReadPosOnReference(uintSeqLen& start_pos, uintSeqLen& end_pos,
                               const seqan::BamAlignmentRecord& record) const;
    inline bool InProperDirection(const seqan::BamAlignmentRecord& record_first, uintSeqLen end_pos_second,
                                  uintSeqLen start_pos_first, uintSeqLen start_pos_second) const {
        // Proper forward reverse direction or read sized pair (proper direction is only checked if they are on same
        // scaffold, so no check needed for that)
        return (!hasFlagRC(record_first) && hasFlagNextRC(record_first) &&
                maximum_insert_length_ >= end_pos_second - start_pos_first) ||
               start_pos_first == start_pos_second;
    }

    [[nodiscard]] bool ReadBam(const char* bam_file, const char* adapter_file, const char* adapter_matrix,
                               const std::string& variant_file, uintSeqLen max_ref_seq_bin_size,
                               uintNumThreads num_threads,
                               bool calculate_bias = true); // Fill the class with the information from a bam file

    [[nodiscard]] bool Load(const char* archive_file);
    [[nodiscard]] bool Save(const char* archive_file, bool text_format = false) const;

    void PrepareProcessing();
    void PreparePlotting();
    void PrepareTesting();

    void CalculateMaxLenDeletion() { errors_.PrepareSimulation(); }
};

} // namespace reseq

#endif // DATASTATS_H
