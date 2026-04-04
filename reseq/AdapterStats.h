#ifndef ADAPTERSTATS_H
#define ADAPTERSTATS_H

#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

#include <boost/serialization/vector.hpp>

#include <seqan/bam_io.h>

#include "container_types.hpp"
#include "Reference.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq {
class AdapterStats {
  private:
    // Definitions
    const uintReadLen kMinimumAdapterLengthUnmapped =
        15; // Minimum adapter length required to insert not completely mapped read pair into statistics (15 should stop
            // bowtie2 from mapping correctly)
    const double kMinFractionOfMaximumForSimulation =
        0.1; // Minimum fraction of the highest adapter count that an adapter needs to have to be simulated (Removes
             // garbage from adapter reading stage)
    static constexpr uintReadLen kKmerLength = 10;

    static constexpr const char* kAdapterSearchInfoFile = nullptr; // "adapter.csv";

    // Temporary variables
    PerSegment<KmerCount<kKmerLength>> adapter_kmers_;          // adapter_kmers_[first/second] = KmerCounts
    PerSegment<std::vector<uintNucCount>> adapter_start_kmers_; // adapter_start_kmers_[first/second] = StartKmerCount

    std::vector<std::vector<std::vector<AtomicVec<uintFragCount>>>>
        tmp_counts_; // counts_[AdapterID][secondAdapterID][firstAdapterLength][secondAdapterLength] = #adapters
    PerSegment<std::vector<AtomicVec<uintFragCount>>>
        tmp_start_cut_; // start_cut_[templateSegment][AdapterID][basesCutFromTheBeginningOfTheAdapterAtPos0]
                        // = #adapters
    AtomicVec<uintFragCount> tmp_polya_tail_length_; // polya_tail_length_[lengthOfPolyATailAfterAdapter] = #adapters
    PerBaseN<std::atomic<uintNucCount>>
        tmp_overrun_bases_; // overrun_bases_[nucleotide] = #basesAfterPolyATailWithThisNucleotide

    // Adapter infos loaded/detected
    PerSegment<std::vector<std::string>> names_;     // names_[templateSegment][AdapterID] = adapterNameFromFasta
    PerSegment<std::vector<seqan::DnaString>> seqs_; // seqs_[templateSegment][AdapterID] = adapterSequenceFromFasta
    std::vector<std::vector<bool>>
        combinations_; // combinations_[AdapterID][secondAdapterID] = notValid(false)/Valid(true)

    // Collected statistics
    std::vector<std::vector<Vect<Vect<uintFragCount>>>>
        counts_; // counts_[AdapterID][secondAdapterID][firstAdapterLength][secondAdapterLength] = #adapters
    PerSegment<std::vector<Vect<uintFragCount>>>
        start_cut_; // start_cut_[templateSegment][AdapterID][basesCutFromTheBeginningOfTheAdapterAtPos0] = #adapters
    Vect<uintFragCount> polya_tail_length_; // polya_tail_length_[lengthOfPolyATailAfterAdapter] = #adapters
    PerBaseN<uintNucCount> overrun_bases_;  // overrun_bases_[nucleotide] = #basesAfterPolyATailWithThisNucleotide

    // Calculated variables
    PerSegment<std::vector<uintFragCount>> count_sum_; // count_sum_[templateSegment][AdapterID] = #adapters
    PerSegment<std::vector<uintFragCount>>
        significant_count_; // significant_count_[templateSegment][AdapterID] = #adapters (All adapters not appearing
                            // often are set to zero)

    // Helper functions
    static inline uintSeqLen GetStartPosOnReference(const seqan::BamAlignmentRecord& record);
    static uintReadLen CountErrors(uintReadLen& last_error_pos, const seqan::BamAlignmentRecord& record,
                                   const Reference& reference);
    bool VerifyInsert(const seqan::CharString& read1, const seqan::CharString& read2, intSeqShift pos1,
                      intSeqShift pos2);
    static bool AdaptersAmbigous(const seqan::DnaString& adaper1, const seqan::DnaString& adaper2,
                                 uintReadLen max_length);

    // Boost archive functions
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned int UNUSED(version)) {
        ar & names_;
        ar & combinations_;

        ar & counts_;
        ar & start_cut_;
        ar & polya_tail_length_;
        ar & overrun_bases_;

        // Work around so that seqan>>DnaString can be stored as no serialization function exists for that class
        std::array<std::vector<std::string>, kTemplateSegments> seqs_archive;
        if (seqs_.at(0).size()) {
            // Serialization writing
            for (uintTempSeq template_segment = kTemplateSegments; template_segment--;) {
                seqs_archive.at(template_segment).reserve(seqs_.at(template_segment).size());
                for (auto& adapter : seqs_.at(template_segment)) {
                    seqs_archive.at(template_segment)
                        .emplace_back(seqan::toCString(static_cast<seqan::CharString>(adapter)));
                }
            }
        }

        ar & seqs_archive;

        if (!seqs_.at(0).size()) {
            // Serialization loading
            for (uintTempSeq template_segment = kTemplateSegments; template_segment--;) {
                seqs_.at(template_segment).clear();
                seqs_.at(template_segment).reserve(seqs_archive.at(template_segment).size());
                for (auto& adapter : seqs_archive.at(template_segment)) {
                    seqs_.at(template_segment).emplace_back(adapter.c_str());
                }
            }
        }
    }

    // Google test
    friend class AdapterStatsTest;

  public:
    //  Constructor
    AdapterStats();

    // Getter functions
    const std::string& Name(uintTempSeq template_segment, uintAdapterId id) const {
        return names_.at(template_segment).at(id);
    }
    const seqan::DnaString& Sequence(uintTempSeq template_segment, uintAdapterId id) const {
        return seqs_.at(template_segment).at(id);
    }
    const std::vector<uintFragCount>& Counts(uintTempSeq template_segment) const {
        return count_sum_.at(template_segment);
    }
    const std::vector<uintFragCount>& SignificantCounts(uintTempSeq template_segment) const {
        return significant_count_.at(template_segment);
    }
    const Vect<uintFragCount>& StartCut(uintTempSeq template_segment, uintAdapterId id) const {
        return start_cut_.at(template_segment).at(id);
    }
    const Vect<uintFragCount>& PolyATailLength() const { return polya_tail_length_; }
    const std::array<uintNucCount, kNumBasesN>& OverrunBases() const { return overrun_bases_; }

    // Main functions
    [[nodiscard]] bool LoadAdapters(const char* adapter_file, const char* adapter_matrix);
    void PrepareAdapterPrediction();
    void ExtractAdapterPart(const seqan::BamAlignmentRecord& record);
    [[nodiscard]] bool PredictAdapters();
    void PrepareAdapters(uintReadLen size_read_length, uintQual phred_quality_offset);
    [[nodiscard]] bool Detect(uintReadLen& adapter_position_first, uintReadLen& adapter_position_second,
                              const seqan::BamAlignmentRecord& record_first,
                              const seqan::BamAlignmentRecord& record_second, const Reference& reference,
                              bool properly_mapped = false);
    void Finalize();
    void SumCounts();
    void Shrink();
    void PrepareSimulation();
};
} // namespace reseq

#endif // ADAPTERSTATS_H
