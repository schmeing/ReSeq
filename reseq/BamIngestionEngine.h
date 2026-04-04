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
#include "FragmentDistributionStats.h"
#include "utilities.hpp"

namespace reseq {

// Forward declaration
class DataStats;

class BamIngestionEngine {
  public:
    // Operator Definitions (moved from DataStats)
    struct RecordHasher {
        size_t operator()(CoverageStats::FullRecord* const& record) const;
    };
    struct RecordEqual {
        bool operator()(CoverageStats::FullRecord* const& lhs, CoverageStats::FullRecord* const& rhs) const;
    };

    [[nodiscard]] bool Run(const char* bam_file, const char* adapter_file, const char* adapter_matrix,
                           const std::string& variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
                           bool calculate_bias, DataStats& target);

    // Static geometry utilities (unchanged signatures)
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record);
    static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord& record, uintReadLen& max_indel);

  private:
    struct ThreadData {
        std::vector<std::pair<CoverageStats::FullRecord*, CoverageStats::FullRecord*>> rec_store_;

        CoverageStats::ThreadData coverage_;
        FragmentDistributionStats::ThreadData fragment_distribution_;

        uintSeqLen last_exclusion_region_id_;
        uintRefSeqId last_exclusion_ref_seq_;

        ThreadData(uintSeqLen maximum_insert_length, uintSeqLen start_exclusion_length)
            : fragment_distribution_(maximum_insert_length, start_exclusion_length), last_exclusion_region_id_(0),
              last_exclusion_ref_seq_(0) {}
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
    inline bool PotentiallyValidGeneral(const seqan::BamAlignmentRecord& record, const DataStats& target) const;
    inline bool PotentiallyValidFirst(const seqan::BamAlignmentRecord& record_first, const DataStats& target) const;
    inline bool PotentiallyValidSecond(const seqan::BamAlignmentRecord& record_second, const DataStats& target) const;
    inline bool PotentiallyValid(const seqan::BamAlignmentRecord& record, const DataStats& target) const;

    bool IsSecondRead(CoverageStats::FullRecord* record, CoverageStats::FullRecord*& record_first,
                      CoverageStats::CoverageBlock*& block, DataStats& target);

    bool CheckForAdapters(const seqan::BamAlignmentRecord& record_first, const seqan::BamAlignmentRecord& record_second,
                          DataStats& target);
    bool EvalReferenceStatistics(CoverageStats::FullRecord* record, uintTempSeq template_segment,
                                 CoverageStats::CoverageBlock* coverage_block, DataStats& target);
    bool EvalRecord(std::pair<CoverageStats::FullRecord*, CoverageStats::FullRecord*> record, ThreadData& thread,
                    DataStats& target);

    bool SignsOfPairsWithNamesNotIdentical(const DataStats& target);
    void PrepareReadIn(uintQual size_mapping_quality, uintReadLen size_indel, uintSeqLen max_ref_seq_bin_size,
                       uintNumThreads num_threads, DataStats& target);
    void FinishReadIn(DataStats& target);
    void Shrink(DataStats& target);
    bool Calculate(uintNumThreads num_threads, DataStats& target);
    bool OrderOfBamFileCorrect(const seqan::BamAlignmentRecord& record,
                               std::pair<uintRefSeqId, uintSeqLen> last_record_pos);
    bool PreRun(seqan::BamFileIn& bam, const char* bam_file, seqan::BamHeader& header, uintQual& size_mapping_quality,
                uintReadLen& size_indel, DataStats& target);
    bool ReadRecords(seqan::BamFileIn& bam, bool& not_done, ThreadData& thread_data, DataStats& target);

    static void ReadThread(BamIngestionEngine& engine, seqan::BamFileIn& bam, size_t thread_idx, DataStats& target);

    // Google test — DataStatsTest checks first_read_records_.size() == 0
    friend class DataStatsTest;
};

} // namespace reseq

#endif // BAMINGESTIONENGINE_H
