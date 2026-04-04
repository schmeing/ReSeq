#ifndef COVERAGESTATS_H
#define COVERAGESTATS_H

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <utility>

#include "logging.hpp"

#include <seqan/bam_io.h>

#include "container_types.hpp"
#include "ErrorStats.h"
#include "QualityStats.h"
#include "Reference.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq {
class CoverageStats {
  public:
    struct FullRecord {
        seqan::BamAlignmentRecord record_;
        uintSeqLen from_ref_pos_;
        uintSeqLen to_ref_pos_;
        uintTileId tile_id_;
        uintQual sequence_quality_;
        uintPercent reference_gc_;
        uintSeqLen fragment_length_;

        FullRecord()
            : from_ref_pos_(0), to_ref_pos_(0),
              sequence_quality_(255) { // Mark them as completely unprocessed (pair not found yet)
        }
    };

    struct CoveragePosition {
        std::array<std::atomic<uintCovCount>, kNumBasesN> coverage_forward_;
        std::array<std::atomic<uintCovCount>, kNumBasesN> coverage_reverse_;

        std::array<seqan::Dna5, kStrands> dom_error_;  // dom_error_[strand] = domError
        std::array<uintPercent, kStrands> error_rate_; // error_rate[strand] = errorRate

        bool valid_;
        std::array<bool, kStrands> coverage_sufficient_;

        CoveragePosition() : dom_error_({4, 4}), error_rate_({0, 0}), valid_(true), coverage_sufficient_({true, true}) {
            for (int i = kNumBasesN; i--;) {
                coverage_forward_.at(i) = 0;
                coverage_reverse_.at(i) = 0;
            }
        }

        CoveragePosition(const CoveragePosition&) = delete;
        CoveragePosition& operator=(const CoveragePosition&) = delete;

        CoveragePosition(CoveragePosition&& other) noexcept
            : dom_error_(other.dom_error_), error_rate_(other.error_rate_), valid_(other.valid_),
              coverage_sufficient_(other.coverage_sufficient_) {
            for (int i = kNumBasesN; i--;) {
                coverage_forward_.at(i).store(other.coverage_forward_.at(i).load(std::memory_order_relaxed),
                                              std::memory_order_relaxed);
                coverage_reverse_.at(i).store(other.coverage_reverse_.at(i).load(std::memory_order_relaxed),
                                              std::memory_order_relaxed);
            }
        }
    };

    struct ProcessedCoveragePosition {
        std::array<seqan::Dna5, kStrands> dom_error_;  // dom_error_[strand] = domError
        std::array<uintPercent, kStrands> error_rate_; // error_rate[strand] = errorRate

        bool valid_;
        std::array<bool, kStrands> coverage_sufficient_;

        ProcessedCoveragePosition()
            : dom_error_({4, 4}), error_rate_({0, 0}), valid_(false), coverage_sufficient_({false, false}) {}

        ProcessedCoveragePosition(const CoveragePosition& rhs)
            : dom_error_(rhs.dom_error_), error_rate_(rhs.error_rate_), valid_(rhs.valid_),
              coverage_sufficient_(rhs.coverage_sufficient_) {}
    };

    struct CoverageBlock {
        uintRefSeqId sequence_id_;
        uintSeqLen start_pos_;
        std::vector<CoveragePosition> coverage_;
        std::vector<ProcessedCoveragePosition>
            previous_coverage_; // End of coverage information from previous block (max read length dependent)
        std::atomic<uintFragCount> unprocessed_fragments_;
        std::vector<std::unique_ptr<FullRecord>> reads_;
        intVariantId first_variant_id_;
        std::atomic_flag scheduled_for_processing_;
        std::atomic<bool> processed_;

        size_t block_idx_;      // Index of this block in blocks_ deque
        size_t prev_block_idx_; // Index of previous block (SIZE_MAX = null)
        size_t next_block_idx_; // Index of next block (SIZE_MAX = null)

        CoverageBlock(uintRefSeqId seq_id, uintSeqLen start_pos)
            : sequence_id_(seq_id), start_pos_(start_pos), unprocessed_fragments_(0), first_variant_id_(0),
              block_idx_(SIZE_MAX), prev_block_idx_(SIZE_MAX), next_block_idx_(SIZE_MAX) {
            scheduled_for_processing_.clear();
            processed_ = false;
        }
    };

    class ThreadData {
      private:
        std::vector<std::array<uintCovCount, kStrands>>
            block_coverage_; // block_coverage_[BlockPos][forward/reverse] = coverage
        std::vector<std::pair<double, uintCovCount>> error_rates_sorted_;
        std::vector<std::array<double, kStrands>>
            non_sytematic_probability_; // block_coverage_[BlockPos][forward/reverse] = prob
        std::vector<double> non_sytematic_probability_sorted_;

        friend class CoverageStats;

      public:
        ThreadData() {
            block_coverage_.reserve(CoverageStats::kBlockSize);
            error_rates_sorted_.reserve(static_cast<size_t>(2) * CoverageStats::kBlockSize);
            non_sytematic_probability_.reserve(CoverageStats::kBlockSize);
            non_sytematic_probability_sorted_.reserve(static_cast<size_t>(2) * CoverageStats::kBlockSize);
        }
    };

  private:
    // Definitions
    static constexpr uintSeqLen kBlockSize = 10000;
    static constexpr uintCovCount kMaxCoverage = 10000; // Maximum counted coverage (Only relevant for plotting)
    const double kSystematicErrorFDR = 0.05;
    const uint32_t kPValueHistBins = 100;

    // Data dependent parameter
    uintCovCount coverage_threshold_; // Minimium coverage to accept value from this position into error_rates_ and
                                      // dominant_errors_
    uintSeqLen reset_distance_;       // Maximum distance between two systematic errors to continue error region
    uintSeqLen gc_range_;
    uintReadLen maximum_read_length_on_reference_;
    uintCovCount genome_coverage_average_;

    // Mutex
    std::mutex clean_up_mutex_;
    std::mutex reuse_mutex_;
    std::mutex blocks_growth_mutex_; // Protects blocks_.emplace_back() — deque growth is not thread-safe
    std::mutex variant_loading_mutex_;

    // Temporary variables
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_dominant_errors_by_distance_;
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_dominant_errors_by_gc_;
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_gc_by_distance_de_;
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_dominant_errors_by_start_rates_;
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_start_rates_by_distance_de_;
    PerBase<PerBaseN<PerBase<std::vector<AtomicVec<uintNucCount>>>>> tmp_start_rates_by_gc_de_;

    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_error_rates_by_distance_;
    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_error_rates_by_gc_;
    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_gc_by_distance_er_;
    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_error_rates_by_start_rates_;
    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_start_rates_by_distance_er_;
    PerBase<PerBaseN<std::vector<AtomicVec<uintNucCount>>>> tmp_start_rates_by_gc_er_;

    AtomicVec<uintSurBlockId> tmp_block_error_rate_;
    AtomicVec<uintSurBlockId> tmp_block_percent_systematic_;
    AtomicVec<uintNucCount> tmp_systematic_error_p_values_;

    AtomicVec<uintNucCount> tmp_coverage_;
    PerStrand<AtomicVec<uintNucCount>> tmp_coverage_stranded_;
    PerStrand<AtomicVec<uintNucCount>> tmp_coverage_stranded_percent_;
    PerStrand<AtomicVec<uintNucCount>> tmp_coverage_stranded_percent_min_cov_10_;
    PerStrand<AtomicVec<uintNucCount>> tmp_coverage_stranded_percent_min_cov_20_;
    AtomicVec<uintNucCount> tmp_error_coverage_;
    AtomicVec<uintNucCount> tmp_error_coverage_percent_;
    AtomicVec<uintNucCount> tmp_error_coverage_percent_min_cov_10_;
    AtomicVec<uintNucCount> tmp_error_coverage_percent_min_cov_20_;
    std::vector<AtomicVec<uintNucCount>> tmp_error_coverage_percent_stranded_;
    std::vector<AtomicVec<uintNucCount>> tmp_error_coverage_percent_stranded_min_strand_cov_10_;
    std::vector<AtomicVec<uintNucCount>> tmp_error_coverage_percent_stranded_min_strand_cov_20_;

    // Collected variables for simulation
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        dominant_errors_by_distance_; // dominant_errors_by_distance_[refBase][previousRefBase][domRefBaseLast5][(distanceToStartOfErrorRegion+9)/10][dominantError]
                                      // = #refBases
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        dominant_errors_by_gc_; // dominant_errors_by_gc_[refBase][previousRefBase][domRefBaseLast5][GClastHalfAverageReadLength][dominantError]
                                // = #refBases
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        gc_by_distance_de_; // gc_by_distance_de_[refBase][previousRefBase][domRefBaseLast5][(distanceToStartOfErrorRegion+9)/10][GClastHalfAverageReadLength]
                            // = #refBases
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        dominant_errors_by_start_rates_; // dominant_errors_by_start_rates_[refBase][previousRefBase][domRefBaseLast5][errorRateStart][dominantError]
                                         // = #refBases
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        start_rates_by_distance_de_; // start_rates_by_distance_de_[refBase][previousRefBase][domRefBaseLast5][(distanceToStartOfErrorRegion+9)/10][errorRateStart]
                                     // = #refBases
    PerBase<PerBaseN<PerBase<Vect<Vect<uintNucCount>>>>>
        start_rates_by_gc_de_; // start_rates_by_gc_de_[refBase][previousRefBase][domRefBaseLast5][(distanceToStartOfErrorRegion+9)/10][errorRateStart]
                               // = #refBases

    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        error_rates_by_distance_; // error_rates_by_distance_[refBase][dominantError][(distanceToStartOfErrorRegion+9)/10][errorRate]
                                  // = #refBases
    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        error_rates_by_gc_; // error_rates_by_gc_[refBase][dominantError][GClastHalfAverageReadLength][errorRate]
                            // = #refBases
    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        gc_by_distance_er_; // gc_by_distance_er_[refBase][dominantError][(distanceToStartOfErrorRegion+9)/10][GClastHalfAverageReadLength]
                            // = #refBases
    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        error_rates_by_start_rates_; // error_rates_by_start_rates_[refBase][dominantError][errorRateStart][errorRate]
                                     // = #refBases
    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        start_rates_by_distance_er_; // start_rates_by_distance_[refBase][dominantError][(distanceToStartOfErrorRegion+9)/10][errorRateStart]
                                     // = #refBases
    PerBase<PerBaseN<Vect<Vect<uintNucCount>>>>
        start_rates_by_gc_er_; // start_rates_by_gc_[refBase][dominantError][GClastHalfAverageReadLength][errorRateStart]
                               // = #refBases

    // Collected variables for plotting
    Vect<uintSurBlockId> block_error_rate_;
    Vect<uintSurBlockId> block_percent_systematic_;
    Vect<uintNucCount> systematic_error_p_values_;

    Vect<uintNucCount> coverage_; // coverage_[ coverageDepth ] = #bases
    PerStrand<Vect<uintNucCount>>
        coverage_stranded_; // coverage_stranded_[forward/reverse][ coverageDepthOnStrand ] = #bases
    PerStrand<Vect<uintNucCount>> coverage_stranded_percent_; // coverage_stranded_percent_[forward/reverse][
                                                              // coverageDepthOnStrand/coverageDepth*100 ] = #bases
    PerStrand<Vect<uintNucCount>>
        coverage_stranded_percent_min_cov_10_; // coverage_stranded_percent_min_cov_10_[forward/reverse][
                                               // coverageDepthOnStrand/coverageDepth*100 ] = #bases
    PerStrand<Vect<uintNucCount>>
        coverage_stranded_percent_min_cov_20_; // coverage_stranded_percent_min_cov_20_[forward/reverse][
                                               // coverageDepthOnStrand/coverageDepth*100 ] = #bases
    Vect<uintNucCount> error_coverage_;        // error_coverage_[ #errorsAtReferencePosition ] = #bases
    Vect<uintNucCount>
        error_coverage_percent_; // error_coverage_percent_[ #errorsAtReferencePosition/coverage*100 ] = #bases
    Vect<uintNucCount> error_coverage_percent_min_cov_10_; // error_coverage_percent_min_cov_10_[
                                                           // #errorsAtReferencePosition/coverage*100 ] = #bases
    Vect<uintNucCount> error_coverage_percent_min_cov_20_; // error_coverage_percent_min_cov_20_[
                                                           // #errorsAtReferencePosition/coverage*100 ] = #bases
    Vect<Vect<uintNucCount>>
        error_coverage_percent_stranded_; // error_coverage_percent_stranded_[
                                          // #forwardErrorsAtReferencePosition/coverage*100 ][
                                          // #reverseErrorsAtReferencePosition/coverage*100 ] = #bases
    Vect<Vect<uintNucCount>>
        error_coverage_percent_stranded_min_strand_cov_10_; // error_coverage_percent_stranded_min_strand_cov_10_[
                                                            // #forwardErrorsAtReferencePosition/coverage*100 ][
                                                            // #reverseErrorsAtReferencePosition/coverage*100 ] = #bases
    Vect<Vect<uintNucCount>>
        error_coverage_percent_stranded_min_strand_cov_20_; // error_coverage_percent_stranded_min_strand_cov_20_[
                                                            // #forwardErrorsAtReferencePosition/coverage*100 ][
                                                            // #reverseErrorsAtReferencePosition/coverage*100 ] = #bases

    // Calculated variables for plotting
    Vect<Vect<uintNucCount>>
        error_rates_by_distance_sum_; // error_rates_by_distance_sum_[(distanceToStartOfErrorRegion+9)/10][errorRate]
                                      // = #refBases
    Vect<Vect<uintNucCount>>
        error_rates_by_gc_sum_; // error_rates_by_gc_sum_[GClastHalfAverageReadLength][errorRate] = #refBases

    // Temporary variables for read in
    std::deque<std::unique_ptr<CoverageBlock>> blocks_; // Owns all CoverageBlock instances; deque for pointer stability
    std::vector<size_t> free_indices_;                  // Recycled block slots
    std::atomic<size_t> first_live_idx_{SIZE_MAX};      // Front of live range
    std::atomic<size_t> last_live_idx_{SIZE_MAX};       // End of live range (publication signal)

    uintRefLenCalc zero_coverage_region_;
    uintRefLenCalc excluded_bases_;
    uintRefLenCalc num_exclusion_regions_;

    // Helper functions
    inline bool CoveragePosValid(intSeqShift coverage_pos, CoverageStats::CoverageBlock* coverage_block) const {
        if (0 > coverage_pos) {
            return coverage_block->previous_coverage_.at(-1 * coverage_pos - 1).valid_;
        } else {
            return coverage_block->coverage_.at(coverage_pos).valid_;
        }
    }
    void EvalRead(FullRecord* record, CoverageStats::CoverageBlock* coverage_block, const Reference& reference,
                  QualityStats& qualities, ErrorStats& errors, uintQual phred_quality_offset);

    inline void CountEmptyEndOfSequence(const Reference& reference) {
        auto& lb = *blocks_[last_live_idx_.load(std::memory_order_relaxed)];
        zero_coverage_region_ += reference.SequenceLength(lb.sequence_id_) - lb.start_pos_ - lb.coverage_.size();
    }
    inline void ApplyZeroCoverageRegion();

    double GetPositionProbabilities(CoverageBlock* block, const Reference& reference, ThreadData& thread);
    void UpdateCoverageAtSinglePosition(CoveragePosition& nuc_coverage, std::array<uintCovCount, 2>& coverage,
                                        seqan::Dna5 ref_base);
    inline bool NextBlockWithInSysErrorResetDistance(CoverageBlock* block) {
        if (block->next_block_idx_ == SIZE_MAX) {
            return false;
        }
        auto& next = *blocks_[block->next_block_idx_];
        return next.sequence_id_ == block->sequence_id_ &&
               reset_distance_ - 1 > next.start_pos_ - (block->start_pos_ + block->coverage_.size());
    }
    inline uintSeqLen BasesWithInSysErrorResetDistance(CoverageBlock* block) {
        auto& next = *blocks_[block->next_block_idx_];
        return reset_distance_ - 1 -
               (next.start_pos_ -
                (block->start_pos_ + block->coverage_.size())); // reset_distance_-1 because we don't need the last base
                                                                // that we ignore but the first that we do not reset
    }
    inline void AddSysError(seqan::Dna5 ref_base, seqan::Dna5 prev_ref_base, seqan::Dna5 dom_ref_base,
                            seqan::Dna5 dom_error, uintPercent error_rate, uintSeqLen distance, uintPercent gc_percent,
                            uintPercent start_rate) {
        ++tmp_dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_ref_base).at(distance).at(dom_error);
        ++tmp_dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_ref_base).at(gc_percent).at(dom_error);
        ++tmp_gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base).at(distance).at(gc_percent);
        ++tmp_dominant_errors_by_start_rates_.at(ref_base)
              .at(prev_ref_base)
              .at(dom_ref_base)
              .at(start_rate)
              .at(dom_error);
        ++tmp_start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base).at(distance).at(start_rate);
        ++tmp_start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base).at(gc_percent).at(start_rate);

        ++tmp_error_rates_by_distance_.at(ref_base).at(dom_error).at(distance).at(error_rate);
        ++tmp_error_rates_by_gc_.at(ref_base).at(dom_error).at(gc_percent).at(error_rate);
        ++tmp_gc_by_distance_er_.at(ref_base).at(dom_error).at(distance).at(gc_percent);
        ++tmp_error_rates_by_start_rates_.at(ref_base).at(dom_error).at(start_rate).at(error_rate);
        ++tmp_start_rates_by_distance_er_.at(ref_base).at(dom_error).at(distance).at(start_rate);
        ++tmp_start_rates_by_gc_er_.at(ref_base).at(dom_error).at(gc_percent).at(start_rate);
    }

    CoverageBlock* CreateBlock(uintRefSeqId seq_id, uintSeqLen start_pos);
    void UpdateFirstVariant(CoverageBlock& block, const Reference& reference);
    bool IsVariantPosition(intVariantId& cur_var, uintRefSeqId ref_seq_id, uintSeqLen pos,
                           const Reference& reference) const;
    void InitBlock(CoverageBlock& block, const Reference& reference);
    bool NonSystematicError(uintCovCount coverage, uintCovCount errors, double prob, double max_prob) {
        return 1 >= errors || 0 == utilities::Percent(errors, coverage) || prob > max_prob;
    }
    void ProcessBlock(CoverageBlock* block, const Reference& reference, ThreadData& thread);
    void CountBlock(CoverageBlock* block, const Reference& reference);
    CoverageBlock* RemoveBlock(CoverageBlock* block);

    // Boost archive functions
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive& ar, const unsigned int UNUSED(version)) {
        ar & coverage_threshold_;
        ar & reset_distance_;

        ar & dominant_errors_by_distance_;
        ar & dominant_errors_by_gc_;
        ar & gc_by_distance_de_;
        ar & dominant_errors_by_start_rates_;
        ar & start_rates_by_distance_de_;
        ar & start_rates_by_gc_de_;

        ar & error_rates_by_distance_;
        ar & error_rates_by_gc_;
        ar & gc_by_distance_er_;
        ar & error_rates_by_start_rates_;
        ar & start_rates_by_distance_er_;
        ar & start_rates_by_gc_er_;

        ar & block_error_rate_;
        ar & block_percent_systematic_;
        ar & systematic_error_p_values_;

        ar & coverage_;
        ar & coverage_stranded_;
        ar & coverage_stranded_percent_;
        ar & coverage_stranded_percent_min_cov_10_;
        ar & coverage_stranded_percent_min_cov_20_;
        ar & error_coverage_;
        ar & error_coverage_percent_;
        ar & error_coverage_percent_min_cov_10_;
        ar & error_coverage_percent_min_cov_20_;
        ar & error_coverage_percent_stranded_;
        ar & error_coverage_percent_stranded_min_strand_cov_10_;
        ar & error_coverage_percent_stranded_min_strand_cov_20_;
    }

    // Google test
    friend class CoverageStatsTest;
    friend class ProbabilityEstimatesTest;

  public:
    CoverageStats()
        : coverage_threshold_(100), zero_coverage_region_(0), excluded_bases_(0), num_exclusion_regions_(0) {}

    // Getter functions
    inline const Vect<Vect<uintNucCount>>& DominantErrorsByDistance(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                                    seqan::Dna5 dom_last5) const {
        return dominant_errors_by_distance_.at(ref_base).at(last_ref_base).at(dom_last5);
    }
    inline const Vect<Vect<uintNucCount>>& DominantErrorsByGC(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                              seqan::Dna5 dom_last5) const {
        return dominant_errors_by_gc_.at(ref_base).at(last_ref_base).at(dom_last5);
    }
    inline const Vect<Vect<uintNucCount>>& GCByDistance(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                        seqan::Dna5 dom_last5) const {
        return gc_by_distance_de_.at(ref_base).at(last_ref_base).at(dom_last5);
    }
    inline const Vect<Vect<uintNucCount>>& DominantErrorsByStartRates(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                                      seqan::Dna5 dom_last5) const {
        return dominant_errors_by_start_rates_.at(ref_base).at(last_ref_base).at(dom_last5);
    }
    inline const Vect<Vect<uintNucCount>>& StartRatesByDistance(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                                seqan::Dna5 dom_last5) const {
        return start_rates_by_distance_de_.at(ref_base).at(last_ref_base).at(dom_last5);
    }
    inline const Vect<Vect<uintNucCount>>& StartRatesByGC(seqan::Dna ref_base, seqan::Dna5 last_ref_base,
                                                          seqan::Dna5 dom_last5) const {
        return start_rates_by_gc_de_.at(ref_base).at(last_ref_base).at(dom_last5);
    }

    inline const Vect<Vect<uintNucCount>>& ErrorRatesByDistance(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return error_rates_by_distance_.at(ref_base).at(dom_error);
    }
    inline const Vect<Vect<uintNucCount>>& ErrorRatesByGC(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return error_rates_by_gc_.at(ref_base).at(dom_error);
    }
    inline const Vect<Vect<uintNucCount>>& GCByDistance(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return gc_by_distance_er_.at(ref_base).at(dom_error);
    }
    inline const Vect<Vect<uintNucCount>>& ErrorRatesByStartRates(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return error_rates_by_start_rates_.at(ref_base).at(dom_error);
    }
    inline const Vect<Vect<uintNucCount>>& StartRatesByDistance(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return start_rates_by_distance_er_.at(ref_base).at(dom_error);
    }
    inline const Vect<Vect<uintNucCount>>& StartRatesByGC(seqan::Dna ref_base, seqan::Dna5 dom_error) const {
        return start_rates_by_gc_er_.at(ref_base).at(dom_error);
    }

    inline const Vect<Vect<uintNucCount>>& ErrorRatesByDistanceSum() const { return error_rates_by_distance_sum_; }
    inline const Vect<Vect<uintNucCount>>& ErrorRatesByGCSum() const { return error_rates_by_gc_sum_; }

    inline const Vect<uintSurBlockId>& BlockErrorRate() const { return block_error_rate_; }
    inline const Vect<uintSurBlockId>& BlockPercentSystematic() const { return block_percent_systematic_; }
    inline const Vect<uintNucCount>& SystematicErrorPValues() const { return systematic_error_p_values_; }

    inline const Vect<uintNucCount>& Coverage() const { return coverage_; }
    inline const Vect<uintNucCount>& CoverageStranded(bool reverse_strand) const {
        return coverage_stranded_.at(reverse_strand);
    }
    inline const Vect<uintNucCount>& CoverageStrandedPercent(bool reverse_strand) const {
        return coverage_stranded_percent_.at(reverse_strand);
    }
    inline const Vect<uintNucCount>& CoverageStrandedPercentMinCov10(bool reverse_strand) const {
        return coverage_stranded_percent_min_cov_10_.at(reverse_strand);
    }
    inline const Vect<uintNucCount>& CoverageStrandedPercentMinCov20(bool reverse_strand) const {
        return coverage_stranded_percent_min_cov_20_.at(reverse_strand);
    }
    inline const Vect<uintNucCount>& ErrorCoverage() const { return error_coverage_; }
    inline const Vect<uintNucCount>& ErrorCoveragePercent() const { return error_coverage_percent_; }
    inline const Vect<uintNucCount>& ErrorCoveragePercentMinCov10() const { return error_coverage_percent_min_cov_10_; }
    inline const Vect<uintNucCount>& ErrorCoveragePercentMinCov20() const { return error_coverage_percent_min_cov_20_; }
    inline const Vect<Vect<uintNucCount>>& ErrorCoveragePercentStranded() const {
        return error_coverage_percent_stranded_;
    }
    inline const Vect<Vect<uintNucCount>>& ErrorCoveragePercentStrandedMinCov10() const {
        return error_coverage_percent_stranded_min_strand_cov_10_;
    }
    inline const Vect<Vect<uintNucCount>>& ErrorCoveragePercentStrandedMinCov20() const {
        return error_coverage_percent_stranded_min_strand_cov_20_;
    }

    // Setter functions
    inline void AddForward(uintSeqLen pos, CoverageBlock* block, seqan::Dna5 called_base) {
        ++(block->coverage_.at(pos).coverage_forward_.at(called_base));
    }
    inline void AddReverse(uintSeqLen pos, CoverageBlock* block, seqan::Dna5 called_base) {
        ++(block->coverage_.at(pos).coverage_reverse_.at(called_base));
    }

    // Main functions
    void Prepare(uintCovCount average_coverage, uintReadLen average_read_length,
                 uintReadLen maximum_read_length_on_reference);

    CoverageBlock* FindBlock(uintRefSeqId ref_seq_id, uintSeqLen ref_pos);
    bool EnsureSpace(uintRefSeqId ref_seq_id, uintSeqLen start_pos, uintSeqLen end_pos, FullRecord* record,
                     Reference& reference);

    void AddFragment(uintRefSeqId ref_seq_id, uintSeqLen ref_pos, CoverageBlock*& block);
    void RemoveFragment(uintRefSeqId ref_seq_id, uintSeqLen ref_pos, CoverageBlock*& block,
                        uintFragCount& processed_fragments);

    inline uintSeqLen GetStartPos(uintSeqLen pos, CoverageBlock*& start_block) const {
        while (pos >= start_block->start_pos_ + kBlockSize) {
            start_block = blocks_[start_block->next_block_idx_].get();
        }
        return pos - start_block->start_pos_;
    }
    inline void IncrementPos(uintSeqLen& pos, CoverageBlock*& block) const {
        if (++pos >= kBlockSize) {
            pos = 0;
            block = blocks_[block->next_block_idx_].get();
        }
    }
    inline void AddPos(uintSeqLen& pos, CoverageBlock*& block, uintSeqLen value) const {
        pos += value;
        while (pos >= kBlockSize) {
            pos -= kBlockSize;
            block = blocks_[block->next_block_idx_].get();
        }
    }
    inline void DecrementPos(uintSeqLen& pos, CoverageBlock*& block) const {
        if (0 == pos) {
            pos = kBlockSize;
            block = blocks_[block->prev_block_idx_].get();
        }
        --pos;
    }
    inline void SubtractPos(uintSeqLen& pos, CoverageBlock*& block, uintSeqLen value) const {
        if (value <= pos) {
            pos -= value;
        } else {
            value -= pos;
            block = blocks_[block->prev_block_idx_].get();
            while (value > kBlockSize) {
                value -= kBlockSize;
                block = blocks_[block->prev_block_idx_].get();
            }
            pos = kBlockSize - value;
        }
    }
    void UpdateDistances(uintSeqLen& distance_to_start_of_error_region, uintPercent& start_rate,
                         uintPercent error_rate) const;
    uintRefSeqId CleanUp(uintSeqLen& still_needed_position, Reference& reference, QualityStats& qualities,
                         ErrorStats& errors, uintQual phred_quality_offset, CoverageBlock* cov_block,
                         uintFragCount& processed_fragments, ThreadData& thread);
    bool PreLoadVariants(Reference& reference);
    bool Finalize(const Reference& reference, QualityStats& qualities, ErrorStats& errors,
                  uintQual phred_quality_offset, std::mutex& print_mutex, ThreadData& thread);
    inline void SetAllZero(const Reference& reference) {
        zero_coverage_region_ = reference.TotalSize();
        for (auto seq_id = reference.NumberSequences(); seq_id--;) {
            excluded_bases_ += reference.SumExcludedBases(seq_id);
            num_exclusion_regions_ += reference.NumExcludedRegions(seq_id);
        }
    }

    void Shrink();
    void PreparePlotting();
    uintCovCount CoveragePeakPosition() const {
        uintCovCount peak(coverage_.from());
        // As long as coverage is strictly monotonously falling it still belongs to the uncovered parts of the genome
        while (peak + 1 < coverage_.to() && coverage_.at(peak + 1) < coverage_.at(peak)) {
            ++peak;
        }

        for (auto cov = peak; cov < coverage_.to(); ++cov) {
            if (coverage_.at(peak) < coverage_.at(cov)) {
                peak = cov;
            }
        }

        return peak;
    }
    double MeanCoverage() const {
        uintNucCount sum(0), count(0);
        for (auto cov = coverage_.from(); cov < coverage_.to(); ++cov) {
            sum += cov * coverage_.at(cov);
            count += coverage_.at(cov);
        }

        if (count) {
            return static_cast<double>(sum) / count;
        } else {
            return 0.0;
        }
    }
};
} // namespace reseq

#endif // COVERAGESTATS_H
