#ifndef QUALITYSTATS_H
#define QUALITYSTATS_H

#include <array>
#include <stdint.h>
#include <vector>

#include "SeqQualityStats.hpp"
#include "utilities.h"
#include "Vect.hpp"

namespace reseq{
	class QualityStats{
	private:
		// Temporary variables
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>,5>, 4>, 2> tmp_base_quality_stats_per_tile_per_error_reference_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>,5>, 4>, 2> tmp_error_rate_for_position_per_tile_per_error_reference_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>,5>, 4>, 2> tmp_base_quality_for_error_rate_per_tile_per_error_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_base_quality_for_preceding_quality_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_preceding_quality_for_error_rate_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_preceding_quality_for_position_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_base_quality_for_sequence_quality_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_preceding_quality_for_sequence_quality_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_sequence_quality_for_error_rate_per_tile_reference_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 4>, 2> tmp_sequence_quality_for_position_per_tile_reference_;

		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_sequence_quality_mean_for_gc_per_tile_reference_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_sequence_quality_mean_for_mean_error_rate_per_tile_reference_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_sequence_quality_mean_for_fragment_length_per_tile_reference_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_mean_error_rate_for_gc_per_tile_reference_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_mean_error_rate_for_fragment_length_per_tile_reference_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_gc_for_fragment_length_per_tile_reference_;

		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_base_quality_for_sequence_per_tile_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_base_quality_for_preceding_quality_per_tile_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_base_quality_stats_per_tile_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_preceding_quality_for_sequence_per_tile_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_preceding_quality_for_position_per_tile_;
		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_sequence_quality_for_position_per_tile_;

		std::array<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>, 2> tmp_base_quality_stats_per_strand_;

		std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 5>, 2> tmp_sequence_quality_for_base_per_tile_;
		std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>> tmp_sequence_quality_mean_paired_per_tile_;
		std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>, 2> tmp_sequence_quality_mean_for_gc_per_tile_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_probability_mean_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_minimum_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_first_quartile_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_median_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_third_quartile_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_sequence_quality_maximum_;
		std::array<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>, 2> tmp_sequence_quality_content_;

		std::vector<std::vector<utilities::VectorAtomic<uint64_t>>> tmp_homoquality_distribution_;
		std::array<std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 5>, 2> tmp_nucleotide_quality_;

		// Collected variables for estimation (based on reference)
		std::array<std::array<std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 5>, 4>, 2> base_quality_stats_per_tile_per_error_reference_; // base_quality_stats_per_tile_per_error_reference_[first/second][refBase][domError][tileId][readPosition][quality] = #reads
		std::array<std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 4>, 2> error_rate_for_position_per_tile_per_error_reference_; // error_rate_for_position_per_tile_per_error_reference_[first/second][refBase][domError][tileId][position][errorRate];
		std::array<std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 4>, 2> base_quality_for_error_rate_per_tile_per_error_reference_; // base_quality_for_error_rate_per_tile_per_error_reference_[first/second][refBase][domError][tileId][errorRate][baseQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> base_quality_for_preceding_quality_per_tile_reference_; // base_quality_for_preceding_quality_per_tile_reference_[first/second][refBase][tileId][qualityOfPrecedingBase][baseQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> preceding_quality_for_error_rate_per_tile_reference_; // preceding_quality_for_error_rate_per_tile_reference_[first/second][currentRefBase][tileId][errorRate][precedingQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> preceding_quality_for_position_per_tile_reference_; // preceding_quality_for_position_per_tile_reference_[first/second][currentRefBase][tileId][currentPosition][precedingQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> base_quality_for_sequence_quality_per_tile_reference_; // base_quality_for_sequence_quality_per_tile_reference_[first/second][refBase][tileId][sequenceQuality][baseQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> preceding_quality_for_sequence_quality_per_tile_reference_; // preceding_quality_for_sequence_quality_per_tile_reference_[first/second][refBase][tileId][sequenceQuality][qualityOfPrecedingBase] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> sequence_quality_for_error_rate_per_tile_reference_; // sequence_quality_for_error_rate_per_tile_reference_[first/second][currentRefBase][tileId][errorRate][sequenceQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> sequence_quality_for_position_per_tile_reference_; // sequence_quality_for_position_per_tile_reference_[first/second][currentRefBase][tileId][currentPosition][sequenceQuality] = #bases

		std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 2> sequence_quality_mean_for_gc_per_tile_reference_; // sequence_quality_mean_for_gc_per_tile_reference_[first/second][tileId][percentageGC][qualityMean] = #reads
		std::array<Vect<Vect<Vect<uint64_t>>>, 2> sequence_quality_mean_for_mean_error_rate_per_tile_reference_; // sequence_quality_mean_for_mean_error_rate_per_tile_reference_[first/second][tileId][meanErrorRate][qualityMean] = #reads
		std::array<Vect<Vect<Vect<uint64_t>>>, 2> sequence_quality_mean_for_fragment_length_per_tile_reference_; // sequence_quality_mean_for_fragment_length_per_tile_reference_[first/second][tileId][fragmentLength][qualityMean] = #reads
		std::array<Vect<Vect<Vect<uint64_t>>>, 2> mean_error_rate_for_gc_per_tile_reference_; // mean_error_rate_for_gc_per_tile_reference_[first/second][tileId][percentageGC][meanErrorRate] = #reads
		std::array<Vect<Vect<Vect<uint64_t>>>, 2> mean_error_rate_for_fragment_length_per_tile_reference_; // mean_error_rate_for_fragment_length_per_tile_reference_[first/second][tileId][fragmentLength][meanErrorRate] = #reads
		std::array<Vect<Vect<Vect<uint64_t>>>, 2> gc_for_fragment_length_per_tile_reference_; // gc_for_fragment_length_per_tile_reference_[first/second][tileId][fragmentLength][percentageGC] = #reads

		// Collected variables for plotting (based on raw reads)
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 2> base_quality_for_sequence_per_tile_; // base_quality_for_sequence_per_tile_[first/second][tileId][qualityProbabilityMeanOfSequence][baseQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 2> base_quality_for_preceding_quality_per_tile_; // base_quality_for_preceding_quality_per_tile_[first/second][base][tileId][qualityOfPrecedingBase][baseQuality] = #bases
		std::array<std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 5>, 2> base_quality_stats_per_tile_; // base_quality_stats_per_tile_[first/second][base][tileId][readPosition][quality] = #reads
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 2> preceding_quality_for_sequence_per_tile_; // preceding_quality_for_sequence_per_tile_[first/second][currentBase][tileId][qualityProbabilityMeanOfSequence][precedingQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 2> preceding_quality_for_position_per_tile_; // preceding_quality_for_position_per_tile_[first/second][currentBase][tileId][currentPosition][precedingQuality] = #bases
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 5>, 2> sequence_quality_for_position_per_tile_; // sequence_quality_for_position_per_tile_[first/second][base][tileId][position][sequenceQualityProbabilityMean];

		std::array<Vect<SeqQualityStats<uint64_t>>, 2> base_quality_stats_per_strand_; // base_quality_stats_per_strand_[+/-][readPosition][quality] = #reads

		std::array<std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 5>, 2> sequence_quality_for_base_per_tile_; // sequence_quality_for_base_per_tile_[first/second][nucleotide][tileId][percentageNucleotide][qualityProbabilityMean] = #reads
		Vect<Vect<Vect<uint64_t>>> sequence_quality_mean_paired_per_tile_; // sequence_quality_mean_paired_per_tile_[tileId][qualityMeanFirst][qualityMeanSecond] = #reads
		std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 2> sequence_quality_mean_for_gc_per_tile_; // sequence_quality_mean_for_gc_per_tile_[first/second][tileId][percentageGC][qualityMean] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_probability_mean_; // sequence_quality_probability_mean_[first/second][errorProbabilityMeanAsQuality] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_minimum_; // sequence_quality_minimum_[first/second][qualityMinimum] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_first_quartile_; // sequence_quality_first_quartile_[first/second][qualitySecondQuartile] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_median_; // sequence_quality_median_[first/second][qualityMedian] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_third_quartile_; // sequence_quality_third_quartile_[first/second][qualityThirdQuartile] = #reads
		std::array<Vect<uint64_t>, 2> sequence_quality_maximum_; // sequence_quality_maximum_[first/second][qualityMaximum] = #reads
		std::array<Vect<Vect<uint64_t>>, 2> sequence_quality_content_; // sequence_quality_content_[first/second][quality][#qualityOccurences] = #reads

		Vect<Vect<uint64_t>> homoquality_distribution_; // homoquality_distribution_[quality][length] = #homoqualities
		std::array<std::array<SeqQualityStats<uint64_t>, 5>, 2> nucleotide_quality_; // nucleotide_quality_[first/second][A/C/G/T/N][quality] = #bases

		// Calculated variables for estimation (based on reference)
		std::array<std::array<Vect<Vect<SeqQualityStats<uint64_t>>>, 4>, 2> base_quality_stats_per_tile_reference_; // base_quality_stats_per_tile_reference_[first/second][refBase][tileId][readPosition][quality] = #reads
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> error_rate_for_position_per_tile_reference_; // error_rate_for_position_per_tile_reference_[first/second][refBase][tileId][position][errorRate];
		std::array<std::array<Vect<Vect<Vect<uint64_t>>>, 4>, 2> base_quality_for_error_rate_per_tile_reference_; // base_quality_for_error_rate_per_tile_reference_[first/second][refBase][tileId][errorRate][baseQuality] = #bases

		// Calculated variables for plotting from variables for estimation (based on reference)
		std::array<Vect<SeqQualityStats<uint64_t>>, 2> base_quality_stats_reference_; // base_quality_stats_reference_[first/second][readPosition][quality] = #reads
		std::array<Vect<unsigned char>, 2> base_quality_mean_reference_; // base_quality_mean_reference_[first/second][readPosition] = baseQualityMean
		std::array<Vect<unsigned char>, 2> base_quality_minimum_reference_; // base_quality_minimum_reference_[first/second][readPosition] = baseQualityMinimum
		std::array<Vect<unsigned char>, 2> base_quality_first_quartile_reference_; // base_quality_first_quartile_reference_[first/second][readPosition] = baseQualityFirstQuartile
		std::array<Vect<unsigned char>, 2> base_quality_median_reference_; // base_quality_median_reference_[first/second][readPosition] = baseQualityMedian
		std::array<Vect<unsigned char>, 2> base_quality_third_quartile_reference_; // base_quality_third_quartile_reference_[first/second][readPosition] = baseQualityThirdQuartile
		std::array<Vect<unsigned char>, 2> base_quality_maximum_reference_; // base_quality_maximum_reference_[first/second][readPosition] = baseQualityMaximum

		std::array<Vect<unsigned char>, 2> average_sequence_quality_for_gc_; // average_sequence_quality_for_gc_[first/second][percentageGC] = averageSequenceQualityProbabilityMean

		// Calculated variables for plotting (based on raw reads)
		std::array<Vect<SeqQualityStats<uint64_t>>, 2> base_quality_stats_; // base_quality_stats_[first/second][readPosition][quality] = #reads
		std::array<Vect<unsigned char>, 2> base_quality_mean_; // base_quality_mean_[first/second][readPosition] = baseQualityMean
		std::array<Vect<unsigned char>, 2> base_quality_minimum_; // base_quality_minimum_[first/second][readPosition] = baseQualityMinimum
		std::array<Vect<unsigned char>, 2> base_quality_first_quartile_; // base_quality_first_quartile_[first/second][readPosition] = baseQualityFirstQuartile
		std::array<Vect<unsigned char>, 2> base_quality_median_; // base_quality_median_[first/second][readPosition] = baseQualityMedian
		std::array<Vect<unsigned char>, 2> base_quality_third_quartile_; // base_quality_third_quartile_[first/second][readPosition] = baseQualityThirdQuartile
		std::array<Vect<unsigned char>, 2> base_quality_maximum_; // base_quality_maximum_[first/second][readPosition] = baseQualityMaximum
		std::array<Vect<Vect<signed char>>, 2> tile_quality_mean_difference_; //tile_quality_mean_difference_[first/second][tile][readPosition] = meanDifferenceToTotalQualityMean
		std::array<Vect<unsigned char>, 2> base_quality_mean_per_strand_; // base_quality_mean_per_strand_[+/-][readPosition] = baseQualityMean

		std::array<Vect<Vect<uint64_t>>, 2> base_quality_for_sequence_; // base_quality_for_sequence_[first/second][qualityProbabilityMeanOfSequence][baseQuality] = #bases
		std::array<Vect<Vect<uint64_t>>, 2> base_quality_for_preceding_quality_; // base_quality_for_preceding_quality_[first/second][qualityOfPrecedingBase][baseQuality] = #bases

		std::array<Vect<uint64_t>, 2> sequence_quality_mean_; // sequence_quality_mean_[first/second][qualityMean] = #reads
		std::array<Vect<Vect<uint64_t>>, 2> sequence_quality_mean_per_tile_; // sequence_quality_mean_per_tile_[first/second][tileId][qualityMean] = #reads
		Vect<Vect<uint64_t>> sequence_quality_mean_paired_; // sequence_quality_mean_paired_[qualityMeanFirst][qualityMeanSecond] = #reads
		std::array<Vect<uint16_t>, 2> mean_sequence_quality_mean_by_fragment_length_; // mean_sequence_quality_mean_by_fragment_length_[first/second][fragmentLength] = meanQualityMean
		std::array<std::array<Vect<unsigned char>, 5>, 2> average_sequence_quality_for_base_; // average_sequence_quality_for_base_[first/second][nucleotide][percentageNucleotide] = averageSequenceQualityMean

		// Helper functions
		void SplitPairedSequenceQuality();
		void SumTiles();
		void CalculateQualityStats();
		
		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & base_quality_stats_per_tile_per_error_reference_;
			ar & error_rate_for_position_per_tile_per_error_reference_;
			ar & base_quality_for_error_rate_per_tile_per_error_reference_;
			ar & base_quality_for_preceding_quality_per_tile_reference_;
			ar & preceding_quality_for_error_rate_per_tile_reference_;
			ar & preceding_quality_for_position_per_tile_reference_;
			ar & base_quality_for_sequence_quality_per_tile_reference_;
			ar & preceding_quality_for_sequence_quality_per_tile_reference_;
			ar & sequence_quality_for_error_rate_per_tile_reference_;
			ar & sequence_quality_for_position_per_tile_reference_;

			ar & sequence_quality_mean_for_gc_per_tile_reference_;
			ar & sequence_quality_mean_for_mean_error_rate_per_tile_reference_;
			ar & sequence_quality_mean_for_fragment_length_per_tile_reference_;
			ar & mean_error_rate_for_gc_per_tile_reference_;
			ar & mean_error_rate_for_fragment_length_per_tile_reference_;
			ar & gc_for_fragment_length_per_tile_reference_;

			ar & base_quality_for_sequence_per_tile_;
			ar & base_quality_for_preceding_quality_per_tile_;
			ar & base_quality_stats_per_tile_;
			ar & preceding_quality_for_sequence_per_tile_;
			ar & preceding_quality_for_position_per_tile_;
			ar & sequence_quality_for_position_per_tile_;

			ar & base_quality_stats_per_strand_;

			ar & sequence_quality_for_base_per_tile_;
			ar & sequence_quality_mean_paired_per_tile_;
			ar & sequence_quality_mean_for_gc_per_tile_;
			ar & sequence_quality_probability_mean_;
			ar & sequence_quality_minimum_;
			ar & sequence_quality_first_quartile_;
			ar & sequence_quality_median_;
			ar & sequence_quality_third_quartile_;
			ar & sequence_quality_maximum_;
			ar & sequence_quality_content_;

			ar & homoquality_distribution_;
			ar & nucleotide_quality_;
		}

		// Google test
		friend class QualityStatsTest;
		friend class ProbabilityEstimatesTest;

	public:
		QualityStats(){
			for( uint16_t template_segment=2; template_segment--; ){
				// There are no quality values below 2 for Illumina
				sequence_quality_mean_paired_per_tile_.SetOffset(2);
				sequence_quality_probability_mean_.at(template_segment).SetOffset(2);
				sequence_quality_minimum_.at(template_segment).SetOffset(2);
				sequence_quality_first_quartile_.at(template_segment).SetOffset(2);
				sequence_quality_median_.at(template_segment).SetOffset(2);
				sequence_quality_third_quartile_.at(template_segment).SetOffset(2);
				sequence_quality_maximum_.at(template_segment).SetOffset(2);
				sequence_quality_content_.at(template_segment).SetOffset(2);

				for(auto called_base = 5; called_base--; ){
					nucleotide_quality_.at(template_segment).at(called_base).SetOffset(2);
				}
			}
		}

		// Getter functions
		inline const Vect<SeqQualityStats<uint64_t>> &BaseQualityStatsReference(uint16_t template_segment, uint16_t tile_id, uint16_t ref_base, uint16_t dom_error) const{
			return base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &ErrorRateForPositionReference(uint16_t template_segment, uint16_t tile_id, uint16_t base, uint16_t dom_error) const{
			return error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(base).at(dom_error).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &BaseQualityForErrorRateReference(uint16_t template_segment, uint16_t tile_id, uint16_t base, uint16_t dom_error) const{
			return base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(base).at(dom_error)[tile_id];
		}
		inline const Vect<Vect<uint64_t>> &BaseQualityForPrecedingQualityReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &PrecedingQualityForErrorRateReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &PrecedingQualityForPositionReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return preceding_quality_for_position_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &BaseQualityForSequenceQualityReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &PrecedingQualityForSequenceQualityReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityForErrorRateReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityForPositionReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return sequence_quality_for_position_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}

		inline const Vect<SeqQualityStats<uint64_t>> &SequenceQualityMeanForGCPerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityMeanForMeanErrorRatePerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityMeanForFragmentLengthPerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &MeanErrorRateForGCPerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return mean_error_rate_for_gc_per_tile_reference_.at(template_segment).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &MeanErrorRateForFragmentLengthPerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &GCForFragmentLengthPerTileReference(uint16_t template_segment, uint16_t tile_id) const{
			return gc_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id);
		}

		inline const Vect<Vect<uint64_t>> &SequenceQualityMeanPaired(uint16_t tile_id) const{
			return sequence_quality_mean_paired_per_tile_.at(tile_id);
		}
		inline const Vect<uint64_t> &SequenceQualityProbabilityMean(uint16_t template_segment) const{
			return sequence_quality_probability_mean_.at(template_segment);
		}
		inline const Vect<uint64_t> &SequenceQualityMinimum(uint16_t template_segment) const{
			return sequence_quality_minimum_.at(template_segment);
		}
		inline const Vect<uint64_t> &SequenceQualityFirstQuartile(uint16_t template_segment) const{
			return sequence_quality_first_quartile_.at(template_segment);
		}
		inline const Vect<uint64_t> &SequenceQualityMedian(uint16_t template_segment) const{
			return sequence_quality_median_.at(template_segment);
		}
		inline const Vect<uint64_t> &SequenceQualityThirdQuartile(uint16_t template_segment) const{
			return sequence_quality_third_quartile_.at(template_segment);
		}
		inline const Vect<uint64_t> &SequenceQualityMaximum(uint16_t template_segment) const{
			return sequence_quality_maximum_.at(template_segment);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityContent(uint16_t template_segment) const{
			return sequence_quality_content_.at(template_segment);
		}

		inline const Vect<Vect<uint64_t>> &HomoqualityDistribution() const{
			return homoquality_distribution_;
		}
		inline const SeqQualityStats<uint64_t> &NucleotideQuality(uint16_t template_segment, uint16_t nucleotide) const{
			return nucleotide_quality_.at(template_segment).at(nucleotide);
		}

		inline const Vect<SeqQualityStats<uint64_t>> &BaseQualityStatsReference(uint16_t template_segment, uint16_t tile_id, uint16_t ref_base) const{
			return base_quality_stats_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &ErrorRateForPositionReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return error_rate_for_position_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}
		inline const Vect<Vect<uint64_t>> &BaseQualityForErrorRateReference(uint16_t template_segment, uint16_t tile_id, uint16_t base) const{
			return base_quality_for_error_rate_per_tile_reference_.at(template_segment).at(base).at(tile_id);
		}

		inline const Vect<SeqQualityStats<uint64_t>> &BaseQualityStatsReference(uint16_t template_segment) const{
			return base_quality_stats_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMeanReference(uint16_t template_segment) const{
			return base_quality_mean_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMinimumReference(uint16_t template_segment) const{
			return base_quality_minimum_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityFirstQuartileReference(uint16_t template_segment) const{
			return base_quality_first_quartile_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMedianReference(uint16_t template_segment) const{
			return base_quality_median_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityThirdQuartileReference(uint16_t template_segment) const{
			return base_quality_third_quartile_reference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMaximumReference(uint16_t template_segment) const{
			return base_quality_maximum_reference_.at(template_segment);
		}

		inline const Vect<unsigned char> &AverageSequenceQualityForGC(uint16_t template_segment) const{
			return average_sequence_quality_for_gc_.at(template_segment);
		}

		inline const Vect<SeqQualityStats<uint64_t>> &BaseQualityStats(uint16_t template_segment) const{
			return base_quality_stats_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMean(uint16_t template_segment) const{
			return base_quality_mean_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMinimum(uint16_t template_segment) const{
			return base_quality_minimum_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityFirstQuartile(uint16_t template_segment) const{
			return base_quality_first_quartile_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMedian(uint16_t template_segment) const{
			return base_quality_median_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityThirdQuartile(uint16_t template_segment) const{
			return base_quality_third_quartile_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMaximum(uint16_t template_segment) const{
			return base_quality_maximum_.at(template_segment);
		}
		inline const Vect<Vect<signed char>> &TileQualityMeanDifference(uint16_t template_segment) const{
			return tile_quality_mean_difference_.at(template_segment);
		}
		inline const Vect<unsigned char> &BaseQualityMeanPerStrand(uint16_t strand) const{
			return base_quality_mean_per_strand_.at(strand);
		}

		inline const Vect<Vect<uint64_t>>  &BaseQualityForSequence(uint16_t template_segment) const{
			return base_quality_for_sequence_.at(template_segment);
		}
		inline const Vect<Vect<uint64_t>>  &BaseQualityForPrecedingQuality(uint16_t template_segment) const{
			return base_quality_for_preceding_quality_.at(template_segment);
		}

		inline const Vect<uint64_t> &SequenceQualityMean(uint16_t template_segment) const{
			return sequence_quality_mean_.at(template_segment);
		}
		inline const Vect<Vect<uint64_t>> &SequenceQualityMeanPaired() const{
			return sequence_quality_mean_paired_;
		}
		inline const Vect<uint16_t> &MeanSequenceQualityMeanByFragmentLength(uint16_t template_segment) const{
			return mean_sequence_quality_mean_by_fragment_length_.at(template_segment);
		}
		inline const Vect<unsigned char> &AverageSequenceQualityForBase(uint16_t template_segment, uint16_t nucleotide) const{
			return average_sequence_quality_for_base_.at(template_segment).at(nucleotide);
		}
		
		// Setter functions
		inline void AddRefBase(
				uint16_t template_segment,
				uint16_t ref_base,
				uint16_t dom_error,
				uint16_t tile_id,
				uint16_t quality,
				uint16_t error_rate,
				uint16_t last_qual,
				uint16_t seq_qual,
				uint32_t read_pos){
			// For estimation
			++tmp_base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(read_pos).at(quality);
			++tmp_error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(read_pos).at(error_rate);
			++tmp_base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(error_rate).at(quality);
			++tmp_base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(last_qual).at(quality);
			++tmp_preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(error_rate).at(last_qual);
			++tmp_preceding_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(read_pos).at(last_qual);
			++tmp_base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(seq_qual).at(quality);
			++tmp_preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(seq_qual).at(last_qual);
			++tmp_sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(error_rate).at(seq_qual);
			++tmp_sequence_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base).at(tile_id).at(read_pos).at(seq_qual);
		}
		inline void AddRefRead(uint16_t template_segment, uint16_t tile_id, uint16_t gc, uint16_t seq_qual_mean, uint16_t mean_error_rate, uint32_t fragment_length){
			++tmp_sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment).at(tile_id).at(gc).at(seq_qual_mean);
			++tmp_sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment).at(tile_id).at(mean_error_rate).at(seq_qual_mean);
			++tmp_sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id).at(fragment_length).at(seq_qual_mean);
			++tmp_mean_error_rate_for_gc_per_tile_reference_.at(template_segment).at(tile_id).at(gc).at(mean_error_rate);
			++tmp_mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id).at(fragment_length).at(mean_error_rate);
			++tmp_gc_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id).at(fragment_length).at(gc);
		}

		inline void AddRawBase(
				uint16_t template_segment,
				uint16_t called_base,
				uint16_t tile_id,
				uint16_t strand,
				uint16_t quality,
				uint16_t seq_qual,
				uint16_t last_qual,
				uint32_t read_pos){
			// For plotting
			++tmp_base_quality_for_sequence_per_tile_.at(template_segment).at(called_base).at(tile_id).at(seq_qual).at(quality);
			++tmp_base_quality_stats_per_tile_.at(template_segment).at(called_base).at(tile_id).at(read_pos).at(quality);
			++tmp_sequence_quality_for_position_per_tile_.at(template_segment).at(called_base).at(tile_id).at(read_pos).at(seq_qual);
			++tmp_base_quality_for_preceding_quality_per_tile_.at(template_segment).at(called_base).at(tile_id).at(last_qual).at(quality);
			++tmp_preceding_quality_for_sequence_per_tile_.at(template_segment).at(called_base).at(tile_id).at(seq_qual).at(last_qual);
			++tmp_preceding_quality_for_position_per_tile_.at(template_segment).at(called_base).at(tile_id).at(read_pos).at(last_qual);

			++tmp_base_quality_stats_per_strand_.at(strand).at(read_pos).at(quality);
			++tmp_nucleotide_quality_.at(template_segment).at(called_base).at(quality);
		}
		inline void AddRawHomoqualimer(uint16_t quality, uint32_t length){
			++tmp_homoquality_distribution_.at(quality).at(length);
		}
		inline void AddRawRead(
				unsigned char &paired_seq_qual,
				const SeqQualityStats<uint64_t> &seq_qual_stats,
				uint16_t template_segment,
				uint16_t tile_id,
				std::array<uint64_t, 5> &read_bases,
				uint32_t read_length){
			for( auto qual = seq_qual_stats.from(); qual < seq_qual_stats.to(); ++qual){
				if(seq_qual_stats.at(qual)){ // So that qualities that do not exist are not filled in
					++tmp_sequence_quality_content_.at(template_segment).at(qual).at(seq_qual_stats.at(qual));
				}
			}
			++tmp_sequence_quality_minimum_.at(template_segment).at(seq_qual_stats.minimum_);
			++tmp_sequence_quality_first_quartile_.at(template_segment).at(seq_qual_stats.first_quartile_);
			++tmp_sequence_quality_median_.at(template_segment).at(seq_qual_stats.median_);
			++tmp_sequence_quality_third_quartile_.at(template_segment).at(seq_qual_stats.third_quartile_);
			++tmp_sequence_quality_maximum_.at(template_segment).at(seq_qual_stats.maximum_);

			if( paired_seq_qual ){
				if(template_segment){
					++tmp_sequence_quality_mean_paired_per_tile_.at(tile_id).at(paired_seq_qual).at(seq_qual_stats.mean_);
				}
				else{
					++tmp_sequence_quality_mean_paired_per_tile_.at(tile_id).at(seq_qual_stats.mean_).at(paired_seq_qual);
				}
			}
			else{
				paired_seq_qual = seq_qual_stats.mean_;
			}

			for(auto nuc=5; nuc--; ){
				++tmp_sequence_quality_for_base_per_tile_.at(template_segment).at(nuc).at(tile_id).at(reseq::utilities::Percent( read_bases.at(nuc), read_length )).at(seq_qual_stats.mean_);
			}
			++tmp_sequence_quality_mean_for_gc_per_tile_.at(template_segment).at(tile_id).at(reseq::utilities::Percent( read_bases.at(1)+read_bases.at(2)+read_bases.at(4)/2, read_length )).at(seq_qual_stats.mean_);

			++tmp_sequence_quality_probability_mean_.at(template_segment).at(seq_qual_stats.probability_mean_);
		}

		// Main functions
		void Prepare(uint16_t num_tiles, uint8_t size_qual, uint32_t size_pos, uint32_t maximum_fragment_length);
		void Finalize(uint64_t total_number_reads);
		void Shrink();
		void PrepareEstimation();
		void PreparePlotting();
		void PrepareTesting();
	};
}

#endif // QUALITYSTATS_H
