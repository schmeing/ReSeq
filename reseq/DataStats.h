#ifndef DATASTATS_H
#define DATASTATS_H

#include <atomic>
#include <array>
#include <condition_variable>
#include <chrono>
#include <map>
#include <mutex>
#include <stdint.h>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <seqan/bam_io.h>

#include "AdapterStats.h"
#include "CoverageStats.h"
#include "ErrorStats.h"
#include "FragmentDistributionStats.h"
#include "FragmentDuplicationStats.h"
#include "QualityStats.h"
#include "Reference.h"
#include "SeqQualityStats.hpp"
#include "TileStats.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	class DataStats{
	public:
		// Operator Definitions
		struct RecordHasher{
			size_t operator() (CoverageStats::FullRecord * const &record) const;
		};
		struct RecordEqual{
			bool operator() (CoverageStats::FullRecord * const &lhs, CoverageStats::FullRecord * const &rhs) const;
		};

	private:
		struct ThreadData{
			std::vector< std::pair<CoverageStats::FullRecord *, CoverageStats::FullRecord *> > rec_store_;

			CoverageStats::ThreadData coverage_;
			FragmentDistributionStats::ThreadData fragment_distribution_;

			uintSeqLen last_exclusion_region_id_;
			uintRefSeqId last_exclusion_ref_seq_;

			ThreadData( uintSeqLen maximum_insert_length, uintSeqLen max_seq_bin_len, uintSeqLen start_exclusion_length ):
				fragment_distribution_(maximum_insert_length, max_seq_bin_len, start_exclusion_length),
				last_exclusion_region_id_(0),
				last_exclusion_ref_seq_(0)
			{}
		};

		// Definitions
		const uintFragCount kBatchSize = 10000; // Number of reads after which the progress is reported

		// User parameter
		Reference *reference_;

		const uintSeqLen maximum_insert_length_;
		const uintQual minimum_mapping_quality_; // Minimum mapping quality required for both reads to use read pair for insert_lengths_, fragment_duplication_number_, ...

		// Mutex
		std::mutex read_mutex_;
		std::mutex print_mutex_;

		std::atomic<uintNumThreads> running_threads_;
		std::mutex finish_threads_mutex_;
		std::condition_variable finish_threads_cv_;
		bool finish_threads_;

		// Subclasses
		AdapterStats adapters_;
		CoverageStats coverage_;
		FragmentDuplicationStats duplicates_;
		ErrorStats errors_;
		FragmentDistributionStats fragment_distribution_;
		QualityStats qualities_;
		TileStats tiles_;

		// Temporary variables
		std::unordered_set<CoverageStats::FullRecord *, RecordHasher, RecordEqual> first_read_records_;
		std::atomic<bool> reading_success_;
		uintFragCount read_records_;
		
		std::array<std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>>, 2> tmp_read_lengths_by_fragment_length_;
		std::array<std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>>, 2> tmp_non_mapped_read_lengths_by_fragment_length_;

		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_proper_pair_mapping_quality_;
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_improper_pair_mapping_quality_;
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_single_read_mapping_quality_;

		std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_gc_read_content_;
		std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_gc_read_content_reference_;
		std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_gc_read_content_mapped_;

		std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 2> tmp_n_content_;
		std::array<std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5>, 2> tmp_sequence_content_;
		std::array<std::array<std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 4>, 2>, 2> tmp_sequence_content_reference_;

		std::array<std::vector<utilities::VectorAtomic<uintNucCount>>, 5> tmp_homopolymer_distribution_;

		std::vector<uintFragCount> reads_per_frag_len_bin_; // reads_per_frag_len_bin_[BinOfReferenceSequenceBinnedInBinsOfFragmentLength] = #Reads

		// Collected variables for simulation
		uint64_t creation_time_; // Store time when bam file was completelly read, can be used to check whether the stats file was updated
		std::array<Vect<uintFragCount>, 2> read_lengths_; // read_lengths_[first/second][length] = #reads
		std::array<Vect<Vect<uintFragCount>>, 2> read_lengths_by_fragment_length_; // read_lengths_by_fragment_length_[first/second][fragment_length][read_length] = #reads
		std::array<Vect<Vect<uintFragCount>>, 2> non_mapped_read_lengths_by_fragment_length_; // non_mapped_read_lengths_by_fragment_length_[first/second][fragment_length][read_length] = #reads
		uintQual phred_quality_offset_;
		uintQual minimum_quality_;
		uintQual maximum_quality_;
		uintReadLen minimum_read_length_on_reference_;
		uintReadLen maximum_read_length_on_reference_;
		double percentage_high_enough_quality_reads_;

		// Collected variables for plotting
		Vect<uintFragCount> proper_pair_mapping_quality_; // proper_pair_mapping_quality_[mappingQuality] = #readsWithProperPairFlag
		Vect<uintFragCount> improper_pair_mapping_quality_; // improper_pair_mapping_quality_[mappingQuality] = #readsWithPartnerMappedButWithoutProperPairFlag
		Vect<uintFragCount> single_read_mapping_quality_; // single_read_mapping_quality_[mappingQuality] = #readsWithPartnerNotMapped

		std::array<Vect<uintFragCount>, 2> gc_read_content_; // gc_read_content_[first/second][gcContent(%)] = #reads
		std::array<Vect<uintFragCount>, 2> gc_read_content_reference_; // gc_read_content_reference_[first/second][gcContentReference(%)] = #reads
		std::array<Vect<uintFragCount>, 2> gc_read_content_mapped_; // gc_read_content_mapped_[first/second][gcContent(%)] = #reads

		std::array<Vect<uintFragCount>, 2> n_content_; // n_content_[first/second][nContent(%)] = #reads
		std::array<std::array<Vect<uintNucCount>, 5>, 2> sequence_content_; // sequence_content_[first/second][A/C/G/T/N][readPosition] = #(reads with given content at given position)
		std::array<std::array<std::array<Vect<uintNucCount>, 4>, 2>, 2> sequence_content_reference_; // sequence_content_reference_[first/second][forward/reverse][A/C/G/T][readPosition] = #(reads with given reference content at given reference position)

		std::array<Vect<uintNucCount>, 5> homopolymer_distribution_; // homopolymer_distribution_[A/C/G/T/N][length] = #homopolymers
		
		// Collected variables for output
		uintFragCount total_number_reads_;
		std::atomic<uintFragCount> reads_in_unmapped_pairs_without_adapters_;
		std::atomic<uintFragCount> reads_in_unmapped_pairs_with_adapters_;
		std::atomic<uintFragCount> reads_with_low_quality_with_adapters_;
		std::atomic<uintFragCount> reads_with_low_quality_without_adapters_;
		std::atomic<uintFragCount> reads_on_too_short_fragments_;
		std::atomic<uintFragCount> reads_in_excluded_regions_;
		std::atomic<uintFragCount> reads_used_;

		// Private functions
		inline bool PotentiallyValidGeneral( const seqan::BamAlignmentRecord &record ) const;
		inline bool PotentiallyValidFirst( const seqan::BamAlignmentRecord &record_first ) const;
		inline bool PotentiallyValidSecond( const seqan::BamAlignmentRecord &record_second ) const;
		inline bool PotentiallyValid( const seqan::BamAlignmentRecord &record ) const;
		bool IsSecondRead( CoverageStats::FullRecord *record, CoverageStats::FullRecord *&record_first, CoverageStats::CoverageBlock *&block );

		bool CheckForAdapters(const seqan::BamAlignmentRecord &record_first, const seqan::BamAlignmentRecord &record_second);
		void EvalBaseLevelStats( CoverageStats::FullRecord *full_record, uintTempSeq template_segment, uintTempSeq strand, uintTileId tile_id, uintQual &paired_seq_qual );
		bool EvalReferenceStatistics( CoverageStats::FullRecord *record, uintTempSeq template_segment, CoverageStats::CoverageBlock *coverage_block );
		bool EvalRecord( std::pair<CoverageStats::FullRecord *, CoverageStats::FullRecord *> record, ThreadData &thread ); // Don't use a reference hear, so it isn't affected by a pointer switch

		bool SignsOfPairsWithNamesNotIdentical();
		void PrepareReadIn(uintQual size_mapping_quality, uintReadLen size_indel, uintSeqLen max_ref_seq_bin_size);
		bool FinishReadIn();
		void Shrink(); // Reduce all arrays to minimal size by removing unused bins at the beginning and end
		bool Calculate(uintNumThreads num_threads);
		void PrepareGeneral();

		bool OrderOfBamFileCorrect( const seqan::BamAlignmentRecord &record, std::pair< uintRefSeqId, uintSeqLen > last_record_pos );
		bool PreRun(seqan::BamFileIn &bam, const char *bam_file, seqan::BamHeader &header, uintQual &size_mapping_quality, uintReadLen &size_indel);
		bool ReadRecords( seqan::BamFileIn &bam, bool &not_done, ThreadData &thread_data );

		static void ReadThread( DataStats &self, seqan::BamFileIn &bam, uintSeqLen max_seq_bin_len );

		// boost serialization
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
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
			ar & percentage_high_enough_quality_reads_;

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
		}

		// Google test
		friend class DataStatsTest;
		friend class ProbabilityEstimatesTest;
		friend class SimulatorTest;
		FRIEND_TEST(DataStatsTest, Construction);

	public:
		DataStats(Reference *ref, uintSeqLen maximum_insert_length=2000, uintQual minimum_mapping_quality=10);

		// Getter functions
		const AdapterStats &Adapters() const{ return adapters_; }
		const CoverageStats &Coverage() const{ return coverage_; }
		const ErrorStats &Errors() const{ return errors_; }
		const FragmentDuplicationStats &Duplicates() const{ return duplicates_; }
		const FragmentDistributionStats &FragmentDistribution() const{ return fragment_distribution_; }
		FragmentDistributionStats &FragmentDistribution(){ return fragment_distribution_; }
		const QualityStats &Qualities() const{ return qualities_; }
		const TileStats &Tiles() const{ return tiles_; }

		uintQual PhredQualityOffset() const{ return phred_quality_offset_; }
		uintFragCount TotalNumberReads() const{ return total_number_reads_; }

		uint64_t CreationTime() const{
			return creation_time_;
		}

		const Vect<uintFragCount> &ReadLengths(uintTempSeq template_segment) const{ return read_lengths_.at(template_segment); }
		const Vect<Vect<uintFragCount>> &ReadLengthsByFragmentLength(uintTempSeq template_segment) const{ return read_lengths_by_fragment_length_.at(template_segment); }
		const Vect<Vect<uintFragCount>> &NonMappedReadLengthsByFragmentLength(uintTempSeq template_segment) const{ return non_mapped_read_lengths_by_fragment_length_.at(template_segment); }

		double PercentageHighEnoughQualityReads() const{return percentage_high_enough_quality_reads_; }
		const Vect<uintFragCount> &ProperPairMappingQuality() const{ return proper_pair_mapping_quality_; }
		const Vect<uintFragCount> &ImproperPairMappingQuality() const{ return improper_pair_mapping_quality_; }
		const Vect<uintFragCount> &SingleReadMappingQuality() const{ return single_read_mapping_quality_; }

		const Vect<uintFragCount> &GCReadContent(uintTempSeq template_segment) const{ return gc_read_content_.at(template_segment); }
		const Vect<uintFragCount> &GCReadContentReference(uintTempSeq template_segment) const{ return gc_read_content_reference_.at(template_segment); }
		const Vect<uintFragCount> &GCReadContentMapped(uintTempSeq template_segment) const{ return gc_read_content_mapped_.at(template_segment); }

		const Vect<uintFragCount> &NContent(uintTempSeq template_segment) const{ return n_content_.at(template_segment); }
		const Vect<uintNucCount> &SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const{ return sequence_content_.at(template_segment).at(nucleotide); }
		const Vect<uintNucCount> &SequenceContentReference(uintTempSeq template_segment, uintTempSeq strand, uintBaseCall nucleotide) const{ return sequence_content_reference_.at(template_segment).at(strand).at(nucleotide); }

		const Vect<uintNucCount> &HomopolymerDistribution(uintBaseCall nucleotide) const{ return homopolymer_distribution_.at(nucleotide); }

		bool HasReference() const{ return reference_; }

		// Setter functions
		void IgnoreTiles(){ tiles_.IgnoreTiles(); }
		void ClearReference(){ reference_ = NULL; }
		void SetUniformBias(){
			fragment_distribution_.SetUniformBias();
		}

		// Public functions
		bool IsValidRecord( const seqan::BamAlignmentRecord &record );
		inline bool QualitySufficient( const seqan::BamAlignmentRecord &record ) const{
			return minimum_mapping_quality_ <= record.mapQ;
		}
		static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord &record);
		static uintReadLen GetReadLengthOnReference(const seqan::BamAlignmentRecord &record, uintReadLen &max_indel);
		inline void GetReadPosOnReference(uintSeqLen &start_pos, uintSeqLen &end_pos, const seqan::BamAlignmentRecord &record) const;
		inline bool InProperDirection( const seqan::BamAlignmentRecord &record_first, uintSeqLen end_pos_second, uintSeqLen start_pos_first, uintSeqLen start_pos_second  ) const{
			// Proper forward reverse direction or read sized pair (proper direction is only checked if they are on same scaffold, so no check needed for that)
			return (!hasFlagRC(record_first) && hasFlagNextRC(record_first) && maximum_insert_length_ >= end_pos_second - start_pos_first) || start_pos_first == start_pos_second;
		}

		bool ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const std::string &variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads, bool calculate_bias = true ); // Fill the class with the information from a bam file

		bool Load( const char *archive_file );
		bool Save( const char *archive_file ) const;

		void PrepareProcessing();
		void PreparePlotting();
		void PrepareTesting();
	};

}

#endif // DATASTATS_H
