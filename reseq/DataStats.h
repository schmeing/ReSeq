#ifndef DATASTATS_H
#define DATASTATS_H

#include <atomic>
#include <array>
#include <condition_variable>
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
#include "utilities.h"
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

			FragmentDistributionStats::ThreadData fragment_distribution_;

			ThreadData( uint32_t maximum_insert_length, uint32_t max_seq_bin_len ):
				fragment_distribution_(maximum_insert_length, max_seq_bin_len)
			{}
		};

		// Definitions
		Reference *reference_;

		const uint32_t batch_size_; // Number of reads after which the progress is reported
		const uint32_t maximum_insert_length_;
		const uint16_t minimum_mapping_quality_; // Minimum mapping quality required for both reads to use read pair for insert_lengths_, fragment_duplication_number_, ...
		const uint64_t block_size_;
		uint16_t min_dist_to_ref_seq_ends_; // Is set in Reference class

		// Mutex
		std::mutex read_mutex_;
		std::mutex print_mutex_;

		std::atomic<uint16_t> running_threads_;
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
		uint64_t read_records_;
		
		std::array<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>, 2> tmp_read_lengths_by_fragment_length_;

		std::vector<utilities::VectorAtomic<uint64_t>> tmp_proper_pair_mapping_quality_;
		std::vector<utilities::VectorAtomic<uint64_t>> tmp_improper_pair_mapping_quality_;
		std::vector<utilities::VectorAtomic<uint64_t>> tmp_single_read_mapping_quality_;

		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_gc_read_content_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_gc_read_content_reference_;
		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_gc_read_content_mapped_;

		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 2> tmp_n_content_;
		std::array<std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 5>, 2> tmp_sequence_content_;
		std::array<std::array<std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 4>, 2>, 2> tmp_sequence_content_reference_;

		std::array<std::vector<utilities::VectorAtomic<uint64_t>>, 5> tmp_homopolymer_distribution_;

		std::vector<uint64_t> reads_per_ref_seq_bin_; // reads_per_ref_seq_bin_[ReferenceSequenceBin] = #Reads

		// Collected variables for simulation
		std::array<Vect<uint64_t>, 2> read_lengths_; // read_lengths_[first/second][length] = #reads
		std::array<Vect<Vect<uint64_t>>, 2> read_lengths_by_fragment_length_; // read_lengths_by_fragment_length_[first/second][fragment_length][read_length] = #reads
		uint8_t phred_quality_offset_;
		uint8_t minimum_quality_;
		uint8_t maximum_quality_;
		uint32_t minimum_read_length_on_reference_;
		uint32_t maximum_read_length_on_reference_;

		// Collected variables for plotting
		Vect<uint64_t> proper_pair_mapping_quality_; // proper_pair_mapping_quality_[mappingQuality] = #readsWithProperPairFlag
		Vect<uint64_t> improper_pair_mapping_quality_; // improper_pair_mapping_quality_[mappingQuality] = #readsWithPartnerMappedButWithoutProperPairFlag
		Vect<uint64_t> single_read_mapping_quality_; // single_read_mapping_quality_[mappingQuality] = #readsWithPartnerNotMapped

		std::array<Vect<uint64_t>, 2> gc_read_content_; // gc_read_content_[first/second][gcContent(%)] = #reads
		std::array<Vect<uint64_t>, 2> gc_read_content_reference_; // gc_read_content_reference_[first/second][gcContentReference(%)] = #reads
		std::array<Vect<uint64_t>, 2> gc_read_content_mapped_; // gc_read_content_mapped_[first/second][gcContent(%)] = #reads

		std::array<Vect<uint64_t>, 2> n_content_; // n_content_[first/second][nContent(%)] = #reads
		std::array<std::array<Vect<uint64_t>, 5>, 2> sequence_content_; // sequence_content_[first/second][A/C/G/T/N][readPosition] = #(reads with given content at given position)
		std::array<std::array<std::array<Vect<uint64_t>, 4>, 2>, 2> sequence_content_reference_; // sequence_content_reference_[first/second][forward/reverse][A/C/G/T][readPosition] = #(reads with given reference content at given reference position)

		std::array<Vect<uint64_t>, 5> homopolymer_distribution_; // homopolymer_distribution_[A/C/G/T/N][length] = #homopolymers
		
		// Collected variables for output
		uint64_t total_number_reads_;
		std::atomic<uint64_t> reads_in_unmapped_pairs_without_adapters_;
		std::atomic<uint64_t> reads_in_unmapped_pairs_with_adapters_;
		std::atomic<uint64_t> reads_with_low_quality_;
		std::atomic<uint64_t> reads_on_too_short_fragments_;
		std::atomic<uint64_t> reads_to_close_to_ref_seq_ends_;
		std::atomic<uint64_t> reads_used_;

		// Private functions
		inline bool PotentiallyValid( const seqan::BamAlignmentRecord &record_first ) const;
		bool IsSecondRead( CoverageStats::FullRecord *record, CoverageStats::FullRecord *&record_first, CoverageStats::CoverageBlock *&block );

		bool CheckForAdapters(const seqan::BamAlignmentRecord &record_first, const seqan::BamAlignmentRecord &record_second);
		void EvalBaseLevelStats( CoverageStats::FullRecord *full_record, uint16_t template_segment, uint16_t tile_id, uint8_t &paired_seq_qual );
		bool EvalReferenceStatistics( CoverageStats::FullRecord *record, uint16_t template_segment, CoverageStats::CoverageBlock *coverage_block, ThreadData& thread_data, uint32_t paired_read_end );
		bool EvalRecord( std::pair<CoverageStats::FullRecord *, CoverageStats::FullRecord *> &record, ThreadData &thread_data );

		bool SignsOfPairsWithNamesNotIdentical();
		void PrepareReadIn(uint8_t size_mapping_quality, uint32_t size_indel);
		bool FinishReadIn();
		void Shrink(); // Reduce all arrays to minimal size by removing unused bins at the beginning and end
		bool Calculate(uint16_t num_threads);
		void PrepareGeneral();

		bool OrderOfBamFileCorrect( const seqan::BamAlignmentRecord &record, std::pair< uint32_t, uint32_t > last_record_pos );
		bool PreRun(seqan::BamFileIn &bam, const char *bam_file, seqan::BamHeader &header, uint8_t &size_mapping_quality, uint32_t &size_indel, uint64_t max_ref_seq_bin_size);
		bool ReadRecords( seqan::BamFileIn &bam, bool &not_done, ThreadData &thread_data );

		static void ReadThread( DataStats &self, seqan::BamFileIn &bam, uint32_t max_seq_bin_len );

		// boost serialization
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & min_dist_to_ref_seq_ends_;

			ar & adapters_;
			ar & coverage_;
			ar & errors_;
			ar & duplicates_;
			ar & fragment_distribution_;
			ar & qualities_;
			ar & tiles_;

			ar & read_lengths_;
			ar & read_lengths_by_fragment_length_;
			ar & phred_quality_offset_;
			ar & minimum_quality_;
			ar & maximum_quality_;
			ar & minimum_read_length_on_reference_;
			ar & maximum_read_length_on_reference_;

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
		DataStats(Reference *ref, uint32_t maximum_insert_length=2000, uint16_t minimum_mapping_quality=10);

		// Getter functions
		const AdapterStats &Adapters() const{ return adapters_; }
		const CoverageStats &Coverage() const{ return coverage_; }
		const ErrorStats &Errors() const{ return errors_; }
		const FragmentDuplicationStats &Duplicates() const{ return duplicates_; }
		const FragmentDistributionStats &FragmentDistribution() const{ return fragment_distribution_; }
		FragmentDistributionStats &FragmentDistribution(){ return fragment_distribution_; }
		const QualityStats &Qualities() const{ return qualities_; }
		const TileStats &Tiles() const{ return tiles_; }

		uint8_t PhredQualityOffset() const{ return phred_quality_offset_; }
		uint64_t TotalNumberReads() const{ return total_number_reads_; }

		const Vect<uint64_t> &ReadLengths(uint16_t template_segment) const{ return read_lengths_.at(template_segment); }
		const Vect<Vect<uint64_t>> &ReadLengthsByFragmentLength(uint16_t template_segment) const{ return read_lengths_by_fragment_length_.at(template_segment); }

		const Vect<uint64_t> &ProperPairMappingQuality() const{ return proper_pair_mapping_quality_; }
		const Vect<uint64_t> &ImproperPairMappingQuality() const{ return improper_pair_mapping_quality_; }
		const Vect<uint64_t> &SingleReadMappingQuality() const{ return single_read_mapping_quality_; }

		const Vect<uint64_t> &GCReadContent(uint16_t template_segment) const{ return gc_read_content_.at(template_segment); }
		const Vect<uint64_t> &GCReadContentReference(uint16_t template_segment) const{ return gc_read_content_reference_.at(template_segment); }
		const Vect<uint64_t> &GCReadContentMapped(uint16_t template_segment) const{ return gc_read_content_mapped_.at(template_segment); }

		const Vect<uint64_t> &NContent(uint16_t template_segment) const{ return n_content_.at(template_segment); }
		const Vect<uint64_t> &SequenceContent(uint16_t template_segment, uint16_t nucleotide) const{ return sequence_content_.at(template_segment).at(nucleotide); }
		const Vect<uint64_t> &SequenceContentReference(uint16_t template_segment, uint16_t strand, uint16_t nucleotide) const{ return sequence_content_reference_.at(template_segment).at(strand).at(nucleotide); }

		const Vect<uint64_t> &HomopolymerDistribution(uint16_t nucleotide) const{ return homopolymer_distribution_.at(nucleotide); }

		bool HasReference() const{ return reference_; }

		// Setter functions
		void IgnoreTiles(){ tiles_.IgnoreTiles(); }
		void ClearReference(){ reference_ = NULL; }

		// Public functions
		bool IsValidRecord( const seqan::BamAlignmentRecord &record );
		inline bool QualitySufficient( const seqan::BamAlignmentRecord &record ) const{
			return minimum_mapping_quality_ <= record.mapQ;
		}
		static uint32_t GetReadLengthOnReference(const seqan::BamAlignmentRecord &record);
		static uint32_t GetReadLengthOnReference(const seqan::BamAlignmentRecord &record, uint32_t &max_indel);
		inline bool InProperDirection( const seqan::BamAlignmentRecord &record_first, uint32_t end_pos ) const{
			// Proper forward reverse direction or read sized pair (proper direction is only checked if they are on same scaffold, so no check needed for that)
			return !hasFlagRC(record_first) && hasFlagNextRC(record_first) && maximum_insert_length_ >= end_pos - record_first.beginPos || record_first.beginPos == record_first.pNext;
		}

		bool ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const std::string &variant_file, uint64_t max_ref_seq_bin_size, uint16_t num_threads, bool calculate_bias = true ); // Fill the class with the information from a bam file

		bool Load( const char *archive_file );
		bool Save( const char *archive_file ) const;

		void PrepareProcessing();
		void PreparePlotting();
		void PrepareTesting();
	};

}

#endif // DATASTATS_H
