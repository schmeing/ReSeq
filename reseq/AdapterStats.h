#ifndef ADAPTERSTATS_H
#define ADAPTERSTATS_H

#include <atomic>
#include <stdint.h>
#include <string>
#include <vector>

#include <boost/serialization/vector.hpp>

#include <seqan/bam_io.h>

#include "Reference.h"
#include "utilities.h"
#include "Vect.hpp"

namespace reseq{
	class AdapterStats{
	private:
		// Definitions
		const uint16_t minimum_adapter_length_unmapped_; // Minimum adapter length required to insert not completely mapped read pair into statistics
		const double min_fraction_of_maximum_for_simulation_; // Minimum fraction of the highest adapter count that an adapter needs to have to be simulated (Removes garbage from adapter reading stage)

		// Temporary variables
		std::vector<std::vector<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>>> tmp_counts_; // counts_[AdapterID][secondAdapterID][firstAdapterLength][secondAdapterLength] = #adapters
		std::array<std::vector<std::vector<utilities::VectorAtomic<uint64_t>>>, 2> tmp_start_cut_; // start_cut_[templateSegment][AdapterID][basesCutFromTheBeginningOfTheAdapterAtPos0] = #adapters
		std::vector<utilities::VectorAtomic<uint64_t>> tmp_polya_tail_length_; // polya_tail_length_[lengthOfPolyATailAfterAdapter] = #adapters
		std::array<std::atomic<uint64_t>, 5> tmp_overrun_bases_; // overrun_bases_[nucleotide] = #basesAfterPolyATailWithThisNucleotide

		// Adapter infos loaded
		std::array<std::vector<std::string>, 2> names_; // names_[templateSegment][AdapterID] = adapterNameFromFasta
		std::array<std::vector<seqan::DnaString>, 2> seqs_; // seqs_[templateSegment][AdapterID] = adapterSequenceFromFasta
		std::vector<std::vector<bool>> combinations_; // combinations_[AdapterID][secondAdapterID] = notValid(false)/Valid(true)
		
		// Collected statistics
		std::vector<std::vector<Vect<Vect<uint64_t>>>> counts_; // counts_[AdapterID][secondAdapterID][firstAdapterLength][secondAdapterLength] = #adapters
		std::array<std::vector<Vect<uint64_t>>, 2> start_cut_; // start_cut_[templateSegment][AdapterID][basesCutFromTheBeginningOfTheAdapterAtPos0] = #adapters
		Vect<uint64_t> polya_tail_length_; // polya_tail_length_[lengthOfPolyATailAfterAdapter] = #adapters
		std::array<uint64_t, 5> overrun_bases_; // overrun_bases_[nucleotide] = #basesAfterPolyATailWithThisNucleotide

		// Calculated variables
		std::array<std::vector<uint64_t>, 2> count_sum_; // count_sum_[templateSegment][AdapterID] = #adapters
		std::array<std::vector<uint64_t>, 2> significant_count_; // significant_count_[templateSegment][AdapterID] = #adapters (All adapters not appearing often are set to zero)
		
		// Helper functions
		static uint32_t CountErrors(uint32_t &last_error_pos, const seqan::BamAlignmentRecord &record, const Reference &reference);
		bool VerifyInsert(const seqan::CharString &read1, const seqan::CharString &read2, int32_t pos1, int32_t pos2);
		static bool AdaptersAmbigous(const seqan::DnaString &adaper1, const seqan::DnaString &adaper2, uint64_t max_length);
		
		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & names_;
			ar & combinations_;
			
			ar & counts_;
			ar & start_cut_;
			ar & polya_tail_length_;
			ar & overrun_bases_;

			// Work around so that seqan>>DnaString can be stored as no serialization function exists for that class
			std::array<std::vector<std::string>, 2> seqs_archive;
			if(seqs_.at(0).size()){
				// Serialization writing
				for( uint16_t template_segment=2; template_segment--; ){
					seqs_archive.at(template_segment).reserve(seqs_.at(template_segment).size());
					for( auto &adapter : seqs_.at(template_segment) ){
						seqs_archive.at(template_segment).push_back( std::string(seqan::toCString(static_cast<seqan::CharString>(adapter))) );
					}
				}
			}

			ar & seqs_archive;

			if(!seqs_.at(0).size()){
				// Serialization loading
				for( uint16_t template_segment=2; template_segment--; ){
					seqs_.at(template_segment).clear();
					seqs_.at(template_segment).reserve(seqs_archive.at(template_segment).size());
					for( auto &adapter : seqs_archive.at(template_segment) ){
						seqs_.at(template_segment).push_back( seqan::DnaString(adapter.c_str()) );
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
		const std::string &Name(uint16_t template_segment, uint16_t id) const{ return names_.at(template_segment).at(id); }
		const seqan::DnaString &Sequence(uint16_t template_segment, uint16_t id) const{ return seqs_.at(template_segment).at(id); }
		const std::vector<uint64_t> &Counts(uint16_t template_segment) const{ return count_sum_.at(template_segment); }
		const std::vector<uint64_t> &SignificantCounts(uint16_t template_segment) const{ return significant_count_.at(template_segment); }
		const Vect<uint64_t> &StartCut(uint16_t template_segment, uint16_t id) const{ return start_cut_.at(template_segment).at(id); }
		const Vect<uint64_t> &PolyATailLength() const{ return polya_tail_length_; }
		const std::array<uint64_t, 5> &OverrunBases() const{ return overrun_bases_; }
	
		// Main functions
		bool LoadAdapters(const char *adapter_file, const char *adapter_matrix, uint16_t phred_quality_offset, uint32_t size_read_length);
		bool Detect(uint32_t &adapter_position_first, uint32_t &adapter_position_second, const seqan::BamAlignmentRecord &record_first, const seqan::BamAlignmentRecord &record_second, const Reference &reference, bool properly_mapped=false);
		void Finalize();
		void SumCounts();
		void Shrink();
		void PrepareSimulation();
	};
}

#endif // ADAPTERSTATS_H
