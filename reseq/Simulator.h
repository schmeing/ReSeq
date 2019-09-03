#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <atomic>
#include <mutex>
#include <random>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "SeqQualityStats.hpp"
#include "Vect.hpp"
#include "DataStats.h"
#include "FragmentDistributionStats.h"
#include "ProbabilityEstimates.h"
#include "Reference.h"

namespace reseq{
	class Simulator{
	private:
		// Subclasses and structures
		class VariantBiasVarModifiers{
		public:
			// Keeping track of which variants to process
			int32_t first_variant_id_;
			std::vector<int32_t> unhandled_variant_id_;
			std::vector<uint16_t> unhandled_bases_in_variant_;

			// Handle indels at start of simulated fragment
			//const Reference::Variant *start_variant_;
			//std::vector<const Reference::Variant *> start_insertions_;
			uint16_t start_variant_pos_;

			// GC modification
			std::vector<int32_t> gc_mod_; // Change in gc for that allele compared to the reference for the given fragment
			std::vector<int32_t> end_pos_shift_;
			int32_t fragment_length_extension_; // Extension after the ReferenceSequence needed for Variants

			// Surroundings
			std::vector<std::array<uint32_t, Reference::num_surrounding_blocks_>> mod_surrounding_start_, mod_surrounding_end_; // Contain the modified surroundings for the given allele if variation exists

			VariantBiasVarModifiers(uint32_t first_variant_id, uint16_t num_alleles):
					first_variant_id_(first_variant_id),
					start_variant_pos_(0),
					fragment_length_extension_(0){
				if(num_alleles){
					mod_surrounding_start_.resize(num_alleles);
					mod_surrounding_end_.resize(num_alleles);
				}
			}

			std::pair<int32_t, uint16_t> StartVariant() const{
				return {first_variant_id_, start_variant_pos_};
			}
			std::pair<int32_t, uint16_t> EndVariant(const std::vector<Reference::Variant> &variants, uint16_t allele, uint32_t cur_end_position) const{
				if(unhandled_bases_in_variant_.at(allele)){
					return {unhandled_variant_id_.at(allele), length(variants.at(unhandled_variant_id_.at(allele)).var_seq_)-unhandled_bases_in_variant_.at(allele)};
				}
				else if( start_variant_pos_ && variants.at(first_variant_id_).position_ == cur_end_position-1 ){
					return {first_variant_id_, start_variant_pos_-end_pos_shift_.at(allele)+1};
				}
				else{
					auto first_rev_variant = unhandled_variant_id_.at(allele);
					if(first_rev_variant == variants.size()){
						--first_rev_variant;
					}
					while(0 <= first_rev_variant && variants.at(first_rev_variant).position_ >= cur_end_position+end_pos_shift_.at(allele)+fragment_length_extension_){
						--first_rev_variant;
					}
					return {first_rev_variant, 0};
				}
			}
		};

		class SysErrorVariant{
		public:
			uint32_t position_;
			std::vector<std::pair<seqan::Dna5,uint8_t>> var_errors_;
			uint64_t allele_;

			SysErrorVariant(uint32_t position, const std::vector<std::pair<seqan::Dna5,uint8_t>> &var_errors, uint64_t allele):
				position_(position),
				var_errors_(var_errors),
				allele_(allele){
			}

			inline bool InAllele(uint16_t allele) const{
				return (allele_ >> allele) & 1; // Get bit corresponding to allele
			}
		};

		struct SimBlock{
			const uint64_t id_;
			const uint64_t start_pos_;
			std::atomic<bool> finished_; // In forward direction it stores that the simulation is finished and the block can be removed, in reverse direction it stores that a thread is already processing this block pair
			std::atomic<SimBlock *> next_block_;
			SimBlock *partner_block_; // forward direction stores here the corresponding reverse direction; reverse direction stores the previous block, which will be the partner for the next forward block
			uint64_t seed_;
			std::vector<std::pair<seqan::Dna5,uint8_t>> sys_errors_; //sys_errors_[pos] = {domError, errorRate}
			std::vector<SysErrorVariant> err_variants_;
			int32_t first_variant_id_;

			SimBlock(uint64_t id, uint64_t start_pos, SimBlock *partner_block, uint64_t seed):
				id_(id),
				start_pos_(start_pos),
				finished_(false),
				next_block_(NULL),
				partner_block_(partner_block),
				seed_(seed),
				first_variant_id_(0){
			}
		};

		struct SimUnit{
			const seqan::Size< seqan::StringSet<seqan::CharString> >::Type ref_seq_id_;
			SimBlock *first_block_;
			SimBlock *last_block_;
			SimUnit *next_unit_;

			SimUnit(seqan::Size< seqan::StringSet<seqan::CharString> >::Type ref_seq_id):
				ref_seq_id_(ref_seq_id),
				first_block_(NULL),
				last_block_(NULL),
				next_unit_(NULL){
			}
		};

		class GeneralRandomDistributions{
		private:
			const DataStats &stats_;

			// Tile probability
			std::discrete_distribution<uint64_t> tile_;

			// Adapter probabilities
			std::array< std::discrete_distribution<uint64_t>, 2> adapter_;
			std::discrete_distribution<uint64_t> polya_tail_length_;
			std::discrete_distribution<uint64_t> overrun_bases_;

			std::array< std::discrete_distribution<uint64_t>, 2> read_length_;

			// General probability
			std::uniform_real_distribution<double> zero_to_one_;

		public:
			GeneralRandomDistributions(const DataStats &stats):
				stats_(stats),
				tile_(stats.Tiles().Abundance().begin(), stats.Tiles().Abundance().end()),
				polya_tail_length_(stats.Adapters().PolyATailLength().begin(), stats.Adapters().PolyATailLength().end()),
				overrun_bases_(stats.Adapters().OverrunBases().begin(), stats.Adapters().OverrunBases().end()-1), // We don't want the N at the end
				zero_to_one_(0, 1){
				for(uint16_t template_segment=2; template_segment--; ){
					adapter_.at(template_segment) = std::discrete_distribution<uint64_t>(stats.Adapters().SignificantCounts(template_segment).begin(), stats.Adapters().SignificantCounts(template_segment).end());
					read_length_.at(template_segment) = std::discrete_distribution<uint64_t>( stats.ReadLengths(template_segment).begin(), stats.ReadLengths(template_segment).end() );
				}
			}

			uint16_t TileId(std::mt19937_64 &rgen){
				if(1 < stats_.Tiles().Abundance().size() ){
					return tile_(rgen);
				}
				return 0;
			}
			uint32_t Adapter(uint16_t template_segment, std::mt19937_64 &rgen){ return adapter_.at(template_segment)(rgen); }
			uint32_t PolyATailLength(std::mt19937_64 &rgen){ return polya_tail_length_(rgen) + stats_.Adapters().PolyATailLength().from(); }
			uint16_t OverrunBases(std::mt19937_64 &rgen){ return overrun_bases_(rgen); }
			uint32_t ReadLength(uint16_t template_segment, uint32_t fragment_length, std::mt19937_64 &rgen){
				if( 1 == stats_.ReadLengths(template_segment).size() ){
					return stats_.ReadLengths(template_segment).from();
				}

				double random_value = ZeroToOne(rgen) * stats_.FragmentDistribution().InsertLengths().at(fragment_length);
				double counter = 0.0;
				uint32_t read_len = stats_.ReadLengthsByFragmentLength(template_segment).at(fragment_length).to();
				while(counter <= random_value && (read_len-- > stats_.ReadLengthsByFragmentLength(template_segment).at(fragment_length).from()) ){
					counter += stats_.ReadLengthsByFragmentLength(template_segment).at(fragment_length).at(read_len);
				}

				return read_len;
			}
			double ZeroToOne(std::mt19937_64 &rgen){ return zero_to_one_(rgen); }

			void Reset(){
				// Shouldn't do anything for gcc implementations of discrete_distribution, but better be save than sorry
				tile_.reset();
				polya_tail_length_.reset();
				overrun_bases_.reset();
				zero_to_one_.reset();

				for(uint16_t template_segment=2; template_segment--; ){
					adapter_.at(template_segment).reset();
					read_length_.at(template_segment).reset();
				}
			}
		};

		struct ReadFillParameter{
			uint32_t read_length_;
			uint32_t read_pos_;

			uint16_t previous_indel_type_;
			uint16_t indel_pos_;
			uint8_t base_call_;
			uint32_t gc_seq_;

			uint8_t seq_qual_;
			uint8_t qual_;
			uint8_t error_rate_;

			uint16_t num_errors_;

			ReadFillParameter():
				read_pos_(0),
				previous_indel_type_(0),
				indel_pos_(0),
				base_call_(5),
				gc_seq_(0),
				qual_(1),
				error_rate_(0),
				num_errors_(0){
			}
		};

		struct SimPair{
			std::array<seqan::CharString, 2> id_;
			std::array<seqan::DnaString, 2> org_seq_;
			std::array<seqan::Dna5String, 2> seq_;
			std::array<seqan::CharString, 2> qual_;
			std::array<seqan::String<seqan::CigarElement<>>, 2> cigar_;
		};

		// Benchmarking
		struct Benchmark{
			double total_;
			double create_block_;
			double handle_block_;
			double get_counts_;
			double create_reads_;

			Benchmark():
				total_(0),
				create_block_(0),
				handle_block_(0),
				get_counts_(0),
				create_reads_(0){
			}
		};
		Benchmark time_;
		std::mutex time_mutex_;

		// Definitions
		const uint64_t batch_size_;
		const uint64_t block_size_;

		// Mutex
		std::mutex print_mutex_;
		std::mutex output_mutex_;
		std::mutex flush_mutex_;

		std::mutex block_creation_mutex_;
		std::mutex var_read_mutex_;

		// private variables
		std::array<seqan::SeqFileOut, 2> dest_;
		uint64_t written_records_;
		std::string record_base_identifier_;

		SimUnit *first_unit_;
		SimUnit *last_unit_;
		SimUnit *current_unit_;
		SimBlock *current_block_;
		std::mt19937_64 block_seed_gen_;
		std::atomic<bool> simulation_error_;

		bool sys_from_file_;
		seqan::SeqFileIn sys_file_;
		seqan::CharString sys_id_;
		seqan::Dna5String sys_dom_err_;
		seqan::CharString sys_err_perc_;

		seqan::Dna5 sys_last_base_;
		seqan::Dna5 sys_dom_last5_;
		std::array<uint64_t,4> sys_seq_content_last5_;
		std::vector<uint32_t> sys_last_var_pos_per_allele_;
		std::vector<seqan::Dna5> sys_last_base_per_allele_;
		std::vector<seqan::Dna5> sys_dom_last5_per_allele_;
		std::vector<std::array<uint64_t,4>> sys_seq_content_last5_per_allele_;
		std::vector<seqan::DnaString> sys_seq_last5_per_allele_;
		uint32_t sys_gc_;
		uint32_t sys_gc_bases_;
		uint32_t sys_gc_range_;
		uint32_t distance_to_start_of_error_region_;
		std::array<std::vector<std::vector<std::pair<seqan::Dna5,uint8_t>>>, 2> adapter_sys_error_;

		double bias_normalization_;
		double non_zero_threshold_; // Threshold that random number from fragment generation has to reach or fragment count will be zero and does not have to be calculated
		uint64_t total_pairs_;
		uint64_t num_adapter_only_pairs_;
		std::atomic_flag adapter_only_simulated_; // Int instead of bool so atomic increment works

		std::array<seqan::StringSet<seqan::CharString> *, 2> output_ids_;
		std::array<seqan::StringSet<seqan::Dna5String> *, 2> output_seqs_;
		std::array<seqan::StringSet<seqan::CharString> *, 2> output_quals_;

		std::vector<double> tmp_probabilities_;
		std::uniform_real_distribution<double> rdist_zero_to_one_;

		// private functions
		static double CoveragePropLostFromAdapters(const DataStats &stats);
		static uint64_t CoverageToNumberPairs(double coverage, uint64_t total_ref_size, double average_read_length, double adapter_part);
		static double NumberPairsToCoverage(uint64_t total_pairs, uint64_t total_ref_size, double average_read_length, double adapter_part);

		bool Flush();
		bool Output(const SimPair &sim_reads);

		void ReadSystematicErrors(std::vector<std::pair<seqan::Dna5,uint8_t>> &sys_errors, uint32_t start_pos, uint32_t end_pos){
			uint8_t error_rate;
			for(auto pos=start_pos; pos < end_pos; ++pos){
				error_rate = sys_err_perc_[pos] - 33;
				if(86 < error_rate){
					error_rate += error_rate-86; // Decompress percentages again (odd percentages over 85 are removed, as fastq can store only up to 94 quality values)
				}
				sys_errors.emplace_back(sys_dom_err_[pos], error_rate);
			}
		}

		static inline uint8_t GcPercent(uint32_t gc, uint32_t bases){
			if(bases){
				return reseq::utilities::Percent(gc,bases);
			}
			else{
				return 50;
			}
		}

		static inline uint32_t TransformDistanceToStartOfErrorRegion(uint32_t start_dist_error_region){
			return (start_dist_error_region+9)/10;
		}

		inline uint8_t DrawSystematicError(std::vector<std::pair<seqan::Dna5,uint8_t>> &sys_errors, seqan::Dna5 ref_base, seqan::Dna5 last_base, seqan::Dna5 dom_base, uint8_t gc_percent, uint32_t start_dist_error_region, const ProbabilityEstimates &estimates){
			double prob_sum;

			seqan::Dna5 dom_error = estimates.DominantError(ref_base, last_base, dom_base).Draw(tmp_probabilities_, prob_sum, {start_dist_error_region, gc_percent}, rdist_zero_to_one_(block_seed_gen_));
			if(0.0 == prob_sum){
				dom_error = 4;
			}

			uint8_t error_rate = estimates.ErrorRate(ref_base, dom_error).Draw(tmp_probabilities_, prob_sum, {start_dist_error_region, gc_percent}, rdist_zero_to_one_(block_seed_gen_));
			if(0.0 == prob_sum){
				error_rate = 0;
			}
			sys_errors.emplace_back(dom_error, error_rate);

			return error_rate;
		}

		template<typename T> inline void UpdateGC( uint32_t &gc, uint32_t &bases, uint32_t pos, const T &sequence ){
			if( 'C' == sequence[pos] || 'G' == sequence[pos] ){
				++gc;
			}
			if(bases < sys_gc_range_){
				++bases;
			}
			else{
				if( 'C' == sequence[pos - bases] || 'G' == sequence[pos - bases] ){
					--gc;
				}
			}
		}

		template<typename T> void SetSystematicErrors(std::vector<std::pair<seqan::Dna5,uint8_t>> &sys_errors, const T &sequence, uint32_t start_pos, uint32_t end_pos, const DataStats &stats, const ProbabilityEstimates &estimates){
			seqan::Dna5 ref_base;
			uint8_t error_rate;

			for(auto pos=start_pos; pos < end_pos; ++pos){
				ref_base = sequence[pos];
				error_rate = DrawSystematicError(sys_errors, ref_base, sys_last_base_, sys_dom_last5_, GcPercent(sys_gc_, sys_gc_bases_), TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_), estimates);

				// Set values for next iteration
				sys_last_base_ = ref_base;
				utilities::SetDominantLastX( sys_dom_last5_, sys_seq_content_last5_, ref_base, 5, sequence, pos);
				stats.Coverage().UpdateDistances(distance_to_start_of_error_region_, error_rate);
				UpdateGC(sys_gc_, sys_gc_bases_, pos, sequence);
			}
		}

		inline void IncrementBlockPos(uint32_t &block_pos, const SimBlock* &block, int32_t &cur_var);
		void GetSysErrorFromBlock(seqan::Dna5 &dom_error, uint8_t &error_rate, const SimBlock* &block, uint32_t &block_pos, int32_t &cur_var, uint16_t &var_pos, uint16_t allele);
		bool FillReadPart( SimPair &sim_reads, uint16_t template_segment, uint16_t tile_id, const seqan::DnaString &org_sequence, uint32_t org_pos, char base_cigar_element, const SimBlock* block, uint32_t block_pos, uint32_t adapter_id, ReadFillParameter &par, std::pair<int32_t, uint16_t> first_variant, uint16_t allele, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen);
		bool FillRead( SimPair &sim_reads, uint16_t template_segment, uint16_t tile_id, uint32_t fragment_length, const SimBlock* start_block, uint32_t block_start_pos, std::pair<int32_t, uint16_t> first_variant, uint16_t allele, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen);

		void CreateReadId( seqan::CharString &id, uint64_t block_number, uint64_t read_number, uint32_t start_pos, uint32_t end_pos, uint16_t allele, uint32_t tile, uint32_t ref_seq_id, const Reference &ref, const seqan::String<seqan::CigarElement<>> &cigar);
		bool CreateReads( const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, SimPair &sim_reads, std::mt19937_64 &rgen, uint64_t counts, bool strand, uint16_t allele, uint32_t ref_seq_id, uint32_t fragment_length, uint64_t &read_number, const SimBlock *start_block, uint32_t start_position_forward, uint32_t end_position_forward, std::pair<int32_t, uint16_t> start_variant={0,0}, std::pair<int32_t, uint16_t> end_variant={0,0} );

		void ResetSystematicErrorCounters(const Reference &ref);
		bool LoadSysErrorRecord(uint32_t ref_id, const Reference &ref);
		void SetSystematicErrorVariantsReverse(uint32_t &start_dist_error_region, SimBlock &block, uint32_t ref_seq_id, uint32_t end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates);
		bool CreateUnit(seqan::Size< seqan::StringSet<seqan::CharString> >::Type ref_id, uint64_t first_block_id, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&first_reverse_block, SimUnit *&unit);
		void SetSystematicErrorVariantsForward(uint32_t &start_dist_error_region, SimBlock &block, uint32_t ref_seq_id, uint32_t end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates);
		bool CreateBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, Benchmark &time );
		bool GetNextBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&block, SimUnit *&unit, Benchmark &time );

		inline bool AlleleSkipped(const VariantBiasVarModifiers &bias_mod, uint16_t allele, const std::vector<Reference::Variant> &variants, uint32_t cur_start_position) const{
			// Exclude alleles with a deletion at start position and alleles that do not belong to the current insertion run
			if( bias_mod.first_variant_id_ < variants.size() && variants.at(bias_mod.first_variant_id_).position_ == cur_start_position ){
				if( 0 == length( variants.at(bias_mod.first_variant_id_).var_seq_) ){
					return variants.at(bias_mod.first_variant_id_).InAllele(allele);
				}
				else if( bias_mod.start_variant_pos_ ){
					return !variants.at(bias_mod.first_variant_id_).InAllele(allele);
				}
			}

			return false;
		}
		inline bool VariantInsideCurrentFragment( uint32_t cur_var_id, uint32_t cur_end_position, int32_t end_pos_shift, uint32_t ref_seq_id, const Reference &ref ) const;
		void HandleGCModAndEndPosShiftForNewVariants( VariantBiasVarModifiers &bias_mod, uint16_t allele, uint32_t ref_seq_id, const Reference &ref, uint32_t cur_end_position ) const;
		static void ChangeSurroundingBase( std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding, uint16_t pos, seqan::Dna new_base );
		static void DeleteSurroundingBaseShiftingOnRightSide( std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding, uint16_t pos, seqan::Dna new_end_base );
		static void DeleteSurroundingBaseShiftingOnLeftSide( std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding, uint16_t pos, seqan::Dna new_end_base );
		static void InsertSurroundingBasesShiftingOnRightSide( std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding, uint16_t pos, seqan::DnaString new_bases );
		static void InsertSurroundingBasesShiftingOnLeftSide( std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding, uint16_t pos, seqan::DnaString new_bases );
		void HandleSurroundingVariantsBeforeCenter(std::array<uint32_t, Reference::num_surrounding_blocks_> &mod_surrounding, uint32_t center_position, int32_t initial_pos_shift, int32_t center_var, uint16_t allele, uint32_t ref_seq_id, const Reference &ref, bool reverse) const;
		void HandleSurroundingVariantsAfterCenter(std::array<uint32_t, Reference::num_surrounding_blocks_> &mod_surrounding, uint32_t center_position, int32_t initial_pos_shift, int32_t center_var, uint16_t allele, uint32_t ref_seq_id, const Reference &ref, bool reverse) const;
		void VariantModStartSurrounding(VariantBiasVarModifiers &bias_mod, uint16_t allele, uint32_t ref_seq_id, const Reference &ref, uint32_t cur_start_position, const std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding_start) const;
		void PrepareBiasModForCurrentStartPos( VariantBiasVarModifiers &bias_mod, uint32_t ref_seq_id, const Reference &ref, uint32_t cur_start_position, uint32_t cur_end_position, std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding_start ) const;
		void VariantModEndSurrounding(VariantBiasVarModifiers &bias_mod, uint16_t allele, uint32_t ref_seq_id, const Reference &ref, uint32_t cur_end_position, const std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding_end) const;
		void UpdateBiasModForCurrentFragmentLength( VariantBiasVarModifiers &bias_mod, uint32_t ref_seq_id, const Reference &ref, uint32_t cur_start_position, uint32_t cur_end_position, std::array<uint32_t, Reference::num_surrounding_blocks_> &surrounding_end ) const;
		void CheckForInsertedBasesToStartFrom( VariantBiasVarModifiers &bias_mod, uint32_t ref_seq_id, uint32_t cur_start_position, const Reference &ref ) const;
		void CheckForFragmentLengthExtension( VariantBiasVarModifiers &bias_mod, uint32_t fragment_length, uint32_t fragment_length_to, const Reference &ref, const DataStats &stats ) const;
		void GetOrgSeq(SimPair &sim_reads, bool strand, uint16_t allele, uint32_t fragment_length, uint32_t cur_start_position, uint32_t cur_end_position, uint32_t ref_seq_id, const Reference &ref, const DataStats &stats, const VariantBiasVarModifiers &bias_mod) const;
		bool SimulateFromGivenBlock( const SimBlock &block, const SimUnit &unit, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen, Benchmark &time );
		bool SimulateAdapterOnlyPairs( const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen );
		static void SimulationThread( Simulator &self, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates );

		bool WriteOutSystematicErrorProfile(const std::string &id, std::vector<std::pair<seqan::Dna5,uint8_t>> sys_errors, seqan::Dna5String dom_err, seqan::CharString err_perc);

		// Google test
		friend class SimulatorTest;

	public:
		Simulator();

		bool CreateSystematicErrorProfile(const char *destination_file, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, uint64_t seed);
		void Simulate(const char *destination_file_first, const char *destination_file_second, Reference &ref, DataStats &stats, const ProbabilityEstimates &estimates, uint16_t num_threads, uint64_t seed, uint64_t num_read_pairs=0, double coverage=0.0, RefSeqBiasSimulation ref_bias_model=kKeep, const std::string &ref_bias_file=std::string(), const std::string &sys_error_file=std::string(), const std::string &record_base_identifier=std::string(), const std::string &var_file=std::string());
	};
}
#endif // SIMULATOR_H
