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
#include "utilities.hpp"

namespace reseq{
	class Simulator{
	private:
		// Subclasses and structures
		class VariantBiasVarModifiers{
		public:
			// Keeping track of which variants to process
			intVariantId first_variant_id_;
			std::vector<intVariantId> unhandled_variant_id_;
			std::vector<uintSeqLen> unhandled_bases_in_variant_;

			// Handle indels at start of simulated fragment
			//const Reference::Variant *start_variant_;
			//std::vector<const Reference::Variant *> start_insertions_;
			uintSeqLen start_variant_pos_;

			// GC modification
			std::vector<intSeqShift> gc_mod_; // Change in gc for that allele compared to the reference for the given fragment
			std::vector<intSeqShift> end_pos_shift_;
			intSeqShift fragment_length_extension_; // Extension after the ReferenceSequence needed for Variants

			// Surroundings
			std::vector<std::array<intSurrounding, Reference::num_surrounding_blocks_>> mod_surrounding_start_, mod_surrounding_end_; // Contain the modified surroundings for the given allele if variation exists

			VariantBiasVarModifiers(intVariantId first_variant_id, uintAlleleId num_alleles):
					first_variant_id_(first_variant_id),
					start_variant_pos_(0),
					fragment_length_extension_(0){
				if(num_alleles){
					mod_surrounding_start_.resize(num_alleles);
					mod_surrounding_end_.resize(num_alleles);
				}
			}

			std::pair<intVariantId, uintSeqLen> StartVariant() const{
				return {first_variant_id_, start_variant_pos_};
			}
			std::pair<intVariantId, uintSeqLen> EndVariant(const std::vector<Reference::Variant> &variants, uintAlleleId allele, uintSeqLen cur_end_position) const{
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
			uintSeqLen position_;
			std::vector<std::pair<seqan::Dna5,uintPercent>> var_errors_;
			uintAlleleBitArray allele_;

			SysErrorVariant(uintSeqLen position, const std::vector<std::pair<seqan::Dna5,uintPercent>> &var_errors, uintAlleleBitArray allele):
				position_(position),
				var_errors_(var_errors),
				allele_(allele){
			}

			inline bool InAllele(uintAlleleId allele) const{
				return (allele_ >> allele) & 1; // Get bit corresponding to allele
			}
		};

		struct SimBlock{
			const uintRefSeqBin id_;
			const uintSeqLen start_pos_;
			std::atomic<bool> finished_; // In forward direction it stores that the simulation is finished and the block can be removed, in reverse direction it stores that a thread is already processing this block pair
			std::atomic<SimBlock *> next_block_;
			SimBlock *partner_block_; // forward direction stores here the corresponding reverse direction; reverse direction stores the previous block, which will be the partner for the next forward block
			uintSeed seed_;
			std::vector<std::pair<seqan::Dna5,uintPercent>> sys_errors_; //sys_errors_[pos] = {domError, errorRate}
			std::vector<SysErrorVariant> err_variants_;
			intVariantId first_variant_id_;

			SimBlock(uintRefSeqBin id, uintSeqLen start_pos, SimBlock *partner_block, uintSeed seed):
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
			const uintRefSeqId ref_seq_id_;
			SimBlock *first_block_;
			SimBlock *last_block_;
			SimUnit *next_unit_;

			SimUnit(uintRefSeqId ref_seq_id):
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
			std::discrete_distribution<uintTileId> tile_;

			// Adapter probabilities
			std::array< std::discrete_distribution<uintAdapterId>, 2> adapter_;
			std::discrete_distribution<uintReadLen> polya_tail_length_;
			std::discrete_distribution<uintBaseCall> overrun_bases_;

			std::array< std::discrete_distribution<uintReadLen>, 2> read_length_;

			// General probability
			std::uniform_real_distribution<double> zero_to_one_;

		public:
			GeneralRandomDistributions(const DataStats &stats):
				stats_(stats),
				tile_(stats.Tiles().Abundance().begin(), stats.Tiles().Abundance().end()),
				polya_tail_length_(stats.Adapters().PolyATailLength().begin(), stats.Adapters().PolyATailLength().end()),
				overrun_bases_(stats.Adapters().OverrunBases().begin(), stats.Adapters().OverrunBases().end()-1), // We don't want the N at the end
				zero_to_one_(0, 1){
				for(uintTempSeq template_segment=2; template_segment--; ){
					adapter_.at(template_segment) = std::discrete_distribution<uintAdapterId>(stats.Adapters().SignificantCounts(template_segment).begin(), stats.Adapters().SignificantCounts(template_segment).end());
					read_length_.at(template_segment) = std::discrete_distribution<uintReadLen>( stats.ReadLengths(template_segment).begin(), stats.ReadLengths(template_segment).end() );
				}
			}

			uintTileId TileId(std::mt19937_64 &rgen){
				if(1 < stats_.Tiles().Abundance().size() ){
					return tile_(rgen);
				}
				return 0;
			}
			uintSeqLen Adapter(uintTempSeq template_segment, std::mt19937_64 &rgen){ return adapter_.at(template_segment)(rgen); }
			uintReadLen PolyATailLength(std::mt19937_64 &rgen){ return polya_tail_length_(rgen) + stats_.Adapters().PolyATailLength().from(); }
			seqan::Dna OverrunBases(std::mt19937_64 &rgen){ return overrun_bases_(rgen); }
			uintReadLen ReadLength(uintTempSeq template_segment, uintSeqLen fragment_length, std::mt19937_64 &rgen){
				if( 1 == stats_.ReadLengths(template_segment).size() ){
					return stats_.ReadLengths(template_segment).from();
				}

				double random_value = ZeroToOne(rgen) * stats_.FragmentDistribution().InsertLengths().at(fragment_length);
				double counter = 0.0;
				uintReadLen read_len = stats_.ReadLengthsByFragmentLength(template_segment).at(fragment_length).to();
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

				for(uintTempSeq template_segment=2; template_segment--; ){
					adapter_.at(template_segment).reset();
					read_length_.at(template_segment).reset();
				}
			}
		};

		struct ReadFillParameter{
			uintReadLen read_length_;
			uintReadLen read_pos_;

			uintInDelType previous_indel_type_;
			uintReadLen indel_pos_;
			uintBaseCall base_call_;
			uintReadLen gc_seq_;

			uintQual seq_qual_;
			uintQual qual_;
			uintPercent error_rate_;

			uintReadLen num_errors_;

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
			std::array<utilities::CigarString, 2> cigar_;
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
		const uintFragCount kBatchSize = 100000; // 100,000 reads are written in each batch, balance between speed and memory usage
		const uintSeqLen kBlockSize = 1000; // Continuous bases that are handled by a single thread; too high: work load isn't shared properly, memory allocation of big arrays slows program down; too low: overhead like random generator seeding slows program down

		// Mutex
		std::mutex print_mutex_;
		std::mutex output_mutex_;
		std::mutex flush_mutex_;

		std::mutex block_creation_mutex_;
		std::mutex var_read_mutex_;

		// private variables
		std::array<seqan::SeqFileOut, 2> dest_;
		uintFragCount written_records_;
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
		std::array<uintReadLen,4> sys_seq_content_last5_;
		std::vector<uintSeqLen> sys_last_var_pos_per_allele_;
		std::vector<seqan::Dna5> sys_last_base_per_allele_;
		std::vector<seqan::Dna5> sys_dom_last5_per_allele_;
		std::vector<std::array<uintReadLen,4>> sys_seq_content_last5_per_allele_;
		std::vector<seqan::DnaString> sys_seq_last5_per_allele_;
		uintReadLen sys_gc_;
		uintReadLen sys_gc_bases_;
		uintReadLen sys_gc_range_;
		uintSeqLen distance_to_start_of_error_region_;
		std::array<std::vector<std::vector<std::pair<seqan::Dna5,uintPercent>>>, 2> adapter_sys_error_;

		double bias_normalization_;
		double non_zero_threshold_; // Threshold that random number from fragment generation has to reach or fragment count will be zero and does not have to be calculated
		uintFragCount total_pairs_;
		uintFragCount num_adapter_only_pairs_;
		std::atomic_flag adapter_only_simulated_; // Int instead of bool so atomic increment works

		std::array<seqan::StringSet<seqan::CharString> *, 2> output_ids_;
		std::array<seqan::StringSet<seqan::Dna5String> *, 2> output_seqs_;
		std::array<seqan::StringSet<seqan::CharString> *, 2> output_quals_;

		std::vector<double> tmp_probabilities_;
		std::uniform_real_distribution<double> rdist_zero_to_one_;

		// private functions
		static double CoveragePropLostFromAdapters(const DataStats &stats);
		static uintFragCount CoverageToNumberPairs(double coverage, uintRefLenCalc total_ref_size, double average_read_length, double adapter_part);
		static double NumberPairsToCoverage(uintFragCount total_pairs, uintRefLenCalc total_ref_size, double average_read_length, double adapter_part);

		bool Flush();
		bool Output(const SimPair &sim_reads);

		void ReadSystematicErrors(std::vector<std::pair<seqan::Dna5,uintPercent>> &sys_errors, uintSeqLen start_pos, uintSeqLen end_pos){
			uintPercent error_rate;
			for(auto pos=start_pos; pos < end_pos; ++pos){
				error_rate = sys_err_perc_[pos] - 33;
				if(86 < error_rate){
					error_rate += error_rate-86; // Decompress percentages again (odd percentages over 86 are removed, as fastq can store only up to 94 quality values)
				}
				sys_errors.emplace_back(sys_dom_err_[pos], error_rate);
			}
		}

		inline uintPercent DrawSystematicError(std::vector<std::pair<seqan::Dna5,uintPercent>> &sys_errors, seqan::Dna5 ref_base, seqan::Dna5 last_base, seqan::Dna5 dom_base, uintPercent gc_percent, uintSeqLen start_dist_error_region, const ProbabilityEstimates &estimates){
			double prob_sum;

			seqan::Dna5 dom_error = estimates.DominantError(ref_base, last_base, dom_base).Draw(tmp_probabilities_, prob_sum, {start_dist_error_region, gc_percent}, rdist_zero_to_one_(block_seed_gen_));
			if(0.0 == prob_sum){
				dom_error = 4;
			}

			uintPercent error_rate = estimates.ErrorRate(ref_base, dom_error).Draw(tmp_probabilities_, prob_sum, {start_dist_error_region, gc_percent}, rdist_zero_to_one_(block_seed_gen_));
			if(0.0 == prob_sum){
				error_rate = 0;
			}
			sys_errors.emplace_back(dom_error, error_rate);

			return error_rate;
		}

		template<typename T> inline void UpdateGC( uintReadLen &gc, uintReadLen &bases, uintReadLen pos, const T &sequence ){
			if( utilities::IsGC(sequence, pos) ){
				++gc;
			}
			if(bases < sys_gc_range_){
				++bases;
			}
			else{
				if( utilities::IsGC(sequence, pos - bases) ){
					--gc;
				}
			}
		}

		template<typename T> void SetSystematicErrors(std::vector<std::pair<seqan::Dna5,uintPercent>> &sys_errors, const T &sequence, uintSeqLen start_pos, uintSeqLen end_pos, const DataStats &stats, const ProbabilityEstimates &estimates){
			seqan::Dna5 ref_base;
			uintPercent error_rate;

			for(auto pos=start_pos; pos < end_pos; ++pos){
				ref_base = sequence[pos];
				error_rate = DrawSystematicError(sys_errors, ref_base, sys_last_base_, sys_dom_last5_, utilities::SafePercent(sys_gc_, sys_gc_bases_), utilities::TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_), estimates);

				// Set values for next iteration
				sys_last_base_ = ref_base;
				utilities::SetDominantLastX( sys_dom_last5_, sys_seq_content_last5_, ref_base, 5, sequence, pos);
				stats.Coverage().UpdateDistances(distance_to_start_of_error_region_, error_rate);
				UpdateGC(sys_gc_, sys_gc_bases_, pos, sequence);
			}
		}

		inline void IncrementBlockPos(uintSeqLen &block_pos, const SimBlock* &block, intVariantId &cur_var);
		void GetSysErrorFromBlock(seqan::Dna5 &dom_error, uintPercent &error_rate, const SimBlock* &block, uintSeqLen &block_pos, intVariantId &cur_var, uintSeqLen &var_pos, uintAlleleId allele);
		bool FillReadPart( SimPair &sim_reads, uintTempSeq template_segment, uintTileId tile_id, const seqan::DnaString &org_sequence, uintReadLen org_pos, char base_cigar_element, const SimBlock* block, uintSeqLen block_pos, uintAdapterId adapter_id, ReadFillParameter &par, std::pair<intVariantId, uintSeqLen> first_variant, uintAlleleId allele, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen);
		bool FillRead( SimPair &sim_reads, uintTempSeq template_segment, uintTileId tile_id, uintSeqLen fragment_length, const SimBlock* start_block, uintSeqLen block_start_pos, std::pair<intVariantId, uintSeqLen> first_variant, uintAlleleId allele, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen);

		void CreateReadId( seqan::CharString &id, uintRefSeqBin block_number, uintFragCount read_number, uintSeqLen start_pos, uintSeqLen end_pos, uintAlleleId allele, uintTile tile, uintRefSeqId ref_seq_id, const Reference &ref, const utilities::CigarString &cigar);
		bool CreateReads( const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, SimPair &sim_reads, std::mt19937_64 &rgen, uintFragCount counts, bool strand, uintAlleleId allele, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintFragCount &read_number, const SimBlock *start_block, uintSeqLen start_position_forward, uintSeqLen end_position_forward, std::pair<intVariantId, uintSeqLen> start_variant={0,0}, std::pair<intVariantId, uintSeqLen> end_variant={0,0} );

		void ResetSystematicErrorCounters(const Reference &ref);
		bool LoadSysErrorRecord(uintRefSeqId ref_id, const Reference &ref);
		void SetSystematicErrorVariantsReverse(uintSeqLen &start_dist_error_region, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates);
		bool CreateUnit(uintRefSeqId ref_id, uintRefSeqBin first_block_id, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&first_reverse_block, SimUnit *&unit);
		void SetSystematicErrorVariantsForward(uintSeqLen &start_dist_error_region, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates);
		bool CreateBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, Benchmark &time );
		bool GetNextBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&block, SimUnit *&unit, Benchmark &time );

		inline bool AlleleSkipped(const VariantBiasVarModifiers &bias_mod, uintAlleleId allele, const std::vector<Reference::Variant> &variants, uintSeqLen cur_start_position) const{
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
		inline bool VariantInsideCurrentFragment( intVariantId cur_var_id, uintSeqLen cur_end_position, intSeqShift end_pos_shift, uintRefSeqId ref_seq_id, const Reference &ref ) const;
		void HandleGCModAndEndPosShiftForNewVariants( VariantBiasVarModifiers &bias_mod, uintAlleleId allele, uintRefSeqId ref_seq_id, const Reference &ref, uintSeqLen cur_end_position ) const;
		static void ChangeSurroundingBase( std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding, uintSurPos pos, seqan::Dna new_base );
		static void DeleteSurroundingBaseShiftingOnRightSide( std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding, uintSurPos pos, seqan::Dna new_end_base );
		static void DeleteSurroundingBaseShiftingOnLeftSide( std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding, uintSurPos pos, seqan::Dna new_end_base );
		static void InsertSurroundingBasesShiftingOnRightSide( std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding, uintSurPos pos, seqan::DnaString new_bases );
		static void InsertSurroundingBasesShiftingOnLeftSide( std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding, uintSurPos pos, seqan::DnaString new_bases );
		void HandleSurroundingVariantsBeforeCenter(std::array<intSurrounding, Reference::num_surrounding_blocks_> &mod_surrounding, uintSeqLen center_position, intSeqShift initial_pos_shift, intVariantId center_var, uintAlleleId allele, uintRefSeqId ref_seq_id, const Reference &ref, bool reverse) const;
		void HandleSurroundingVariantsAfterCenter(std::array<intSurrounding, Reference::num_surrounding_blocks_> &mod_surrounding, uintSeqLen center_position, intSeqShift initial_pos_shift, intVariantId center_var, uintAlleleId allele, uintRefSeqId ref_seq_id, const Reference &ref, bool reverse) const;
		void VariantModStartSurrounding(VariantBiasVarModifiers &bias_mod, uintAlleleId allele, uintRefSeqId ref_seq_id, const Reference &ref, uintSeqLen cur_start_position, const std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_start) const;
		void PrepareBiasModForCurrentStartPos( VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, const Reference &ref, uintSeqLen cur_start_position, uintSeqLen cur_end_position, std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_start ) const;
		void VariantModEndSurrounding(VariantBiasVarModifiers &bias_mod, uintAlleleId allele, uintRefSeqId ref_seq_id, const Reference &ref, uintSeqLen cur_end_position, const std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_end) const;
		void UpdateBiasModForCurrentFragmentLength( VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, const Reference &ref, uintSeqLen cur_start_position, uintSeqLen cur_end_position, std::array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_end ) const;
		void CheckForInsertedBasesToStartFrom( VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position, const Reference &ref ) const;
		void CheckForFragmentLengthExtension( VariantBiasVarModifiers &bias_mod, uintSeqLen fragment_length, uintSeqLen fragment_length_to, const Reference &ref, const DataStats &stats ) const;
		void GetOrgSeq(SimPair &sim_reads, bool strand, uintAlleleId allele, uintSeqLen fragment_length, uintSeqLen cur_start_position, uintSeqLen cur_end_position, uintRefSeqId ref_seq_id, const Reference &ref, const DataStats &stats, const VariantBiasVarModifiers &bias_mod) const;
		bool SimulateFromGivenBlock( const SimBlock &block, const SimUnit &unit, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen, Benchmark &time );
		bool SimulateAdapterOnlyPairs( const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, GeneralRandomDistributions &rdist, std::mt19937_64 &rgen );
		static void SimulationThread( Simulator &self, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates );

		bool WriteOutSystematicErrorProfile(const std::string &id, std::vector<std::pair<seqan::Dna5,uintPercent>> sys_errors, seqan::Dna5String dom_err, seqan::CharString err_perc);

		// Google test
		friend class SimulatorTest;

	public:
		Simulator();

		bool CreateSystematicErrorProfile(const char *destination_file, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, uintSeed seed);
		void Simulate(const char *destination_file_first, const char *destination_file_second, Reference &ref, DataStats &stats, const ProbabilityEstimates &estimates, uintNumThreads num_threads, uintSeed seed, uintFragCount num_read_pairs=0, double coverage=0.0, RefSeqBiasSimulation ref_bias_model=kKeep, const std::string &ref_bias_file=std::string(), const std::string &sys_error_file=std::string(), const std::string &record_base_identifier=std::string(), const std::string &var_file=std::string());
	};
}
#endif // SIMULATOR_H
