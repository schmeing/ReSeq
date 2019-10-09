#ifndef ERRORSTATS_H
#define ERRORSTATS_H

#include <array>
#include <stdint.h>
#include <vector>

#include <seqan/basic.h>

#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	class ErrorStats{
	public:
		enum InDelDef{
			kNoInDel,
			kDeletion,
			kInsertionA,
			kInsertionC,
			kInsertionG,
			kInsertionT,
			kInsertionN
		};
	private:
		// Temporary variables
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_called_bases_by_base_quality_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_called_bases_by_position_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_called_bases_by_error_num_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_called_bases_by_error_rate_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_error_num_by_quality_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_error_num_by_position_per_tile_;
		std::array<std::array<std::array<std::vector<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>>,5>,4>,2> tmp_error_num_by_error_rate_per_tile_;

		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_indel_by_indel_pos_;
		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_indel_by_position_;
		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_indel_by_gc_;
		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_indel_pos_by_position_;
		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_indel_pos_by_gc_;
		std::array<std::array<std::vector<std::vector<utilities::VectorAtomic<uintNucCount>>>,6>,2> tmp_gc_by_position_;

		std::array<std::vector<utilities::VectorAtomic<uintFragCount>>,2> tmp_errors_per_read_;
		std::array<std::array<std::array<std::array<std::vector<utilities::VectorAtomic<uintNucCount>>,6>,5>,4>,2> tmp_called_bases_by_base_quality_per_previous_called_base_;

		// Collected variables for estimation
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> called_bases_by_base_quality_per_tile_; // called_bases_by_base_quality_per_tile_[first/second][refBase][domError][tileId][calledBase][baseQuality] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> called_bases_by_position_per_tile_; // called_bases_by_position_per_tile_[first/second][refBase][domError][tileId][calledBase][position] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> called_bases_by_error_num_per_tile_; // called_bases_by_error_num_per_tile_[first/second][refBase][domError][tileId][calledBase][numErrors] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> called_bases_by_error_rate_per_tile_; // called_bases_by_error_rate_per_tile_[first/second][refBase][domError][tileId][calledBase][errorRate] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> error_num_by_quality_per_tile_; // error_num_by_quality_per_tile_[first/second][refBase][domError][tileId][numErrors][baseQuality] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> error_num_by_position_per_tile_; // error_num_by_position_per_tile_[first/second][refBase][domError][tileId][numErrors][position] = #bases
		std::array<std::array<std::array<Vect<Vect<Vect<uintNucCount>>>,5>,4>,2> error_num_by_error_rate_per_tile_; // error_num_by_error_rate_per_tile_[first/second][refBase][domError][tileId][numErrors][errorRate] = #bases

		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> indel_by_indel_pos_; // indel_by_indel_pos_[Insertion/DeletionBefore][PreviousRegularCall][NumInDelsDirectlyBefore][Nothing/Deletion/Insertion(A/C/G/T)]
		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> indel_by_position_; // indel_by_position_[Insertion/DeletionBefore][PreviousRegularCall][Position][Nothing/Deletion/Insertion(A/C/G/T)]
		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> indel_by_gc_; // indel_by_gc_[Insertion/DeletionBefore][PreviousRegularCall][GC][Nothing/Deletion/Insertion(A/C/G/T)]
		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> indel_pos_by_position_; // indel_pos_by_position_[Insertion/DeletionBefore][PreviousRegularCall][Position][NumInDelsDirectlyBefore]
		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> indel_pos_by_gc_; // indel_pos_by_gc_[Insertion/DeletionBefore][PreviousRegularCall][GC][NumInDelsDirectlyBefore]
		std::array<std::array<Vect<Vect<uintNucCount>>,6>,2> gc_by_position_; // position_by_gc_[Insertion/DeletionBefore][PreviousRegularCall][Position][GC]

		// Collected variables for plotting
		std::array<Vect<uintFragCount>,2> errors_per_read_; // errors_per_read_[first/second][#errors] = #reads
		std::array<std::array<std::array<std::array<Vect<uintNucCount>,6>,5>,4>,2> called_bases_by_base_quality_per_previous_called_base_; // called_bases_by_base_quality_per_previous_called_base_[first/second][refBase][calledBase][previouslyCalledBase][baseQuality] = #bases

		// Calculated variables for simulation from variables for estimation
		uintReadLen max_len_deletion_; //max_deletion_[first/second] = Maximum length of deletion

		// Calculated variables for plotting from variables for estimation
		std::array<std::array<std::array<Vect<uintNucCount>,5>,4>,2> called_bases_by_base_quality_; // called_bases_by_base_quality_[first/second][refBase][calledBase][baseQuality] = #bases
		std::array<std::array<std::array<Vect<uintNucCount>,5>,4>,2> called_bases_by_position_; // called_bases_by_position_[first/second][refBase][calledBase][position] = #bases
		
		std::array<Vect<uintNucCount>,2> indel_error_by_length_; // indel_error_by_length_[Insertion/Deletion][LengthOfInDel] = #indels
		std::array<Vect<uintNucCount>,2> indel_error_by_position_; // indel_error_by_position_[Insertion/Deletion][ReadPosition] = #indelBases
		std::array<Vect<uintNucCount>,2> indel_error_by_gc_; // indel_error_by_gc_[Insertion/Deletion][GC] = #indelBases

		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & called_bases_by_base_quality_per_tile_;
			ar & called_bases_by_position_per_tile_;
			ar & called_bases_by_error_num_per_tile_;
			ar & called_bases_by_error_rate_per_tile_;
			ar & error_num_by_quality_per_tile_;
			ar & error_num_by_position_per_tile_;
			ar & error_num_by_error_rate_per_tile_;

			ar & indel_by_indel_pos_;
			ar & indel_by_position_;
			ar & indel_by_gc_;
			ar & indel_pos_by_position_;
			ar & indel_pos_by_gc_;
			ar & gc_by_position_;

			ar & errors_per_read_;
			ar & called_bases_by_base_quality_per_previous_called_base_;
		}

		// Google test
		friend class ErrorStatsTest;
		friend class ProbabilityEstimatesTest;

	public:
		// Getter functions
		inline const Vect<Vect<uintNucCount>> &CalledBasesByBaseQuality(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &CalledBasesByPosition(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &CalledBasesByErrorNum(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &CalledBasesByErrorRate(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &ErrorNumByBaseQuality(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &ErrorNumByPosition(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &ErrorNumByErrorRate(uintTempSeq template_segment, uintTileId tile_id, seqan::Dna5 ref_base, seqan::Dna5 dom_error) const{
			return error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id);
		}

		inline const Vect<Vect<uintNucCount>> &InDelByInDelPos(uintInDelType indel_type, uintBaseCall last_call) const{
			return indel_by_indel_pos_.at(indel_type).at(last_call);
		}
		inline const Vect<Vect<uintNucCount>> &InDelByPosition(uintInDelType indel_type, uintBaseCall last_call) const{
			return indel_by_position_.at(indel_type).at(last_call);
		}
		inline const Vect<Vect<uintNucCount>> &InDelByGC(uintInDelType indel_type, uintBaseCall last_call) const{
			return indel_by_gc_.at(indel_type).at(last_call);
		}
		inline const Vect<Vect<uintNucCount>> &InDelPosByPosition(uintInDelType indel_type, uintBaseCall last_call) const{
			return indel_pos_by_position_.at(indel_type).at(last_call);
		}
		inline const Vect<Vect<uintNucCount>> &InDelPosByGC(uintInDelType indel_type, uintBaseCall last_call) const{
			return indel_pos_by_gc_.at(indel_type).at(last_call);
		}
		inline const Vect<Vect<uintNucCount>> &GCByPosition(uintInDelType indel_type, uintBaseCall last_call) const{
			return gc_by_position_.at(indel_type).at(last_call);
		}

		inline const Vect<uintFragCount> &ErrorsPerRead(uintTempSeq template_segment) const{
			return errors_per_read_.at(template_segment);
		}
		inline const Vect<uintNucCount> &CalledBasesByBaseQualityPerPreviousCalledBase(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base, uintBaseCall previous_base) const{
			return called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(called_base).at(previous_base);
		}

		inline uintReadLen MaxLenDeletion() const{
			return max_len_deletion_;
		}

		inline const Vect<uintNucCount> &CalledBasesByBaseQuality(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base) const{
			return called_bases_by_base_quality_.at(template_segment).at(ref_base).at(called_base);
		}
		inline const Vect<uintNucCount> &CalledBasesByPosition(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base) const{
			return called_bases_by_position_.at(template_segment).at(ref_base).at(called_base);
		}
		
		inline const Vect<uintNucCount> &InDelErrorByLength(uintInDelType indel_type) const{
			return indel_error_by_length_.at(indel_type);
		}
		inline const Vect<uintNucCount> &InDelErrorByPosition(uintInDelType indel_type) const{
			return indel_error_by_position_.at(indel_type);
		}
		inline const Vect<uintNucCount> &InDelErrorByGC(uintInDelType indel_type) const{
			return indel_error_by_gc_.at(indel_type);
		}

		// Setter functions
		inline void AddBase(
				uintTempSeq template_segment,
				uintBaseCall ref_base,
				uintBaseCall dom_error,
				uintTileId tile_id,
				uintBaseCall called_base,
				uintQual quality,
				uintReadLen read_pos,
				uintReadLen num_errors,
				uintPercent error_rate){
			// For estimation
			++tmp_called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(called_base).at(quality);
			++tmp_called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(called_base).at(read_pos);
			++tmp_called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(called_base).at(num_errors);
			++tmp_called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(called_base).at(error_rate);
			++tmp_error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(num_errors).at(quality);
			++tmp_error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(num_errors).at(read_pos);
			++tmp_error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).at(tile_id).at(num_errors).at(error_rate);
		}

		inline void AddInDel(
				uintInDelType indel_type,
				uintBaseCall last_call,
				InDelDef indel,
				uintReadLen indel_pos,
				uintReadLen pos,
				uintPercent gc){
			// For estimation
			++tmp_indel_by_indel_pos_.at(indel_type).at(last_call).at(indel_pos).at(indel);
			++tmp_indel_by_position_.at(indel_type).at(last_call).at(pos).at(indel);
			++tmp_indel_by_gc_.at(indel_type).at(last_call).at(gc).at(indel);
			++tmp_indel_pos_by_position_.at(indel_type).at(last_call).at(pos).at(indel_pos);
			++tmp_indel_pos_by_gc_.at(indel_type).at(last_call).at(gc).at(indel_pos);
			++tmp_gc_by_position_.at(indel_type).at(last_call).at(pos).at(gc);
		}

		inline void AddBasePlotting(
				uintTempSeq template_segment,
				uintBaseCall ref_base,
				uintBaseCall called_base,
				uintQual quality,
				uintBaseCall last_called_base){
			// For plotting
			++tmp_called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(called_base).at(last_called_base).at(quality);
		}

		inline void AddRead(uintTempSeq template_segment, uintReadLen num_errors){
			++tmp_errors_per_read_.at(template_segment).at(num_errors);
		}
	
		// Main functions
		void Prepare(uintTileId num_tiles, uintQual size_qual, uintReadLen size_pos, uintReadLen size_indel);
		bool Finalize();
		void Shrink();
		void PreparePlotting();
		void PrepareSimulation();
	};
}

#endif // ERRORSTATS_H
