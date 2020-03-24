#ifndef REFERENCE_H
#define REFERENCE_H

#include <atomic>
#include <array>
#include <fstream>
#include <set>
#include <stdint.h>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

#include "Surrounding.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	struct FragmentSite;

	class Reference{
	public:
#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the nested classes
		class Variant{
		public:
			static const uintAlleleId kMaxAlleles = 64; // One bit per allele is needed and uint64_t is used for allele_

			uintSeqLen position_;
			seqan::DnaString var_seq_;
			uintAlleleBitArray allele_;

			template<typename T> Variant(uintSeqLen position, const T &var_seq, uintAlleleBitArray allele):
				position_(position),
				var_seq_(var_seq),
				allele_(allele){
			}

			inline bool InAllele(uintAlleleId allele) const{
				return (allele_ >> allele) & 1; // Get bit corresponding to allele
			}

			inline uintAlleleId FirstAllele() const{
				if(allele_){
					uintAlleleId first = 0;
					while(!InAllele(first)){
						++first;
					}
					return first;
				}

				printErr << "Variant belonging to no allele at position " << position_ << ": " << var_seq_ << std::endl;
				return 0;
			}

			// Google test
			friend class ReferenceTest;
		};
#endif //SWIG

		static const uintSeqLen kMinDistToContigEnds = 50; // Minimum distance from reference sequence ends to be taken into account for statistics so mapping issues at the border are avoided

	private:
		const uintSeqLen kMaxNInFragmentSite = 50; // The maximum number of n that is allowed in a fragment site to be counted for the bias calculation as it is highly likely that the fragment length and GC from those sites are wrong and highly unlikely that real fragments are found in this regions
		const uintSeqLen kMinNToSplitContigs = 10; // At positions with minimum kMinNToSplitContigs N's scaffolds will be broken into contigs, which means kMinDistToContigEnds applies
		const uintSeqLen kMinNToReplaceNWithRepeat = 100; // If at least kMinNToReplaceNWithRepeat consecutive N's are in a reference sequence the function ReplaceN fills in a short repeated pattern for these N instead of random bases
		const uintErrorCount kMaxErrorsShownPerFile = 50;

		seqan::StringSet<seqan::CharString> reference_ids_;
		seqan::StringSet<seqan::Dna5String> reference_sequences_;

		inline seqan::CharString &RefId( uintRefSeqId n ){
			return utilities::at(reference_ids_, n);
		}

		inline seqan::Dna5String &RefSeq( uintRefSeqId n ){
			return utilities::at(reference_sequences_, n);
		}

#ifndef SWIG // This part is not needed for the python plotting
		// Exclusion regions variables
		std::atomic<uintRefSeqId> obtained_exclusion_regions_for_num_sequences_;
		std::atomic<uintRefSeqId> cleared_exclusion_regions_for_num_sequences_;
		std::vector<std::vector<std::pair<uintSeqLen,uintSeqLen>>> excluded_regions_; // excluded_regions_[RefSeqId][PositionSortedExclusionRegionID] = {Start,End}
		std::vector<uintSeqLen> sum_of_excluded_bases_;

		// Variants variables
		uintAlleleId num_alleles_;
		seqan::VcfFileIn vcf_file_;
		seqan::VcfRecord cur_vcf_record_; // Already read but not yet handled vcf record (mostly the first of the next reference sequence, that is not yet prepared)
		std::atomic<uintRefSeqId> read_variation_for_num_sequences_;
		std::atomic<uintRefSeqId> cleared_variation_for_num_sequences_;
		std::vector<std::vector<Reference::Variant>> variants_; // variants_[RefSeqId][PositionSortedVariantID] = {Position,VarSeq,[Yes/No per allele]}
		std::vector<std::vector<uintSeqLen>> variant_positions_; // variant_positions_[RefSeqId][PositionSortedVariantID] = VariantPosition

		// Methylation bisulfite conversion variables
		std::ifstream methylation_file_;
		std::string cur_methylation_line_;
		std::string cur_methylation_sequence_;
		std::vector<std::vector<std::vector<double>>> unmethylation_; // unmethylation_[RefSeqId][Allele][PositionSortedMethylationID] = 1 - Methylation
		std::vector<std::vector<std::pair<uintSeqLen,uintSeqLen>>> unmethylated_regions_; // unmethylated_regions_[RefSeqId][PositionSortedMethylationID] = {UnMethylatedRegionStart, UnMethylatedRegionEnd}
		std::atomic<uintRefSeqId> read_methylation_for_num_sequences_;
		std::atomic<uintRefSeqId> cleared_methylation_for_num_sequences_;

		// Exclusion regions private functions
		inline void PushBackExclusionRegion( std::pair<uintSeqLen, uintSeqLen> region, uintRefSeqId ref_seq, uintSeqLen maximum_fragment_length);

		// Variants private functions
		bool CheckVcf() const;

		template<typename T> void InsertVariant(uintRefSeqId ref_seq_id, uintSeqLen position, const T &var_seq, uintAlleleBitArray allele){
			if(0 == variants_.at(ref_seq_id).size()){
				variants_.at(ref_seq_id).emplace_back(position, var_seq, allele);
			}
			else{
				// Check if variant is already added
				auto insert_it=variants_.at(ref_seq_id).end(); // Have variant at same position sorted by length: deletion/substitution/insertion(by length)
				auto var=variants_.at(ref_seq_id).size();
				while(0 < var && variants_.at(ref_seq_id).at(--var).position_ == position){
					if(variants_.at(ref_seq_id).at(var).var_seq_ == var_seq){
						// Variation is already in, so only adjust allele
						variants_.at(ref_seq_id).at(var).allele_ = variants_.at(ref_seq_id).at(var).allele_ | allele;
						return;
					}
					else if(length(variants_.at(ref_seq_id).at(var).var_seq_) > length(var_seq)){
						--insert_it;
					}
				}

				// If variant has not been added yet
				variants_.at(ref_seq_id).emplace(insert_it, position, var_seq, allele);
			}
		}

		inline bool ReadFirstVcfRecord();
		bool ReadVariants(uintRefSeqId end_ref_seq_id, bool positions_only);
		void CloseVcfFile();

		// Methylation bisulfite conversion private functions
		void CloseMethylationFile();
#endif //SWIG

		inline void UpdateGC( uintSeqLen &gc, uintSeqLen &n_count, const seqan::Dna5String &ref_seq, uintSeqLen old_pos, uintSeqLen new_pos ) const{
			// Account for new position in fragment
			if( utilities::IsGC(ref_seq, new_pos) ){
				gc += 1;
			}
			else if( utilities::IsN(ref_seq, new_pos) ){
				++n_count;
			}

			// Account for removed position from fragment
			if( utilities::IsGC(ref_seq, old_pos) ){
				gc -= 1;
			}
			else if( utilities::IsN(ref_seq, old_pos) ){
				--n_count;
			}
		}

		static inline double Bias( double general, double gc, double start_sur, double end_sur ){
			return general * gc * start_sur * end_sur;
		}

		void AddFragmentSite(std::vector<FragmentSite> &sites, uintSeqLen fragment_length, uintSeqLen gc, uintSeqLen n_count, const Surrounding &start_sur, const Surrounding &end_sur) const;

		// Google test
		friend class ReferenceTest;
		friend class FragmentDistributionStatsTest;
		friend class SimulatorTest;
		friend class SurroundingTest;
	public:
		Reference();

		inline uintRefSeqId NumberSequences() const{
			return length(reference_ids_);
		}

		inline const seqan::CharString &ReferenceId( uintRefSeqId n ) const{
			return utilities::at(reference_ids_, n);
		}

		const seqan::Prefix<const seqan::CharString>::Type ReferenceIdFirstPart( uintRefSeqId n ) const;

		inline const seqan::Dna5String &ReferenceSequence( uintRefSeqId n ) const{
			return utilities::at(reference_sequences_, n);
		}

		void ReferenceSequence(seqan::DnaString &insert_string, uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen frag_length, bool reversed=false) const;
#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the nested classes
		void ReferenceSequence(seqan::DnaString &insert_string, uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen frag_length, bool reversed, const std::vector<Variant> &variants, std::pair<intVariantId, uintSeqLen> first_variant, uintAlleleId allele) const;
#endif //SWIG

		inline const seqan::Dna5String &operator[]( uintRefSeqId n ) const{
			return ReferenceSequence(n);
		}

		uint32_t SequenceLength( uintRefSeqId n ) const{
			return length(ReferenceSequence(n));
		}

		uintSeqLen MaxSequenceLength() const{
			uintSeqLen max(0);
			for(uintRefSeqId seq=0; seq < NumberSequences(); ++seq){
				utilities::SetToMax(max, SequenceLength(seq));
			}
			return max;
		}

		uintRefLenCalc TotalSize() const;

		void RefSeqsSortedByNxx(std::vector<std::pair<double, uintRefSeqId>> &ref_seqs) const;
		uintRefSeqId RefSeqsInNxx(std::vector<bool> &ref_seqs, double nxx_ref_seqs) const;
		uintRefSeqBin NumRefSeqBinsInNxx(const std::vector<bool> &ref_seqs, uintSeqLen max_ref_seq_bin_length) const;

		inline bool GC( uintRefSeqId seq_id, uintSeqLen pos ) const{
			return utilities::IsGC( ReferenceSequence(seq_id), pos );
		}

		inline bool N( uintRefSeqId seq_id, uintSeqLen pos ) const{
			return utilities::IsN( ReferenceSequence(seq_id), pos );
		}

		inline uintSeqLen GCContentAbsolut( uintSeqLen &n_count, uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen end_pos ) const{
			uintSeqLen gc = 0;
			n_count = 0;

			const auto &seq(ReferenceSequence(seq_id));
			for( auto i = start_pos; i < end_pos; ++i ){
				if( utilities::IsGC(seq, i) ){
					++gc;
				}
				else if( utilities::IsN(seq, i) ){
					++n_count;
				}
			}

			return gc;
		}

		inline uintSeqLen GCContentAbsolut( uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen end_pos ) const{
			uintSeqLen n_count(0);
			return GCContentAbsolut( n_count, seq_id, start_pos, end_pos );
		}

		inline void UpdateGC( uintSeqLen &gc, uintSeqLen &n_count, uintRefSeqId seq_id, uintSeqLen old_pos, uintSeqLen new_pos ) const{
			UpdateGC(gc, n_count, ReferenceSequence(seq_id), old_pos, new_pos);
		}

		inline unsigned char GCContent( uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen end_pos ) const{ // end_pos like basically everywhere else is the position after the last element that still should be counted
			uintSeqLen n_count(0);
			auto gc = GCContentAbsolut(n_count, seq_id, start_pos, end_pos);

			return utilities::SafePercent(gc, end_pos-start_pos-n_count);
		}

		void ForwardSurrounding(Surrounding &surrounding, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			surrounding.Forward(ReferenceSequence(ref_seq_id), pos);
		}

		void ReverseSurrounding(Surrounding &surrounding, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			surrounding.Reverse(ReferenceSequence(ref_seq_id), pos);
		}

		void ForwardSurroundingWithN(Surrounding &surrounding, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			surrounding.ForwardWithN(ReferenceSequence(ref_seq_id), pos);
		}

		void ReverseSurroundingWithN(Surrounding &surrounding, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			surrounding.ReverseWithN(ReferenceSequence(ref_seq_id), pos);
		}

#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the decltype in Vect and SeqQualityStats
		static inline double Bias( double normalization, double ref_seq, double fragment_length, double gc, double start_sur, double end_sur ){
			return Bias(ref_seq*fragment_length*normalization, gc, start_sur, end_sur);
		}
		static inline double Bias( double ref_seq, double fragment_length, double gc, double start_sur, double end_sur ){
			return Bias(ref_seq*fragment_length, gc, start_sur, end_sur);
		}
		double SumBias(double &max_bias, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, double general_bias, const Vect<double> &gc_bias, const SurroundingBias &sur_bias) const;
		double SumBias(uintRefSeqId ref_seq_id, uintSeqLen fragment_length, double general_bias, const Vect<double> &gc_bias, const SurroundingBias &sur_bias) const;
		void GetFragmentSites(std::vector<FragmentSite> &sites, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintSeqLen start, uintSeqLen end) const;
#endif //SWIG

		bool ReadFasta(const char *fasta_file);
		void ReplaceN( uintSeed seed );
		bool HasN() const;
		bool WriteFasta(const char *fasta_file) const;

#ifndef SWIG // This part is not needed for the python plotting
		// Exclusion regions public functions
		void PrepareExclusionRegions();
		void ObtainExclusionRegions( uintRefSeqId end_ref_seq_id, uintSeqLen maximum_fragment_length );
		inline bool ObtainedExclusionRegionsForSequence(uintRefSeqId ref_seq_id) const{ return obtained_exclusion_regions_for_num_sequences_ > ref_seq_id; }
		inline bool ExclusionRegionsCompletelyObtained() const{ return obtained_exclusion_regions_for_num_sequences_ >= NumberSequences(); }
		bool FragmentExcluded( uintSeqLen &last_region_id, uintRefSeqId &last_ref_seq, uintRefSeqId ref_seq_id, uintSeqLen fragment_start, uintSeqLen fragment_end ) const;
		bool ReferenceSequenceExcluded( uintRefSeqId ref_seq_id ) const{ return 1 == excluded_regions_.at(ref_seq_id).size(); }
		uintRefSeqId FirstNotExcludedSequence() const{
			uintRefSeqId first = 0;
			while( first < NumberSequences() && ReferenceSequenceExcluded(first) ){
				++first;
			}
			return first;
		}
		uintSeqLen SumExcludedBases( uintRefSeqId ref_seq_id ) const{
			return sum_of_excluded_bases_.at(ref_seq_id);
		}
		uintRefLenCalc TotalExcludedBases() const{
			uintRefLenCalc total(0);
			for( auto excluded : sum_of_excluded_bases_ ){
				total += excluded;
			}
			return total;
		}
		uintSeqLen NumExcludedRegions( uintRefSeqId ref_seq_id ) const{
			return excluded_regions_.at(ref_seq_id).size();
		}
		uintSeqLen StartExclusion( uintRefSeqId ref_seq_id ) const{
			return excluded_regions_.at(ref_seq_id).front().second;
		}
		uintSeqLen EndExclusion( uintRefSeqId ref_seq_id ) const{
			return excluded_regions_.at(ref_seq_id).back().first;
		}
		uintSeqLen NextExclusionRegion(uintRefSeqId ref_seq_id, uintSeqLen position) const{
			uintSeqLen next_exclusion_region(1);
			while(position > excluded_regions_.at(ref_seq_id).at(next_exclusion_region).first){
				++next_exclusion_region;
			}
			return next_exclusion_region;
		}
		void CorrectPositionAddingExcludedRegions(uintSeqLen &start_pos, uintSeqLen &next_exclusion_region, uintRefSeqId ref_seq) const{
			while( start_pos >= excluded_regions_.at(ref_seq).at(next_exclusion_region).first ){ // We don't reach the end of the list because we don't test for the last bin (start of last bin is tested on second to last bin)
				start_pos += excluded_regions_.at(ref_seq).at(next_exclusion_region).second - excluded_regions_.at(ref_seq).at(next_exclusion_region).first;
				++next_exclusion_region;
			}
		}
		void CorrectPositionRemovingExcludedRegions(uintSeqLen &correction, uintSeqLen &cur_exclusion_region, uintRefSeqId ref_seq_id, uintSeqLen position) const{
			while( cur_exclusion_region+1 < NumExcludedRegions(ref_seq_id) && position > excluded_regions_.at(ref_seq_id).at(cur_exclusion_region+1).first ){
				++cur_exclusion_region;
				correction += excluded_regions_.at(ref_seq_id).at(cur_exclusion_region).second - excluded_regions_.at(ref_seq_id).at(cur_exclusion_region).first;
			}

			while( cur_exclusion_region && position < excluded_regions_.at(ref_seq_id).at(cur_exclusion_region).second ){
				correction -= excluded_regions_.at(ref_seq_id).at(cur_exclusion_region).second - excluded_regions_.at(ref_seq_id).at(cur_exclusion_region).first;
				--cur_exclusion_region;
			}
		}
		inline void ClearExclusionRegions( uintRefSeqId end_ref_seq_id ){
			while( cleared_exclusion_regions_for_num_sequences_ < end_ref_seq_id ){
				excluded_regions_.at(cleared_exclusion_regions_for_num_sequences_).clear();
				excluded_regions_.at(cleared_exclusion_regions_for_num_sequences_++).shrink_to_fit();
			}
		}
		inline void ClearAllExclusionRegions(){
			excluded_regions_.clear();
			excluded_regions_.shrink_to_fit();
		}

		// Variants public functions
		inline uintAlleleId NumAlleles() const{ return num_alleles_; }
		inline bool VariantsLoaded() const{ return variants_.size(); }
		inline bool VariantPositionsLoaded() const{ return variant_positions_.size(); }
		inline bool VariantsLoadedForSequence(uintRefSeqId ref_seq_id) const{ return read_variation_for_num_sequences_ > ref_seq_id; }
		inline bool VariantPositionsLoadedForSequence(uintRefSeqId ref_seq_id) const{ return read_variation_for_num_sequences_ >= ref_seq_id; }
		inline bool VariantsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= variants_.size(); }
		inline bool VariantPositionsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= variant_positions_.size(); }
		inline const std::vector<Reference::Variant> &Variants(uintRefSeqId ref_seq_id) const{ return variants_.at(ref_seq_id); }
		inline const std::vector<uintRefSeqId> &VariantPositions(uintRefSeqId ref_seq_id) const{ return variant_positions_.at(ref_seq_id); }

		bool PrepareVariantFile(const std::string &var_file);
		inline bool ReadVariants(uintRefSeqId end_ref_seq_id){ return ReadVariants(end_ref_seq_id, false); };
		inline bool ReadVariantPositions(uintRefSeqId end_ref_seq_id){ return ReadVariants(end_ref_seq_id, true); };
		bool ReadFirstVariants();
		bool ReadFirstVariantPositions();
		void ClearVariants(uintRefSeqId end_ref_seq_id);
		void ClearVariantPositions(uintRefSeqId end_ref_seq_id);
		void ClearAllVariants();
		void ClearAllVariantPositions();

		// Methylation bisulfite conversion public functions
		inline bool MethylationLoaded() const{ return unmethylation_.size(); }
		inline bool MethylationLoadedForSequence(uintRefSeqId ref_seq_id) const{ return read_methylation_for_num_sequences_ > ref_seq_id; }
		inline bool MethylationCompletelyLoaded() const{ return read_methylation_for_num_sequences_ >= unmethylation_.size(); }
		inline const std::vector<double> &Unmethylation(uintRefSeqId ref_seq_id, uintAlleleId allele) const{
			if(1 < unmethylation_.at(ref_seq_id).size()){
				return unmethylation_.at(ref_seq_id).at(allele);
			}
			else{
				return unmethylation_.at(ref_seq_id).at(0);
			}
		}
		inline const std::vector<std::pair<uintSeqLen,uintSeqLen>> &UnmethylatedRegions(uintRefSeqId ref_seq_id) const{
			return unmethylated_regions_.at(ref_seq_id);
		}

		bool PrepareMethylationFile(const std::string &methylation_file);
		bool ReadMethylation(uintRefSeqId end_ref_seq_id);
		void ClearMethylation(uintRefSeqId end_ref_seq_id);
		void ClearAllMethylation();
#endif //SWIG
	};

	struct FragmentSite{
		uintPercent gc_;
		Surrounding start_surrounding_;
		Surrounding end_surrounding_;
		uintDupCount count_forward_;
		uintDupCount count_reverse_;
		double bias_;

		FragmentSite(uintPercent gc, const Surrounding &start_sur, const Surrounding &end_sur):
			gc_(gc),
			start_surrounding_(start_sur),
			end_surrounding_(end_sur),
			count_forward_(0),
			count_reverse_(0),
			bias_(1.0)
		{}

		FragmentSite(): // Placeholder
			gc_(0),
			count_forward_(0),
			count_reverse_(0),
			bias_(0.0)
		{}

		static bool LowerBias(const FragmentSite &lhs, const FragmentSite &rhs){
			return lhs.bias_ < rhs.bias_;
		}
	};
}

#endif // REFERENCE_H
