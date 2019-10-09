#ifndef REFERENCE_H
#define REFERENCE_H

#include <atomic>
#include <array>
#include <set>
#include <stdint.h>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

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

		static const uintSurBlockId num_surrounding_blocks_ = 3; // Number of surrounding blocks
		static const uintSeqLen surrounding_range_ = 10; // Length of a single surrounding block
		static const int16_t surrounding_start_pos_ = -10; // Position where the surrounding should start relative to the first base in the fragment (so -10 is 10 bases before the fragment) [-num_surrounding_blocks_*surrounding_range_ < surrounding_start_pos_ <= 0]
		static uintSeqLen SurroundingSize(){
			return 1 << 2*surrounding_range_;
		}

		static const uintSeqLen kMinDistToRefSeqEnds = 50; // Minimum distance from reference sequence ends to be taken into account for statistics so mapping issues at the border are avoided

	private:
		const uintSeqLen kMaxNInFragmentSite = 50; // The maximum number of n that is allowed in a fragment site to be counted for the bias calculation as it is highly likely that the fragment length and GC from those sites are wrong and highly unlikely that real fragments are found in this regions
		const uintErrorCount kMaxErrorsShownPerFile = 50;

		static seqan::CharString dummy_charstring_;
		static seqan::Dna5String dummy_dna5string_;

		seqan::StringSet<seqan::CharString> reference_ids_;
		seqan::StringSet<seqan::Dna5String> reference_sequences_;

		inline seqan::CharString &RefId( uintRefSeqId n ){
			if( length(reference_ids_) > n ){
				return reference_ids_[n];
			}
			else{
				printErr << "Called reference id " << n << ". Size is " << length(reference_ids_) << std::endl;
				throw std::out_of_range( "Non-existing reference id called" );
				return dummy_charstring_;
			}
		}

		inline seqan::Dna5String &RefSeq( uintRefSeqId n ){
			if( length(reference_sequences_) > n ){
				return reference_sequences_[n];
			}
			else{
				printErr << "Called reference sequence " << n << ". Size is " << length(reference_sequences_) << std::endl;
				throw std::out_of_range( "Non-existing reference sequence called" );
				return dummy_dna5string_;
			}
		}

#ifndef SWIG // This part is not needed for the python plotting
		uintAlleleId num_alleles_;
		seqan::VcfFileIn vcf_file_;
		seqan::VcfRecord cur_vcf_record_; // Already read but not yet handled vcf record (mostly the first of the next reference sequence, that is not yet prepared)
		std::atomic<uintRefSeqId> read_variation_for_num_sequences_;
		std::atomic<uintRefSeqId> cleared_variation_for_num_sequences_;
		std::vector<std::vector<Reference::Variant>> variants_; // variants_[RefSeqId][PositionSortedVariantID] = {Position,VarSeq,[Yes/No per allele]}
		std::vector<std::vector<uintSeqLen>> variant_positions_; // variant_positions_[RefSeqId][PositionSortedVariantID] = VariantPosition

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

		intSurrounding ForwardSurroundingBlock( uintRefSeqId ref_seq_id, uintSeqLen sur_start_pos, uintSurBlockId block ) const;
		intSurrounding ReverseSurroundingBlock( uintRefSeqId ref_seq_id, uintSeqLen sur_start_pos, uintSurBlockId block ) const;

		inline uintBaseCall NewForwardSurroundingBase(const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_start, uintSurBlockId block) const{
			return ref_seq[(length(ref_seq) + (block+1)*surrounding_range_ + new_fragment_start-1 + surrounding_start_pos_)%length(ref_seq)];
		}

		inline uintBaseCall NewReverseSurroundingBase(const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_end, uintSurBlockId block) const{
			return utilities::Complement::Dna5(ref_seq[(length(ref_seq) + new_fragment_end - surrounding_start_pos_ - block*surrounding_range_)%length(ref_seq)]);
		}

		static inline double Bias( double general, double gc, const std::array<double, num_surrounding_blocks_> &start_sur, const std::array<double, num_surrounding_blocks_> &end_sur ){
			double bias = general * gc;
			double sur_start(0.0), sur_end(0.0);
			for(auto block = num_surrounding_blocks_; block--;){
				sur_start += start_sur.at(block);
				sur_end += end_sur.at(block);
			}
			bias *= utilities::InvLogit2(sur_start) * utilities::InvLogit2(sur_end);
			return  bias;
		}

		void AddFragmentSite(std::vector<FragmentSite> &sites, uintSeqLen fragment_length, uintSeqLen gc, uintSeqLen n_count, const std::array<intSurrounding, num_surrounding_blocks_> &start_sur, const std::array<intSurrounding, num_surrounding_blocks_> &end_sur) const;

		// Google test
		friend class ReferenceTest;
		friend class SimulatorTest;
		FRIEND_TEST(ReferenceTest, Functionality);
	public:
		Reference();

		inline uint16_t SurroundingRange() const{
			return surrounding_range_;
		}
		inline uint16_t SurroundingStartPos() const{
			return surrounding_start_pos_;
		}

		inline uintRefSeqId NumberSequences() const{
			return length(reference_ids_);
		}

		inline const seqan::CharString &ReferenceId( uintRefSeqId n ) const{
			if( length(reference_ids_) > n ){
				return reference_ids_[n];
			}
			else{
				printErr << "Called reference id " << n << ". Size is " << length(reference_ids_) << std::endl;
				throw std::out_of_range( "Non-existing reference id called" );
				return dummy_charstring_;
			}
		}

		const seqan::Prefix<const seqan::CharString>::Type ReferenceIdFirstPart( uintRefSeqId n ) const;

		inline const seqan::Dna5String &ReferenceSequence( uintRefSeqId n ) const{
			if( length(reference_sequences_) > n ){
				return reference_sequences_[n];
			}
			else{
				printErr << "Called reference sequence " << n << ". Size is " << length(reference_sequences_) << std::endl;
				throw std::out_of_range( "Non-existing reference sequence called" );
				return dummy_dna5string_;
			}
		}

		void ReferenceSequence(seqan::DnaString &insert_string, uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen min_length, bool reversed=false, uintSeqLen seq_size=0) const;
#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the nested classes
		void ReferenceSequence(seqan::DnaString &insert_string, uintRefSeqId seq_id, uintSeqLen start_pos, uintSeqLen min_length, bool reversed, const std::vector<Variant> &variants, std::pair<intVariantId, uintSeqLen> first_variant, uintAlleleId allele, uintSeqLen seq_size=0) const;
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

			return reseq::utilities::SafePercent(gc, end_pos-start_pos-n_count);
		}

		void ForwardSurrounding(std::array<intSurrounding, num_surrounding_blocks_> &surroundings, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			// Other blocks
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				// Add SequenceLength(ref_seq_id) to avoid negative positions if block reaches over lower sequence ends so it wraps around and start at the end
				// The higher sequence end is covered by the called function as long as the given position is still a valid position and not higher or equal to the length
				surroundings.at(block) = ForwardSurroundingBlock(ref_seq_id, pos, block);
			}
		}

		void ForwardSurroundingWithN(std::array<intSurrounding, num_surrounding_blocks_> &surroundings, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			const seqan::Dna5String &ref_seq(ReferenceSequence(ref_seq_id));
			uintBaseCall new_value;
			uintSeqLen block_pos(pos+surrounding_start_pos_);
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				surroundings.at(block) = 0;
				for( auto ref_pos = block_pos; ref_pos < block_pos+surrounding_range_; ++ref_pos ){
					new_value = ref_seq[ref_pos];
					if(utilities::IsN(new_value)){
						surroundings.at(block) = static_cast<intSurrounding>(block_pos)-ref_pos-1;
						for( auto ref_pos2 = block_pos+surrounding_range_; --ref_pos2 > ref_pos; ){
							if(utilities::IsN(ref_seq, ref_pos2)){
								surroundings.at(block) = static_cast<intSurrounding>(block_pos)-ref_pos2-1;
								break; // Last N is only interesting one
							}
						}
						break; // We don't need to calculate surrounding anymore as it is invalidated by N
					}
					else{
						surroundings.at(block) = (surroundings.at(block) << 2) + new_value;
					}
				}
				block_pos += surrounding_range_;
			}
		}

		void ReverseSurrounding(std::array<intSurrounding, num_surrounding_blocks_> &surroundings, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				surroundings.at(block) = ReverseSurroundingBlock(ref_seq_id, pos, block);
			}
		}

		void ReverseSurroundingWithN(std::array<intSurrounding, num_surrounding_blocks_> &surroundings, uintRefSeqId ref_seq_id, uintSeqLen pos ) const{
			utilities::ComplementedConstDna5String ref_seq(ReferenceSequence(ref_seq_id));
			uintBaseCall new_value;
			uintSeqLen block_pos(pos+1-surrounding_start_pos_-surrounding_range_);
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				surroundings.at(block) = 0;
				for( auto ref_pos = block_pos+surrounding_range_; ref_pos-- > block_pos; ){
					new_value = ref_seq[ref_pos];
					if(utilities::IsN(new_value)){
						surroundings.at(block) = static_cast<intSurrounding>(block_pos)-ref_pos-1;
						break; // We don't need to calculate surrounding anymore as it is invalidated by N and we already have last N
					}
					else{
						surroundings.at(block) = (surroundings.at(block) << 2) + new_value;
					}
				}
				block_pos -= surrounding_range_;
			}
		}

		inline void UpdateForwardSurrounding( std::array<intSurrounding, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_start ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = (sur.at(block) << 2)%SurroundingSize() + NewForwardSurroundingBase(ref_seq, new_fragment_start, block);
			}
		}

		inline void UpdateForwardSurroundingWithN( std::array<intSurrounding, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_start, uintSeqLen ref_seq_id ) const{
			uintBaseCall new_base;
			for(auto block=sur.size(); block--; ){
				new_base = NewForwardSurroundingBase(ref_seq, new_fragment_start, block);

				if(utilities::IsN(new_base)){
					sur.at(block) = -surrounding_range_; // Wait surrounding_range_ positions until calculating surrounding again
				}
				else if(0 > sur.at(block)){
					if(-1 == sur.at(block)){
						// Calculate current surrounding
						sur.at(block) = ForwardSurroundingBlock(ref_seq_id, new_fragment_start, block);
					}
					else{
						++(sur.at(block));
					}
				}
				else{
					// No N: Update surrounding normally
					sur.at(block) = (sur.at(block) << 2)%SurroundingSize() + new_base;
				}
			}
		}

		inline void UpdateReverseSurrounding( std::array<intSurrounding, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_end ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = (sur.at(block) + NewReverseSurroundingBase(ref_seq, new_fragment_end, block)*SurroundingSize())/4;
			}
		}

		inline void UpdateReverseSurroundingWithN( std::array<intSurrounding, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_end, uintSeqLen ref_seq_id ) const{
			uintBaseCall new_base;
			for(auto block=sur.size(); block--; ){
				new_base = NewReverseSurroundingBase(ref_seq, new_fragment_end, block);

				if(utilities::IsN(new_base)){
					sur.at(block) = -surrounding_range_; // Wait surrounding_range_ positions until calculating surrounding again
				}
				else if(0 > sur.at(block)){
					if(-1 == sur.at(block)){
						// Calculate current surrounding
						sur.at(block) = ReverseSurroundingBlock(ref_seq_id, new_fragment_end, block);
					}
					else{
						++(sur.at(block));
					}
				}
				else{
					// No N: Update surrounding normally
					sur.at(block) = (sur.at(block) + new_base*SurroundingSize())/4;
				}
			}
		}

		inline void RollBackReverseSurrounding( std::array<intSurrounding, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, uintSeqLen new_fragment_end ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = ((sur.at(block)<<2) + static_cast<uintBaseCall>(utilities::Complement::Dna5(ref_seq[(length(ref_seq)+new_fragment_end-surrounding_start_pos_-(block+1)*surrounding_range_+1)%length(ref_seq)])))%SurroundingSize();
			}
		}

#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the decltype in Vect and SeqQualityStats
		static inline double Bias( double normalization, double ref_seq, double fragment_length, double gc, const std::array<double, num_surrounding_blocks_> &start_sur, const std::array<double, num_surrounding_blocks_> &end_sur ){
			return Bias(ref_seq*fragment_length*normalization, gc, start_sur, end_sur);
		}
		static inline double Bias( double ref_seq, double fragment_length, double gc, const std::array<double, num_surrounding_blocks_> &start_sur, const std::array<double, num_surrounding_blocks_> &end_sur ){
			return Bias(ref_seq*fragment_length, gc, start_sur, end_sur);
		}
		double SumBias(std::vector<utilities::VectorAtomic<uintFragCount>> &gc_sites, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, double general_bias, const Vect<double> &gc_bias, const std::array<std::vector<double>, num_surrounding_blocks_> &sur_bias) const;
		double SumBias(double &max_bias, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, double general_bias, const Vect<double> &gc_bias, const std::array<std::vector<double>, num_surrounding_blocks_> &sur_bias) const;
		void GetFragmentSites(std::vector<FragmentSite> &sites, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintSeqLen start, uintSeqLen end) const;
#endif //SWIG

		bool ReadFasta(const char *fasta_file);
		void ReplaceN( uintSeed seed );
		bool HasN() const;
		bool WriteFasta(const char *fasta_file) const;

#ifndef SWIG // This part is not needed for the python plotting
		inline uintAlleleId NumAlleles() const{ return num_alleles_; }
		inline bool VariantsLoaded() const{ return variants_.size(); }
		inline bool VariantPositionsLoaded() const{ return variant_positions_.size(); }
		inline bool VariantsLoadedForSequence(uintRefSeqId ref_seq_id) const{ return read_variation_for_num_sequences_ >= ref_seq_id; }
		inline bool VariantPositionsLoadedForSequence(uintRefSeqId ref_seq_id) const{ return read_variation_for_num_sequences_ >= ref_seq_id; }
		inline bool VariantsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= NumberSequences(); }
		inline bool VariantPositionsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= NumberSequences(); }
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
#endif //SWIG
	};

	struct FragmentSite{
		uintPercent gc_;
		std::array<intSurrounding, Reference::num_surrounding_blocks_> start_surrounding_;
		std::array<intSurrounding, Reference::num_surrounding_blocks_> end_surrounding_;
		uintDupCount count_forward_;
		uintDupCount count_reverse_;
		double bias_;

		FragmentSite(uintPercent gc, const std::array<intSurrounding, Reference::num_surrounding_blocks_> &start_sur, const std::array<intSurrounding, Reference::num_surrounding_blocks_> &end_sur):
			gc_(gc),
			start_surrounding_(start_sur),
			end_surrounding_(end_sur),
			count_forward_(0),
			count_reverse_(0),
			bias_(1.0)
		{}

		FragmentSite(): // Placeholder
			gc_(0),
			start_surrounding_({0,0,0}),
			end_surrounding_({0,0,0}),
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
