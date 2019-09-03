#ifndef REFERENCE_H
#define REFERENCE_H

#include <atomic>
#include <array>
#include <set>
#include <stdint.h>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

#include "utilities.h"
#include "Vect.hpp"

namespace reseq{
	struct FragmentSite;

	class Reference{
	public:
#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the nested classes
		class Variant{
		public:
			uint32_t position_;
			seqan::DnaString var_seq_;
			uint64_t allele_;

			template<typename T> Variant(uint32_t position, const T &var_seq, uint64_t allele):
				position_(position),
				var_seq_(var_seq),
				allele_(allele){
			}

			inline bool InAllele(uint16_t allele) const{
				return (allele_ >> allele) & 1; // Get bit corresponding to allele
			}

			inline uint16_t FirstAllele() const{
				if(allele_){
					uint16_t first = 0;
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

		static const uint16_t num_surrounding_blocks_ = 3; // Number of surrounding blocks
		static const uint16_t surrounding_range_ = 10; // Length of a single surrounding block
		static const int16_t surrounding_start_pos_ = -10; // Position where the surrounding should start relative to the first base in the fragment (so -10 is 10 bases before the fragment) [-num_surrounding_blocks_*surrounding_range_ < surrounding_start_pos_ <= 0]
		static uint32_t SurroundingSize(){
			return 1 << 2*surrounding_range_;
		}

#ifndef SWIG // This part is not needed for the python plotting
		static seqan::FunctorComplement<seqan::Dna5> Complementor;
		static seqan::Dna5String ReverseComplementor(seqan::Dna5String org){
			return seqan::ModifiedString< seqan::ModifiedString<const seqan::Dna5String, seqan::ModComplementDna5>, seqan::ModReverse>(org);
		}
#endif //SWIG

	private:
		const uint16_t min_dist_to_ref_seq_ends_; // Minimum distance from reference sequence ends to be taken into account for statistics so mapping issues at the border are avoided
		const uint32_t max_n_in_fragment_site_; // The maximum number of n that is allowed in a fragment site to be counted for the bias calculation as it is highly likely that the fragment length and GC from those sites are wrong and highly unlikely that real fragments are found in this regions

		seqan::StringSet<seqan::CharString> reference_ids_;
		seqan::StringSet<seqan::Dna5String> reference_sequences_;

#ifndef SWIG // This part is not needed for the python plotting
		uint16_t num_alleles_;
		seqan::VcfFileIn vcf_file_;
		seqan::VcfRecord cur_vcf_record_; // Already read but not yet handled vcf record (mostly the first of the next reference sequence, that is not yet prepared)
		std::atomic<uint32_t> read_variation_for_num_sequences_;
		std::atomic<uint32_t> cleared_variation_for_num_sequences_;
		std::vector<std::vector<Reference::Variant>> variants_; // variants_[RefSeqId][PositionSortedVariantID] = {Position,VarSeq,[Yes/No per allele]}
		std::vector<std::vector<uint32_t>> variant_positions_; // variant_positions_[RefSeqId][PositionSortedVariantID] = VariantPosition

		bool CheckVcf() const;

		template<typename T> void InsertVariant(uint32_t ref_seq_id, uint32_t position, const T &var_seq, uint64_t allele){
			if(0 == variants_.at(ref_seq_id).size()){
				variants_.at(ref_seq_id).emplace_back(position, var_seq, allele);
			}
			else{
				// Check if variant is already added
				auto insert_it=variants_.at(ref_seq_id).end(); // Have variant at same position sorted by length: deletion/substitution/insertion(by length)
				uint32_t var=variants_.at(ref_seq_id).size();
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
		bool ReadVariants(uint32_t end_ref_seq_id, bool positions_only);
		void CloseVcfFile();
#endif //SWIG

		inline void UpdateGC( seqan::Size<seqan::Dna5String>::Type &gc, uint32_t &n_count, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type old_pos, seqan::Size<seqan::Dna5String>::Type new_pos ) const{
			// Account for new position in fragment
			if( 1 == ref_seq[new_pos] || 2 == ref_seq[new_pos] ){
				gc += 1;
			}
			else if( 3 < ref_seq[new_pos]){
				++n_count;
			}

			// Account for removed position from fragment
			if( 1 == ref_seq[old_pos] || 2 == ref_seq[old_pos] ){
				gc -= 1;
			}
			else if( 3 < ref_seq[old_pos]){
				--n_count;
			}
		}

		uint32_t ForwardSurroundingBlock( seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, uint64_t pos ) const;
		uint32_t ReverseSurroundingBlock( seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, uint64_t pos ) const;

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

		void AddFragmentSite(std::vector<FragmentSite> &sites, uint32_t fragment_length, uint32_t gc, uint32_t n_count, const std::array<int32_t, num_surrounding_blocks_> &start_sur, const std::array<int32_t, num_surrounding_blocks_> &end_sur) const;

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
		inline uint16_t MinDistToRefSeqEnds() const{
			return min_dist_to_ref_seq_ends_;
		}

		inline seqan::Size<decltype(reference_ids_)>::Type NumberSequences() const{
			return length(reference_ids_);
		}

		inline const seqan::CharString &ReferenceId( seqan::Size<decltype(reference_ids_)>::Type n ) const{
			return reference_ids_[n];
		}

		const seqan::Prefix<const seqan::CharString>::Type ReferenceIdFirstPart( seqan::Size<decltype(reference_ids_)>::Type n ) const;

		inline const seqan::Dna5String &ReferenceSequence( seqan::Size<decltype(reference_sequences_)>::Type n ) const{
			return reference_sequences_[n];
		}

		void ReferenceSequence(seqan::DnaString &insert_string, seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type start_pos, seqan::Size<seqan::Dna5String>::Type min_length, bool reversed=false, seqan::Size<seqan::Dna5String>::Type seq_size=0) const;
#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the nested classes
		void ReferenceSequence(seqan::DnaString &insert_string, seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type start_pos, seqan::Size<seqan::Dna5String>::Type min_length, bool reversed, const std::vector<Variant> &variants, std::pair<int32_t, uint16_t> first_variant, uint16_t allele, seqan::Size<seqan::Dna5String>::Type seq_size=0) const;
#endif //SWIG

		inline const seqan::Dna5String &operator[]( seqan::Size<decltype(reference_sequences_)>::Type n ) const{
			return reference_sequences_[n];
		}

		seqan::Size<seqan::Dna5String>::Type SequenceLength( seqan::Size<decltype(reference_sequences_)>::Type n ) const{
			return length(reference_sequences_[n]);
		}

		seqan::Size<seqan::Dna5String>::Type MaxSequenceLength() const{
			seqan::Size<seqan::Dna5String>::Type max(0);
			for(uint32_t seq=0; seq < NumberSequences(); ++seq){
				utilities::SetToMax(max, SequenceLength(seq));
			}
			return max;
		}

		uint64_t TotalSize() const;

		void RefSeqsSortedByNxx(std::vector<std::pair<double, uint32_t>> &ref_seqs) const;
		uint32_t RefSeqsInNxx(std::vector<bool> &ref_seqs, double nxx_ref_seqs) const;
		uint32_t NumRefSeqBinsInNxx(const std::vector<bool> &ref_seqs, uint64_t max_ref_seq_bin_length) const;

		inline bool GC( seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			return ( 'G' == reference_sequences_[seq_id][pos] || 'C' == reference_sequences_[seq_id][pos] );
		}

		inline bool N( seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			return ( 'N' == reference_sequences_[seq_id][pos] );
		}

		inline uint64_t GCContentAbsolut( uint32_t &n_count, seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type start_pos, seqan::Size<seqan::Dna5String>::Type end_pos ) const{
			uint64_t gc = 0;

			for( auto i = start_pos; i < end_pos; ++i ){
				if( GC(seq_id, i) ){
					++gc;
				}
				else if( N(seq_id, i) ){
					++n_count;
				}
			}

			return gc;
		}

		inline void UpdateGC( seqan::Size<seqan::Dna5String>::Type &gc, uint32_t &n_count, seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type old_pos, seqan::Size<seqan::Dna5String>::Type new_pos ) const{
			// Account for new position in fragment
			if( GC(seq_id, new_pos) ){
				gc += 1;
			}
			else if( N(seq_id, new_pos) ){
				++n_count;
			}

			// Account for removed position from fragment
			if( GC(seq_id, old_pos) ){
				gc -= 1;
			}
			else if( N(seq_id, old_pos) ){
				--n_count;
			}
		}

		inline unsigned char GCContent( seqan::Size<decltype(reference_sequences_)>::Type seq_id, seqan::Size<seqan::Dna5String>::Type start_pos, seqan::Size<seqan::Dna5String>::Type end_pos ) const{ // end_pos like basically everywhere else is the position after the last element that still should be counted
			uint32_t n_count(0);
			auto gc = GCContentAbsolut(n_count, seq_id, start_pos, end_pos);
			if( end_pos-start_pos > n_count){
				return reseq::utilities::Percent(gc, end_pos-start_pos-n_count);
			}
			else{
				return 50;
			}
		}

		void ForwardSurrounding(std::array<uint32_t, num_surrounding_blocks_> &surroundings, seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			// Other blocks
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				// Add SequenceLength(ref_seq_id) to avoid negative positions if block reaches over lower sequence ends so it wraps around and start at the end
				// The higher sequence end is covered by the called function as long as the given position is still a valid position and not higher or equal to the length
				surroundings.at(block) = ForwardSurroundingBlock(ref_seq_id, (SequenceLength(ref_seq_id)+pos+block*surrounding_range_+surrounding_start_pos_)%SequenceLength(ref_seq_id));
			}
		}

		void ForwardSurroundingWithN(std::array<int32_t, num_surrounding_blocks_> &surroundings, seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			const seqan::Dna5String &ref_seq(reference_sequences_[ref_seq_id]);
			uint16_t new_value;
			uint64_t block_pos(pos+surrounding_start_pos_);
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				surroundings.at(block) = 0;
				for( auto ref_pos = block_pos; ref_pos < block_pos+surrounding_range_; ++ref_pos ){
					new_value = static_cast<uint16_t>(ref_seq[ref_pos]);
					if(3 < new_value){
						// N
						surroundings.at(block) = static_cast<int32_t>(block_pos)-ref_pos-1;
						for( auto ref_pos2 = block_pos+surrounding_range_; --ref_pos2 > ref_pos; ){
							new_value = static_cast<uint16_t>(ref_seq[ref_pos2]);
							if(3 < new_value){
								surroundings.at(block) = static_cast<int32_t>(block_pos)-ref_pos2-1;
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

		void ReverseSurrounding(std::array<uint32_t, num_surrounding_blocks_> &surroundings,  seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				if(pos+1-surrounding_start_pos_ < (block+1)*surrounding_range_){
					surroundings.at(block) = ReverseSurroundingBlock(ref_seq_id, SequenceLength(ref_seq_id)+pos+1-surrounding_start_pos_-(block+1)*surrounding_range_);
				}
				else{
					surroundings.at(block) = ReverseSurroundingBlock(ref_seq_id, pos+1-surrounding_start_pos_-(block+1)*surrounding_range_);
				}
			}
		}

		void ReverseSurroundingWithN(std::array<int32_t, num_surrounding_blocks_> &surroundings,  seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type pos ) const{
			const seqan::ModifiedString<const seqan::Dna5String, seqan::ModComplementDna5> ref_seq(reference_sequences_[ref_seq_id]);
			uint16_t new_value;
			uint64_t block_pos(pos+1-surrounding_start_pos_-surrounding_range_);
			for(auto block=0; block < num_surrounding_blocks_; ++block){
				surroundings.at(block) = 0;
				for( auto ref_pos = block_pos+surrounding_range_; ref_pos-- > block_pos; ){
					new_value = static_cast<uint16_t>(ref_seq[ref_pos]);
					if(3 < new_value){
						// N
						surroundings.at(block) = static_cast<int32_t>(block_pos)-ref_pos-1;
						break; // We don't need to calculate surrounding anymore as it is invalidated by N and we already have last N
					}
					else{
						surroundings.at(block) = (surroundings.at(block) << 2) + new_value;
					}
				}
				block_pos -= surrounding_range_;
			}
		}

		inline void UpdateForwardSurrounding( std::array<uint32_t, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type new_fragment_start ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = (sur.at(block) << 2)%SurroundingSize() + static_cast<uint16_t>(ref_seq[(new_fragment_start-1+surrounding_start_pos_+(block+1)*surrounding_range_)%length(ref_seq)]);
			}
		}

		inline void UpdateForwardSurroundingWithN( std::array<int32_t, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type new_fragment_start, seqan::Size<seqan::Dna5String>::Type ref_seq_id ) const{
			uint16_t new_base;
			for(auto block=sur.size(); block--; ){
				new_base = ref_seq[new_fragment_start-1+surrounding_start_pos_+(block+1)*surrounding_range_];

				if(3 < new_base){
					// N
					sur.at(block) = -surrounding_range_; // Wait surrounding_range_ positions until calculating surrounding again
				}
				else if(0 > sur.at(block)){
					if(-1 == sur.at(block)){
						// Calculate current surrounding
						sur.at(block) = static_cast<int32_t>( ForwardSurroundingBlock(ref_seq_id, new_fragment_start+surrounding_start_pos_+block*surrounding_range_) );
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

		inline void UpdateReverseSurrounding( std::array<uint32_t, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type new_fragment_end ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = (sur.at(block) + static_cast<uint16_t>(Complementor(ref_seq[(length(ref_seq)+new_fragment_end-surrounding_start_pos_-block*surrounding_range_)%length(ref_seq)]))*SurroundingSize())/4;
			}
		}

		inline void UpdateReverseSurroundingWithN( std::array<int32_t, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type new_fragment_end, seqan::Size<seqan::Dna5String>::Type ref_seq_id ) const{
			uint16_t new_base;
			for(auto block=sur.size(); block--; ){
				new_base = ref_seq[new_fragment_end-surrounding_start_pos_-block*surrounding_range_];

				if(3 < new_base){
					// N
					sur.at(block) = -surrounding_range_; // Wait surrounding_range_ positions until calculating surrounding again
				}
				else if(0 > sur.at(block)){
					if(-1 == sur.at(block)){
						// Calculate current surrounding
						sur.at(block) = static_cast<int32_t>( ReverseSurroundingBlock(ref_seq_id, new_fragment_end+1-surrounding_start_pos_-(block+1)*surrounding_range_) );
					}
					else{
						++(sur.at(block));
					}
				}
				else{
					// No N: Update surrounding normally
					sur.at(block) = (sur.at(block) + static_cast<uint16_t>(Complementor(new_base))*SurroundingSize())/4;
				}
			}
		}

		inline void RollBackReverseSurrounding( std::array<uint32_t, num_surrounding_blocks_> &sur, const seqan::Dna5String &ref_seq, seqan::Size<seqan::Dna5String>::Type new_fragment_end ) const{
			for(auto block=sur.size(); block--; ){
				sur.at(block) = ((sur.at(block)<<2) + static_cast<uint16_t>(Complementor(ref_seq[(length(ref_seq)+new_fragment_end-surrounding_start_pos_-(block+1)*surrounding_range_+1)%length(ref_seq)])))%SurroundingSize();
			}
		}

#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the decltype in Vect and SeqQualityStats
		static inline double Bias( double normalization, double ref_seq, double fragment_length, double gc, const std::array<double, num_surrounding_blocks_> &start_sur, const std::array<double, num_surrounding_blocks_> &end_sur ){
			return Bias(ref_seq*fragment_length*normalization, gc, start_sur, end_sur);
		}
		static inline double Bias( double ref_seq, double fragment_length, double gc, const std::array<double, num_surrounding_blocks_> &start_sur, const std::array<double, num_surrounding_blocks_> &end_sur ){
			return Bias(ref_seq*fragment_length, gc, start_sur, end_sur);
		}
		double SumBias(std::vector<utilities::VectorAtomic<uint64_t>> &gc_sites, seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type fragment_length, double general_bias, const Vect<double> &gc_bias, const std::array<std::vector<double>, num_surrounding_blocks_> &sur_bias) const;
		double SumBias(double &max_bias, seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type fragment_length, double general_bias, const Vect<double> &gc_bias, const std::array<std::vector<double>, num_surrounding_blocks_> &sur_bias) const;
		void GetFragmentSites(std::vector<FragmentSite> &sites, seqan::Size<decltype(reference_sequences_)>::Type ref_seq_id, seqan::Size<seqan::Dna5String>::Type fragment_length, uint32_t start, uint32_t end) const;
#endif //SWIG

		bool ReadFasta(const char *fasta_file);
		void ReplaceN( uint64_t seed );
		bool hasN() const;
		bool WriteFasta(const char *fasta_file) const;

#ifndef SWIG // This part is not needed for the python plotting
		inline uint16_t NumAlleles() const{ return num_alleles_; }
		inline bool VariantsLoaded() const{ return variants_.size(); }
		inline bool VariantPositionsLoaded() const{ return variant_positions_.size(); }
		inline bool VariantsLoadedForSequence(uint32_t ref_seq_id) const{ return read_variation_for_num_sequences_ >= ref_seq_id; }
		inline bool VariantPositionsLoadedForSequence(uint32_t ref_seq_id) const{ return read_variation_for_num_sequences_ >= ref_seq_id; }
		inline bool VariantsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= NumberSequences(); }
		inline bool VariantPositionsCompletelyLoaded() const{ return read_variation_for_num_sequences_ >= NumberSequences(); }
		inline const std::vector<Reference::Variant> &Variants(uint32_t ref_seq_id) const{ return variants_.at(ref_seq_id); }
		inline const std::vector<uint32_t> &VariantPositions(uint32_t ref_seq_id) const{ return variant_positions_.at(ref_seq_id); }

		bool PrepareVariantFile(const std::string &var_file);
		inline bool ReadVariants(uint32_t end_ref_seq_id){ return ReadVariants(end_ref_seq_id, false); };
		inline bool ReadVariantPositions(uint32_t end_ref_seq_id){ return ReadVariants(end_ref_seq_id, true); };
		bool ReadFirstVariants();
		bool ReadFirstVariantPositions();
		void ClearVariants(uint32_t end_ref_seq_id);
		void ClearVariantPositions(uint32_t end_ref_seq_id);
		void ClearAllVariants();
		void ClearAllVariantPositions();
#endif //SWIG
	};

	struct FragmentSite{
		uint16_t gc_;
		std::array<int32_t, Reference::num_surrounding_blocks_> start_surrounding_;
		std::array<int32_t, Reference::num_surrounding_blocks_> end_surrounding_;
		uint32_t count_forward_;
		uint32_t count_reverse_;
		double bias_;

		FragmentSite(uint16_t gc, const std::array<int32_t, Reference::num_surrounding_blocks_> &start_sur, const std::array<int32_t, Reference::num_surrounding_blocks_> &end_sur):
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
