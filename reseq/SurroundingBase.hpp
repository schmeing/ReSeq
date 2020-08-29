#ifndef SURROUNDINGBASE_HPP
#define SURROUNDINGBASE_HPP

#include <array>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "utilities.hpp"

namespace reseq{
	template<uintSurBlockId N,uintSurPos R, int16_t S, typename I> class SurroundingBase{
	protected:
		static const uintSurBlockId kNumBlocks = N; // Number of surrounding blocks
		static const uintSurPos kRange = R; // Length of a single surrounding block
		static constexpr size_t Size(){
			return 1 << 2*static_cast<size_t>(kRange);
		}

	public:
		typedef I intType;
		static const int16_t kStartPos = S; // Position where the surrounding should start relative to the first base in the fragment (so 10 is 10 bases before the fragment) [0 <= kStartPos < kNumBlocks*kRange]
		static constexpr uintSurPos Length(){
			return kNumBlocks * kRange;
		}

		static uintSeqLen StartPos(uintSeqLen pos, uintSeqLen seq_len){
			return pos+seq_len - kStartPos;
		}
		static uintSeqLen EndPos(uintSeqLen pos, uintSeqLen seq_len){
			return StartPos(pos, seq_len) + Length();
		}

		static uintSeqLen StartPosReverse(uintSeqLen pos, uintSeqLen seq_len){
			return pos+seq_len + kStartPos;
		}
		static uintSeqLen EndPosReverse(uintSeqLen pos, uintSeqLen seq_len){
			return StartPosReverse(pos, seq_len) - Length();
		}

	protected:
		std::array<intType, kNumBlocks> sur_;

		uintSeqLen BlockStartPos(uintSeqLen pos, uintSeqLen seq_len, uintSurBlockId block) const{
			return StartPos(pos, seq_len) + block*kRange;
		}

		uintSeqLen BlockLastPos(uintSeqLen pos, uintSeqLen seq_len, uintSurBlockId block) const{
			return StartPos(pos, seq_len) + (block+1)*kRange - 1;
		}

		uintSeqLen BlockStartPosReverse(uintSeqLen pos, uintSeqLen seq_len, uintSurBlockId block) const{
			return StartPosReverse(pos, seq_len) - block*kRange;
		}

		uintSeqLen BlockLastPosReverse(uintSeqLen pos, uintSeqLen seq_len, uintSurBlockId block) const{
			return StartPosReverse(pos, seq_len) - (block+1)*kRange + 1;
		}

		template<typename T> intType GetBlock(const T &sequence, uintSeqLen start_pos) const{
			intType sur(0);

			for(auto pos = start_pos; pos < start_pos+kRange; ++pos){
				sur <<= 2;
				sur += static_cast<uintBaseCall>(utilities::at(sequence, pos%length(sequence)));
			}

			return sur;
		}

		template<typename T> void Set(const T &sequence, uintSeqLen pos){
			for(auto block=0; block < kNumBlocks; ++block){
				// Add sequence length to avoid negative positions if block reaches over lower sequence ends so it wraps around and start at the end
				// The higher sequence end is covered by the called function
				sur_.at(block) = GetBlock(sequence, BlockStartPos(pos, length(sequence), block));
			}
		}

		bool BlockHasN(intType &invalid_sur, const seqan::Dna5String &sequence, uintSeqLen start_pos){
			for( auto cur_pos = start_pos+kRange; cur_pos-- > start_pos; ){
				if(utilities::IsN(sequence, cur_pos%length(sequence))){
					invalid_sur = static_cast<intType>(start_pos) - cur_pos - 1; // Negative counts of base increments needed to not have an N anymore
					return true; // Last N is only interesting one
				}
			}

			return false;
		}

		void UpdateOnRightSide(seqan::Dna new_base){
			for(auto block=0; block < kNumBlocks-1; ++block){
				sur_.at(block) <<= 2;
				sur_.at(block) %= Size();
				sur_.at(block) += sur_.at(block+1) >> 2*(kRange-1);
			}

			sur_.at(kNumBlocks-1) <<= 2;
			sur_.at(kNumBlocks-1) %= Size();
			sur_.at(kNumBlocks-1) += static_cast<uintBaseCall>(new_base);
		}

		void UpdateOnLeftSide(seqan::Dna new_base){
			for(auto block=kNumBlocks; block-- > 1; ){
				sur_.at(block) >>= 2;
				sur_.at(block) += sur_.at(block-1)%4 << 2*(kRange-1);
			}

			sur_.at(0) >>= 2;
			sur_.at(0) += static_cast<uintBaseCall>(new_base) << 2*(kRange-1);
		}

		void UpdateWithN(const seqan::Dna5String &sequence, uintSeqLen new_pos, bool forward){
			uintBaseCall new_base;
			for(auto block=kNumBlocks; block--; ){
				if(forward){
					new_base = utilities::at( sequence, BlockLastPos(new_pos, length(sequence), block)%length(sequence) );
				}
				else{
					new_base = utilities::Complement::Dna5( utilities::at( sequence, BlockStartPosReverse(new_pos, length(sequence), block)%length(sequence) ) );
				}

				if(utilities::IsN(new_base)){
					sur_.at(block) = -kRange; // Wait kRange positions until calculating surrounding again
				}
				else if(0 > sur_.at(block)){
					if(-1 == sur_.at(block)){
						// Calculate current surrounding
						if(forward){
							sur_.at(block) = GetBlock(sequence, BlockStartPos(new_pos, length(sequence), block));
						}
						else{
							sur_.at(block) = GetBlock(utilities::ConstDna5StringReverseComplement(sequence), BlockStartPos(length(sequence)-new_pos-1, length(sequence), block));
						}
					}
					else{
						++(sur_.at(block));
					}
				}
				else{
					// No N: Update surrounding normally
					if(forward){
						sur_.at(block) = (sur_.at(block) << 2)%Size() + new_base;
					}
					else{
						sur_.at(block) = (sur_.at(block) + new_base*Size()) >> 2;
					}
				}
			}
		}

	public:
		SurroundingBase(){
			sur_.fill(-1);
		}

		SurroundingBase(const SurroundingBase &origin){
			sur_ = origin.sur_;
		}

#ifndef SWIG // To avoid the warning that the operator= will be ignored
		SurroundingBase& operator=(const SurroundingBase& x){
			for(auto block=kNumBlocks; block--; ){
				sur_.at(block) = x.sur_.at(block);
			}
			return *this;
		}

		bool operator==(const SurroundingBase& right) const{
			for(auto block=kNumBlocks; block--; ){
				if( sur_.at(block) != right.sur_.at(block) ){
					return false;
				}
			}
			return true;
		}
#endif //SWIG

		inline uintBaseCall BaseAt(uintSurPos pos) const{
			return (sur_.at(pos/kRange) >> 2*(kRange - pos%kRange - 1)) % 4;
		}

		bool Valid() const{
			for(auto block=kNumBlocks; block--; ){
				if(0 > sur_.at(block)){
					return false;
				}
			}

			return true;
		}

		void Forward(const seqan::Dna5String &sequence, uintSeqLen pos){
			Set(sequence, pos);
		}

		void Reverse(const seqan::Dna5String &sequence, uintSeqLen pos){
			Set(utilities::ConstDna5StringReverseComplement(sequence), length(sequence)-pos-1);
		}

		void ForwardWithN(const seqan::Dna5String &sequence, uintSeqLen pos){
			for(auto block=0; block < kNumBlocks; ++block){
				if( !BlockHasN(sur_.at(block), sequence, BlockStartPos(pos, length(sequence), block)) ){
					sur_.at(block) = GetBlock(sequence, BlockStartPos(pos, length(sequence), block));
				}
			}
		}

		void ReverseWithN(const seqan::Dna5String &sequence, uintSeqLen pos){
			for(auto block=0; block < kNumBlocks; ++block){
				if( !BlockHasN(sur_.at(block), sequence, BlockLastPosReverse(pos, length(sequence), block)) ){
					sur_.at(block) = GetBlock(utilities::ConstDna5StringReverseComplement(sequence), BlockStartPos(length(sequence)-pos-1, length(sequence), block));
				}
			}
		}

		inline void UpdateForward(const seqan::Dna5String &sequence, uintSeqLen new_pos){
			UpdateOnRightSide( utilities::at(sequence, (EndPos(new_pos, length(sequence))-1)%length(sequence)) );
		}

		inline void UpdateReverse(const seqan::Dna5String &sequence, uintSeqLen new_pos){
			UpdateOnLeftSide( utilities::Complement::Dna5( utilities::at(sequence, (StartPosReverse(new_pos, length(sequence)))%length(sequence)) ) );
		}

		inline void RollBackReverse(const seqan::Dna5String &sequence, uintSeqLen new_pos){
			UpdateOnRightSide( utilities::Complement::Dna5( utilities::at(sequence, (EndPosReverse(new_pos, length(sequence))+1)%length(sequence)) ) );
		}

		inline void UpdateForwardWithN(const seqan::Dna5String &sequence, uintSeqLen new_pos){
			UpdateWithN(sequence, new_pos, true);
		}

		inline void UpdateReverseWithN(const seqan::Dna5String &sequence, uintSeqLen new_pos){
			UpdateWithN(sequence, new_pos, false);
		}
	};
}

#endif // SURROUNDINGBASE_HPP
