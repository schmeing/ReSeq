#include "Surrounding.h"
using reseq::Surrounding;
using reseq::SurroundingBias;

#include <algorithm>
using std::min;
//include <array>
using std::array;
//include<vector>
using std::vector;

#include "reportingUtils.hpp"
//include <seqan/sequence.h>
using seqan::Dna;
using seqan::DnaString;

//include "utilities.hpp"
using reseq::uintSurBlockId;
using reseq::uintSurPos;
using reseq::utilities::at;

void Surrounding::ChangeSurroundingBase(uintSurPos pos, Dna new_base){
	auto bit_in_block = 2*(kRange - (pos%kRange) - 1);
	// ~ = bit-wise not
	sur_.at(pos/kRange) = ( sur_.at(pos/kRange) & ~(3 << bit_in_block) ) + ( static_cast<intType>(new_base) << bit_in_block );
}

void Surrounding::DeleteSurroundingBaseShiftingOnRightSide(uintSurPos pos, Dna new_end_base){
	intType new_base = new_end_base;
	auto del_block = pos/kRange;

	// Add base at end and shift bases from end block to block with deletion
	for(auto block=kNumBlocks; --block > del_block; ){
		sur_.at(block) = (sur_.at(block) << 2) + new_base;
		new_base = sur_.at(block) / Size();
		sur_.at(block) = sur_.at(block) % Size();
	}

	// Remove deleted base and add shifted base
	auto bit_in_block = 2*(kRange - (pos%kRange) - 1);
	new_base += ( sur_.at(del_block) % (1 << bit_in_block) ) << 2;
	sur_.at(del_block) = ( sur_.at(del_block) >> (bit_in_block+2) << (bit_in_block+2) ) + new_base;
}

void Surrounding::DeleteSurroundingBaseShiftingOnLeftSide(uintSurPos pos, Dna new_end_base){
	intType new_base = new_end_base;
	auto del_block = pos/kRange;

	// Add base in front and shift bases from first block to block with deletion
	for(uintSurBlockId block=0; block < del_block; ++block){
		sur_.at(block) += new_base*Size();
		new_base = sur_.at(block)%4;
		sur_.at(block) = sur_.at(block) >> 2;
	}

	// Add shifted base and remove deleted base
	auto bit_in_block = 2*(kRange - (pos%kRange) - 1);
	sur_.at(del_block) += new_base*Size();
	sur_.at(del_block) = ( sur_.at(del_block) >> (bit_in_block+2) << (bit_in_block) ) + sur_.at(del_block) % (1 << bit_in_block);
}

void Surrounding::InsertSurroundingBasesShiftingOnRightSide(uintSurPos pos, DnaString new_bases){
	auto block = pos/kRange;
	uintSurPos bases_to_insert = min(length(new_bases), static_cast<size_t>(Length()-pos));

	// Shift right from current block to make room for new bases
	auto shift_blocks = bases_to_insert/kRange; // For every kRange in bases_to_insert we can shift whole blocks
	auto shift_bases = bases_to_insert%kRange;
	for(auto cur_block=kNumBlocks-shift_blocks; cur_block-- > block+1; ){
		// First shift only the bases and ignore the blocks at the end that we completely remove
		sur_.at(cur_block) >>= 2*shift_bases; // Remove bases not longer needed
		sur_.at(cur_block) += sur_.at(cur_block-1) % (1 << 2*shift_bases) * (1 << 2*(kRange-shift_bases)); // Add new bases from previous block
	}

	uintSurPos inv_pos_in_block = kRange - (pos%kRange);
	intType tmp_sur = sur_.at(block) % (1 << 2*inv_pos_in_block); // Store everything that needs to be shifted out of the start block
	sur_.at(block) >>= 2*inv_pos_in_block; // And remove it from the block
	tmp_sur >>= 2*shift_bases; // Remove what was already shifted

	if(0 < shift_blocks && kNumBlocks > block+shift_blocks){
		for(auto cur_block=kNumBlocks; cur_block-- > block+shift_blocks+1; ){
			// Now shift the blocks
			sur_.at(cur_block) = sur_.at(cur_block-shift_blocks);
		}

		// Shift everything that has to be shifted right from starting block and was not already shifted above (Happens if the bases need to go into an empty block freed by the block shift)
		sur_.at(block+shift_blocks) = tmp_sur;
	}

	// Insert new bases in starting block
	uintSurPos ins_pos=0;
	uintSurPos bases_to_insert_into_this_block = min(bases_to_insert, inv_pos_in_block);
	bases_to_insert -= bases_to_insert_into_this_block;
	for(; bases_to_insert_into_this_block--; ){
		sur_.at(block) <<= 2;
		sur_.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos++));
	}

	if(0 == shift_blocks && inv_pos_in_block > shift_bases){
		// If the complete insertion happens in the starting block add the stuff back that was shifted out too much
		sur_.at(block) <<= 2*(inv_pos_in_block-shift_bases);
		sur_.at(block) += tmp_sur;
	}
	else{
		while(0 < bases_to_insert){
			bases_to_insert_into_this_block = min(bases_to_insert, (uintSurPos)kRange);
			bases_to_insert -= bases_to_insert_into_this_block;

			// Store the part of the block that is not overwritten by insertions
			tmp_sur = sur_.at(++block) % (1 << 2*(kRange-bases_to_insert_into_this_block));

			// Insert new bases
			sur_.at(block) = 0;
			for(auto i=bases_to_insert_into_this_block; i--; ){
				sur_.at(block) <<= 2;
				sur_.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos++));
			}

			// Restore the bases that should not be overwritten
			sur_.at(block) <<= 2*(kRange-bases_to_insert_into_this_block);
			sur_.at(block) += tmp_sur;
		}
	}
}

void Surrounding::InsertSurroundingBasesShiftingOnLeftSide(uintSurPos pos, DnaString new_bases){
	auto block = pos/kRange;
	uintSurPos bases_to_insert = min(length(new_bases), static_cast<size_t>(pos+1));

	// Shift left from current block to make room for new bases
	auto shift_blocks = bases_to_insert/kRange; // For every kRange in bases_to_insert we can shift whole blocks
	auto shift_bases = bases_to_insert%kRange;
	for(uintSurBlockId cur_block=shift_blocks; cur_block < block; ++cur_block ){
		// First shift only the bases and ignore the blocks at the beginning that we completely remove
		sur_.at(cur_block) %= 1 << 2*(kRange-shift_bases);// Remove bases not longer needed
		sur_.at(cur_block) <<= 2*shift_bases; // Shift remaining bases (do this after the removal to avoid overflow)
		sur_.at(cur_block) += sur_.at(cur_block+1) >> 2*(kRange-shift_bases); // Add new bases from next block
	}

	uintSurPos pos_in_block = pos%kRange+1; // 1-bases position in block
	intType tmp_sur = sur_.at(block) % (1 << 2*(kRange-pos_in_block)); // Store everything that needs to stay in the start block
	if(pos_in_block > shift_bases){
		// Keep the stuff that is after the shift on the left before the insertion
		sur_.at(block) >>= 2*(kRange-pos_in_block);
		sur_.at(block) %= 1 << 2*(pos_in_block-shift_bases); // Remove what was already shifted
	}
	else{
		sur_.at(block) = 0; // If everything was already shifted, we only need what's right of the insertion and that is in tmp_sur
	}

	if(0 < shift_blocks){
		if(block >= shift_blocks){
			for(uintSurBlockId cur_block=0; cur_block+shift_blocks < block; ++cur_block){
				// Now shift the blocks
				sur_.at(cur_block) = sur_.at(cur_block+shift_blocks);
			}

			// Shift everything that has to be shifted left from starting block and was not already shifted above (Happens if the bases need to go into an empty block freed by the block shift)
			sur_.at(block-shift_blocks) = sur_.at(block) << 2*(kRange-(pos_in_block-shift_bases));
		}
		sur_.at(block) = 0; // The not-shifted stuff to keep was moved to the proper shifted block, so now we can clean the start block surrounding
	}

	// Insert new bases in starting block
	uintSurPos bases_to_insert_into_this_block = min(bases_to_insert, pos_in_block);
	uintSurPos ins_pos_to = length(new_bases); // In case we do not insert all bases from new_bases, we need to keep track of this separately to bases_to_insert
	for(uintSurPos ins_pos=ins_pos_to-bases_to_insert_into_this_block; ins_pos<ins_pos_to; ++ins_pos){
		sur_.at(block) <<= 2;
		sur_.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos));
	}
	bases_to_insert -= bases_to_insert_into_this_block;
	ins_pos_to -= bases_to_insert_into_this_block;

	// Add back what was right of insertion
	sur_.at(block) <<= 2*(kRange-pos_in_block);
	sur_.at(block) += tmp_sur;

	while(0 < bases_to_insert){
		bases_to_insert_into_this_block = min(bases_to_insert, (uintSurPos)kRange);

		// Leave only the part of the block that should not be overwritten by insertions
		sur_.at(--block) >>= 2*bases_to_insert_into_this_block;

		// Insert new bases
		for(uintSurPos ins_pos=ins_pos_to-bases_to_insert_into_this_block; ins_pos<ins_pos_to; ++ins_pos){
			sur_.at(block) <<= 2;
			sur_.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos));
		}
		bases_to_insert -= bases_to_insert_into_this_block;
		ins_pos_to -= bases_to_insert_into_this_block;
	}
}

void SurroundingBias::CombinePositions(const array<double, 4*Surrounding::Length()> &separated){
	for( auto &block : bias_ ){
		block.clear();
		block.resize(Surrounding::Size(), 0.0);
	}

	vector<uintBaseCall> bases;
	uintSeqLen pos;
	for( uintSurBlockId block=0; block < Surrounding::kNumBlocks; ++block ){
		bases.clear();
		bases.resize(Surrounding::kRange+1, 0); // We need a buffer to catch the final increase in the loop so +1 here
		for(Surrounding::intType sur=0; sur < Surrounding::Size(); ++sur){
			// Bias for a surrounding is the sum/product of all biases from the given bases at its positions
			for(pos=0; pos < Surrounding::kRange; ++pos){
				bias_.at(block).at(sur) += separated.at(bases.at(Surrounding::kRange-1-pos) + pos*4 + block*Surrounding::kRange*4);
			}

			// Update the current bases we are at: If a position exceeds valid bases, set it back to A and increase next position
			// bases[0] is for pos=surrounding_range_-1
			pos=0;
			while( ++(bases.at(pos)) > 3 ){
				bases.at(pos++) = 0;
			}
		}
	}
}

void SurroundingBias::SeparatePositions(array<double, 4*Surrounding::Length()> &separated) const{
	separated.fill(0.0);

	std::vector<uintBaseCall> bases;
	uintSeqLen pos;
	for( uintSurBlockId block=0; block < Surrounding::kNumBlocks; ++block ){
		bases.clear();
		bases.resize(Surrounding::kRange+1, 0); // We need a buffer to catch the final increase in the loop so +1 here
		for(Surrounding::intType sur=0; sur < Surrounding::Size(); ++sur){
			// Add surrounding bias to the corresponding base at each position
			for(pos=0; pos < Surrounding::kRange; ++pos){
				separated.at(bases.at(Surrounding::kRange-1-pos) + pos*4 + block*Surrounding::kRange*4) += bias_.at(block).at(sur);
			}

			// Update the current bases we are at: If a position exceeds valid bases, set it back to A and increase next position
			// bases[0] is for pos=surrounding_range_-1
			pos=0;
			while( ++(bases.at(pos)) > 3 ){
				bases.at(pos++) = 0;
			}
		}
	}

	// Normalize separated values
	for(auto sur_pos = 0; sur_pos < Surrounding::Length(); ++sur_pos){
		double sur_sum(0.0);
		for(auto sur = 4*sur_pos; sur < 4*sur_pos+4; ++sur){
			sur_sum += separated.at(sur);
		}
		sur_sum /= 4;
		for(auto sur = 4*sur_pos; sur < 4*sur_pos+4; ++sur){
			separated.at(sur) -= sur_sum;
			separated.at(sur) /= Surrounding::Size()/4; // The base of the current position is defined so divide by 4 to get number of averaged values
		}
	}
}

