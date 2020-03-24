#include "FragmentDuplicationStats.h"
using reseq::FragmentDuplicationStats;

//include <algorithm>
using std::sort;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

void FragmentDuplicationStats::AddDuplicates( vector<uintSeqLen> &fragment_positions ){
	sort(fragment_positions.begin(), fragment_positions.end());

	uintSeqLen cur_pos(0);
	uintDupCount count = 0;
	for( auto pos : fragment_positions ){
		if(pos == cur_pos){
			++count;
		}
		else{
			AddDuplication(count);

			cur_pos = pos;
			count = 1;
		}
	}
	// Insert last position
	AddDuplication(count);
}

void FragmentDuplicationStats::FinalizeDuplicationVector(){
	// Multiply counts by duplication number to go from #FragmentSites to #Fragments
	for(uintDupCount dup=0; dup<tmp_duplication_number_.size(); ++dup){
		tmp_duplication_number_.at(dup) = tmp_duplication_number_.at(dup) * dup;
	}

	duplication_number_.Acquire(tmp_duplication_number_);
}
