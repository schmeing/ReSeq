#include "FragmentDuplicationStats.h"
using reseq::FragmentDuplicationStats;

//include <algorithm>
using std::max;
using std::sort;
#include<iterator>
using std::distance;
#include <list>
using std::list;
//include <cmath>
using std::sqrt;
//include <set>
using std::set;
#include <utility>
using std::pair;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::Dna5;

//include "utilities.hpp"
using reseq::utilities::SetToMax;
using reseq::utilities::VectorAtomic;

inline void FragmentDuplicationStats::SumUpDuplicationNumber(){
	duplication_number_.Clear();
	duplication_number_.SetOffset(1);

	for( auto &frag_dups_by_length : duplication_number_by_gc_insert_length_ ){
		for( auto &frag_dups : frag_dups_by_length ){
			if( frag_dups.size() ){
				for( auto dup=frag_dups.to(); dup-- > max(static_cast<size_t>(1),frag_dups.from()); ){ // Ignore the zeros
					duplication_number_[dup] += frag_dups.at(dup);
				}
			}
		}
	}
	duplication_number_.Shrink();
}

void FragmentDuplicationStats::AddDuplicates( vector<uintSeqLen> &fragment_positions, uintRefSeqId ref_seq_id, uintSeqLen insert_length, uintSeqLen base_pos, const Reference &reference ){
	sort(fragment_positions.begin(), fragment_positions.end());

	uintSeqLen cur_pos(0), last_pos(0);
	auto corrected_pos = max(reference.StartExclusion(ref_seq_id), base_pos);
	auto next_exclusion_region = reference.NextExclusionRegion(ref_seq_id, corrected_pos);
	uintDupCount count = 0;
	for( auto pos : fragment_positions ){
		if(pos == cur_pos){
			++count;
		}
		else{
			if(count){
				corrected_pos += cur_pos/2-last_pos;
				reference.CorrectPositionAddingExcludedRegions(corrected_pos, next_exclusion_region, ref_seq_id);
				AddDuplication(insert_length, reference.GCContent(ref_seq_id, corrected_pos, corrected_pos+insert_length), count);
				last_pos = cur_pos/2;
			}

			cur_pos = pos;
			count = 1;
		}
	}
	// Insert last position
	corrected_pos += cur_pos/2-last_pos;
	reference.CorrectPositionAddingExcludedRegions(corrected_pos, next_exclusion_region, ref_seq_id);
	AddDuplication(insert_length, reference.GCContent(ref_seq_id, corrected_pos, corrected_pos+insert_length), count);
}

void FragmentDuplicationStats::FinalizeDuplicationVector(const vector<vector<VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc){
	uintFragCount site_sum;
	for(uintSeqLen insert_length=tmp_duplication_number_by_gc_insert_length_.size(); insert_length-- > 1; ){
		for(uintPercent gc=tmp_duplication_number_by_gc_insert_length_.at(insert_length).size(); gc-- ; ){
			site_sum = 0;
			for(uintDupCount dup=tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).size(); dup-- > 1; ){
				duplication_number_by_gc_insert_length_[gc][insert_length][dup] += tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(dup);
				site_sum += tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(dup);
			}
			if(insert_length < site_count_by_insert_length_gc.size()){
				duplication_number_by_gc_insert_length_[gc][insert_length][0] = 2*site_count_by_insert_length_gc.at(insert_length).at(gc) - site_sum;
			}
		}
	}

	ShrinkVect(duplication_number_by_gc_insert_length_);
	tmp_duplication_number_by_gc_insert_length_.clear();
	tmp_duplication_number_by_gc_insert_length_.shrink_to_fit();
}

void FragmentDuplicationStats::CalculateDispersionPlot(){
	uintFragCount count;
	double mean, var;

	uintNumFits needed_size=0;
	for(auto &dup_by_len : duplication_number_by_gc_insert_length_ ){
		needed_size += dup_by_len.size();
	}

	mean_list_.reserve(needed_size);
	dispersion_list_.reserve(needed_size);
	for(auto gc=duplication_number_by_gc_insert_length_.to(); gc-- > duplication_number_by_gc_insert_length_.from(); ){
		for(auto len=duplication_number_by_gc_insert_length_.at(gc).to(); len-- > duplication_number_by_gc_insert_length_.at(gc).from(); ){
			if(2 < duplication_number_by_gc_insert_length_.at(gc).at(len).to()){ // For distributions that only have counts at 0 or 1 dispersion is always smaller than mean
				// Calculate mean u
				count = 0; // Count of non-zero sites
				mean = 0.0;
				for(auto dup=duplication_number_by_gc_insert_length_.at(gc).at(len).to(); dup-- > duplication_number_by_gc_insert_length_.at(gc).at(len).from(); ){
					count += duplication_number_by_gc_insert_length_.at(gc).at(len).at(dup);
					mean += dup * duplication_number_by_gc_insert_length_.at(gc).at(len).at(dup);
				}
				mean /= count;

				// Calculate variance o^2
				var = 0.0;
				for(auto dup=duplication_number_by_gc_insert_length_.at(gc).at(len).to(); dup-- > duplication_number_by_gc_insert_length_.at(gc).at(len).from(); ){
					var += (mean-dup) * (mean-dup) * duplication_number_by_gc_insert_length_.at(gc).at(len).at(dup);
				}
				var /= count;

				// Dispersion: r = u^2/(o^2-u)
				if( var > mean ){
					mean_list_.push_back(mean);
					dispersion_list_.push_back(mean*mean/(var-mean));
				}
			}
		}
	}

	mean_list_.shrink_to_fit();
	dispersion_list_.shrink_to_fit();
}

void FragmentDuplicationStats::PreparePlotting(){
	SumUpDuplicationNumber();
	CalculateDispersionPlot();
}

void FragmentDuplicationStats::PrepareTesting(){
	SumUpDuplicationNumber();
	CalculateDispersionPlot();
}
