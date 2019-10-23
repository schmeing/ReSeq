#include "CoverageStats.h"
using reseq::CoverageStats;

#include <algorithm>
using std::min;
//include <array>
using std::array;
//include <mutex>
using std::mutex;
using std::lock_guard;
//include <utility>
using std::pair;

//include <seqan/bam_io.h>
using seqan::Dna;
using seqan::Dna5;
using seqan::hasFlagRC;
using seqan::hasFlagLast;

//include "utilities.hpp"
using reseq::utilities::ConstDna5StringReverseComplement;
using reseq::utilities::ConstIupacStringReverseComplement;
using reseq::utilities::ReversedConstCharString;
using reseq::utilities::ReversedConstCigarString;
using reseq::utilities::Complement;
using reseq::utilities::at;
using reseq::utilities::Divide;
using reseq::utilities::GetDominantLastX;
using reseq::utilities::IsN;
using reseq::utilities::Percent;
using reseq::utilities::SafePercent;
using reseq::utilities::SetDominantLastX;
using reseq::utilities::SetToMin;
using reseq::utilities::TransformDistanceToStartOfErrorRegion;

bool CoverageStats::IsVariantPosition( intVariantId &cur_var, uintRefSeqId ref_seq_id, uintSeqLen pos, const Reference &reference, bool reverse_direction ){
	if( reference.VariantPositionsLoaded() ){
		if(reverse_direction){
			while(0 <= cur_var && reference.VariantPositions(ref_seq_id).at(cur_var) > pos){
				--cur_var;
			}

			return 0 <= cur_var && reference.VariantPositions(ref_seq_id).at(cur_var) == pos;
		}
		else{
			while(cur_var < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(cur_var) < pos){
				++cur_var;
			}

			return cur_var < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(cur_var) == pos;
		}
	}

	return false;
}

void CoverageStats::PrepareVariantPositionCheck( intVariantId &cur_var, uintRefSeqId ref_seq_id, uintSeqLen pos, const Reference &reference, bool reverse_direction ){
	if( reference.VariantPositionsLoaded() ){
		if(reverse_direction){
			--cur_var; // To guarantee we are not at the end
			while(cur_var+1 < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(cur_var+1) <= pos){
				++cur_var;
			}
		}
		else{
			while(cur_var && reference.VariantPositions(ref_seq_id).at(cur_var-1) >= pos){
				--cur_var;
			}
		}
	}
}

void CoverageStats::EvalRead( FullRecord *record, CoverageStats::CoverageBlock *coverage_block, const Reference &reference, QualityStats &qualities, ErrorStats &errors, uintQual phred_quality_offset){
	if(0 == record->num_pointers_to_element_){
		printErr << "Accessing removed FullRecord: seqid(" << record->record_.rID << ") pos(" << record->record_.beginPos << ") name(" << record->record_.qName << ")" << std::endl;
		throw std::out_of_range( "Accessing removed FullRecord" );
		return;
	}
	// 1: Fragments that did not made the filters have to_ref_pos_=0 and even if they are from valid fragments are reads added to the coverage blocks by maximum read length and not by its actual read length, so it might happen that they do not really reach the given block
	// 2: Because of potential soft-clipping reads are added to coverage blocks before beginPos, which they might not belong to
	if(record->to_ref_pos_ > coverage_block->start_pos_ && record->from_ref_pos_ < coverage_block->start_pos_+coverage_block->coverage_.size()){
		ConstIupacStringReverseComplement reversed_seq(record->record_.seq);
		ReversedConstCharString reversed_qual(record->record_.qual);
		ReversedConstCigarString rev_cigar(record->record_.cigar);

		uintReadLen read_pos(0), ref_pos(0), ref_pos_start(0), ref_pos_end;
		uintSeqLen coverage_pos;
		Dna ref_base;
		Dna5 base;
		uintQual qual(0), last_qual(1);
		auto cur_var = coverage_block->first_variant_id_;

		// Set ref_pos_start to the first ref_pos(0 is read start) that belongs to this block
		// Set ref_pos_end to the ref_pos(0 is read start) that is one past the last in this block and or in the read(adapters are excluded from the read already), whatever is first
		if( hasFlagRC(record->record_) ){
			coverage_pos = record->to_ref_pos_-coverage_block->start_pos_-1;

			if(record->to_ref_pos_ > coverage_block->start_pos_+coverage_block->coverage_.size()){
				ref_pos_start = record->to_ref_pos_-coverage_block->start_pos_-coverage_block->coverage_.size();
			}

			if(record->from_ref_pos_ < coverage_block->start_pos_){
				ref_pos_end = record->to_ref_pos_-coverage_block->start_pos_;
			}
			else{
				ref_pos_end = record->to_ref_pos_-record->from_ref_pos_;
			}

			PrepareVariantPositionCheck( cur_var, record->record_.rID, record->to_ref_pos_-1, reference, true );
		}
		else{
			coverage_pos = record->from_ref_pos_-coverage_block->start_pos_; // If we have an overrun because it goes below 0 we won't use it until it's back at 0 anyways (due to ref_pos_start), so no problem

			if(record->from_ref_pos_ < coverage_block->start_pos_){
				ref_pos_start = coverage_block->start_pos_-record->from_ref_pos_;
			}

			if(record->to_ref_pos_ > coverage_block->start_pos_+coverage_block->coverage_.size()){
				ref_pos_end = coverage_block->start_pos_+coverage_block->coverage_.size()-record->from_ref_pos_;
			}
			else{
				ref_pos_end = record->to_ref_pos_-record->from_ref_pos_;
			}

			PrepareVariantPositionCheck( cur_var, record->record_.rID, record->from_ref_pos_, reference, false );
		}

		uintReadLen num_errors = 0;
		uintPercent error_rate;
		uintSeqLen full_ref_pos;
		for( const auto &cigar_element : (hasFlagRC(record->record_)?rev_cigar:record->record_.cigar) ){
			switch( cigar_element.operation ){
			case 'M':
			case '=':
			case 'X':
			case 'S': // Treat soft-clipping as match, so that bwa and bowtie2 behave the same
				for(auto i=cigar_element.count; i--; ){
					if(ref_pos < ref_pos_end){
						// We are currently not comparing the adapter to the reference
						if( hasFlagRC(record->record_) ){
							base = at(reversed_seq, read_pos);
							qual = at(reversed_qual, read_pos)-phred_quality_offset;
							full_ref_pos = record->to_ref_pos_-1-ref_pos;
							ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), full_ref_pos));
						}
						else{
							base = at(record->record_.seq, read_pos);
							qual = at(record->record_.qual, read_pos)-phred_quality_offset;
							full_ref_pos = record->from_ref_pos_+ref_pos;
							ref_base = at(reference.ReferenceSequence(record->record_.rID), full_ref_pos);
						}

						if(ref_pos >= ref_pos_start){
							if( IsN(ref_base) || IsVariantPosition(cur_var, record->record_.rID, full_ref_pos, reference, hasFlagRC(record->record_))){
								++(record->not_considered_ref_bases_);
							}
							else{
								error_rate = coverage_block->coverage_.at(coverage_pos).error_rate_.at(hasFlagRC(record->record_));
								qualities.AddRefBase(hasFlagLast(record->record_), ref_base, coverage_block->coverage_.at(coverage_pos).dom_error_.at(hasFlagRC(record->record_)), record->tile_id_, qual, error_rate, last_qual, record->sequence_quality_, read_pos);
								errors.AddBase(hasFlagLast(record->record_), ref_base, coverage_block->coverage_.at(coverage_pos).dom_error_.at(hasFlagRC(record->record_)), record->tile_id_, base, qual, read_pos, num_errors, error_rate);
								record->error_rate_sum_ += error_rate;
							}
						}

						if(ref_base != base && !IsN(ref_base) && !IsVariantPosition(cur_var, record->record_.rID, full_ref_pos, reference, hasFlagRC(record->record_))){
							++num_errors;
						}

						last_qual = qual;

						++ref_pos;
						++read_pos;

						if( hasFlagRC(record->record_) ){
							--coverage_pos;
						}
						else{
							++coverage_pos;
						}
					}
				}
				break;
			case 'N':
			case 'D':
				for(auto i=cigar_element.count; i--; ){
					if(ref_pos < ref_pos_end){
						// We are currently not comparing the adapter to the reference
						if( hasFlagRC(record->record_) ){
							full_ref_pos = record->to_ref_pos_-1-ref_pos;
							ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), full_ref_pos));
						}
						else{
							full_ref_pos = record->from_ref_pos_+ref_pos;
							ref_base = at(reference.ReferenceSequence(record->record_.rID), full_ref_pos);
						}

						if(ref_pos >= ref_pos_start){
							if(IsVariantPosition(cur_var, record->record_.rID, full_ref_pos, reference, hasFlagRC(record->record_))){
								++(record->not_considered_ref_bases_);
							}
							else{
								record->error_rate_sum_ += coverage_block->coverage_.at(coverage_pos).error_rate_.at(hasFlagRC(record->record_));
							}
						}

						if( !IsN(ref_base) && !IsVariantPosition(cur_var, record->record_.rID, full_ref_pos, reference, hasFlagRC(record->record_))){
							++num_errors;
						}

						++ref_pos;

						if( hasFlagRC(record->record_) ){
							--coverage_pos;
						}
						else{
							++coverage_pos;
						}
					}
				}
				break;
			case 'I':
				if(ref_pos < ref_pos_end){
					if( hasFlagRC(record->record_) ){
						last_qual = at(reversed_qual, read_pos+cigar_element.count-1)-phred_quality_offset; // Directly fill last_qual as qual is not needed here
						full_ref_pos = record->to_ref_pos_-1-ref_pos;
						ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), full_ref_pos));
					}
					else{
						last_qual = at(record->record_.qual, read_pos+cigar_element.count-1)-phred_quality_offset; // Directly fill last_qual as qual is not needed here
						full_ref_pos = record->from_ref_pos_+ref_pos;
						ref_base = at(reference.ReferenceSequence(record->record_.rID), full_ref_pos);
					}

					if( !IsN(ref_base) && !IsVariantPosition(cur_var, record->record_.rID, full_ref_pos, reference, hasFlagRC(record->record_))){
						num_errors += cigar_element.count;
					}

					read_pos += cigar_element.count;
				}
				break;
			}
		}
	}

	DeleteRecord( record, qualities );
}

void CoverageStats::UpdateZeroCoverageRegion(intFragCountShift zero_coverage_length){
	coverage_[ 0 ] += zero_coverage_length;
	for( int strand=2; strand--; ){
		coverage_stranded_.at(strand)[0] += zero_coverage_length;
	}
}

void CoverageStats::UpdateCoverageAtSinglePosition(CoveragePosition &coverage, Dna5 ref_base, bool is_variant_position){
	Dna5 reverse_base = Complement::Dna5(ref_base);

	// Stored error rate for the simulation is the rate of the dominant error, real error rate is calculated further below for plotting stats
	uintCovCount errors(0);
	Dna5 dom_error(4);
	uintCovCount coverage_forward = 0;
	for( int i=5; i--; ){
		coverage_forward += coverage.coverage_forward_.at(i);
		if(i != ref_base){
			if(coverage.coverage_forward_.at(i) > errors){
				errors = coverage.coverage_forward_.at(i);
				dom_error = i;
			}
		}
	}
	if(coverage_forward){
		tmp_errors_forward_.emplace_back(dom_error, Percent(errors,coverage_forward), coverage_threshold_ <= coverage_forward);
		coverage.dom_error_.at(0) = dom_error;
		coverage.error_rate_.at(0) = Percent(errors,coverage_forward);
	}
	else{
		tmp_errors_forward_.emplace_back(4, 0, false);
		coverage.dom_error_.at(0) = 4;
		coverage.error_rate_.at(0) = 0;
	}

	errors = 0;
	dom_error = 4;
	uintCovCount coverage_reverse = 0;
	for( int i=5; i--; ){
		coverage_reverse += coverage.coverage_reverse_.at(i);
		if(i != reverse_base){
			if(coverage.coverage_reverse_.at(i) > errors){
				errors = coverage.coverage_reverse_.at(i);
				dom_error = i;
			}
		}
	}
	if(coverage_reverse){
		tmp_errors_reverse_.emplace_back(dom_error, Percent(errors,coverage_reverse), coverage_threshold_ <= coverage_reverse);
		coverage.dom_error_.at(1) = dom_error;
		coverage.error_rate_.at(1) = Percent(errors,coverage_reverse);
	}
	else{
		tmp_errors_reverse_.emplace_back(4, 0, false);
		coverage.dom_error_.at(1) = 4;
		coverage.error_rate_.at(1) = 0;
	}

	auto errors_forward = coverage_forward - coverage.coverage_forward_.at(ref_base);
	auto errors_reverse = coverage_reverse - coverage.coverage_reverse_.at(reverse_base);

	auto cov_sum = coverage_forward + coverage_reverse;
	auto error_sum = errors_forward + errors_reverse;

	// Account for coverage and error_coverage at that position
	++coverage_[ cov_sum ];
	++coverage_stranded_.at(0)[ coverage_forward ];
	++coverage_stranded_.at(1)[ coverage_reverse ];
	if(cov_sum){
		++coverage_stranded_percent_.at(0)[ Percent(coverage_forward, cov_sum) ];
		++coverage_stranded_percent_.at(1)[ Percent(coverage_reverse, cov_sum) ];

		bool errors_valid = !is_variant_position && !IsN(ref_base);
		if(errors_valid){
			++error_coverage_[ error_sum ];
			++error_coverage_percent_[ Percent(error_sum, cov_sum) ];
			if(coverage_forward && coverage_reverse){
				++error_coverage_percent_stranded_[ Percent(errors_forward, coverage_forward) ][ Percent(errors_reverse, coverage_reverse) ];
			}
		}

		if( 10 <= cov_sum ){
			++coverage_stranded_percent_min_cov_10_.at(0)[ Percent(coverage_forward, cov_sum) ];
			++coverage_stranded_percent_min_cov_10_.at(1)[ Percent(coverage_reverse, cov_sum) ];
			if(errors_valid){
				++error_coverage_percent_min_cov_10_[ Percent(error_sum, cov_sum) ];
			}

			if( 20 <= cov_sum ){
				++coverage_stranded_percent_min_cov_20_.at(0)[ Percent(coverage_forward, cov_sum) ];
				++coverage_stranded_percent_min_cov_20_.at(1)[ Percent(coverage_reverse, cov_sum) ];
				if(errors_valid){
					++error_coverage_percent_min_cov_20_[ Percent(error_sum, cov_sum) ];

					if(10 <=coverage_forward && 10 <= coverage_reverse){
						++error_coverage_percent_stranded_min_strand_cov_10_[ Percent(errors_forward, coverage_forward) ][ Percent(errors_reverse, coverage_reverse) ];
						if(20 <=coverage_forward && 20 <= coverage_reverse){
							++error_coverage_percent_stranded_min_strand_cov_20_[ Percent(errors_forward, coverage_forward) ][ Percent(errors_reverse, coverage_reverse) ];
						}
					}
				}
			}
		}
	}
}

void CoverageStats::UpdateDistances(uintSeqLen &distance_to_start_of_error_region, uintPercent error_rate) const{
	if( distance_to_start_of_error_region ){
		if( ++distance_to_start_of_error_region >= reset_distance_ ){
			distance_to_start_of_error_region = 0;
		}
	}
	else{
		if(error_rate >= kErrorRateThreshold){
			distance_to_start_of_error_region = 1;
		}
	}
}

CoverageStats::CoverageBlock *CoverageStats::CreateBlock(uintRefSeqId seq_id, uintSeqLen start_pos){
	CoverageBlock *new_block;
	if( reuse_mutex_.try_lock() ){
		if( reusable_blocks_.size() ){
			new_block = reusable_blocks_.back();
			reusable_blocks_.pop_back();
			reuse_mutex_.unlock();

			new_block->sequence_id_ = seq_id;
			new_block->start_pos_ = start_pos;
			new_block->previous_block_ = last_block_;
			new_block->next_block_ = NULL;
			new_block->coverage_.clear();
			new_block->reads_.clear();
			new_block->first_variant_id_ = 0;
		}
		else{
			reuse_mutex_.unlock();
			new_block = new CoverageBlock(seq_id, start_pos, last_block_);
			new_block->reads_.reserve( (*last_block_).reads_.capacity() );
		}
	}
	else{
		new_block = new CoverageBlock(seq_id, start_pos, last_block_);
		new_block->reads_.reserve( (*last_block_).reads_.capacity() );
	}

	new_block->coverage_.resize( kBlockSize );
	(*last_block_).next_block_ = new_block;

	return new_block;
}

void CoverageStats::CountBlock(CoverageBlock *block, const Reference &reference){
	auto cur_var = block->first_variant_id_;
	for( uintSeqLen pos = 0; pos < block->coverage_.size(); ++pos ){
		UpdateCoverageAtSinglePosition( block->coverage_.at(pos), at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos), IsVariantPosition( cur_var, block->sequence_id_, block->start_pos_ + pos, reference, false ) );
	}

	// Enter error_rates_ and dominant_errors_ in forward direction
	Dna5 ref_base, prev_ref_base(4), dom_ref_base5(4);
	array<uintReadLen,5> seq_content_reference_last5 = {0,0,0,0,0};
	if(block->start_pos_){
		prev_ref_base = at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ - 1);
		GetDominantLastX( dom_ref_base5, seq_content_reference_last5, 5, reference.ReferenceSequence(block->sequence_id_), block->start_pos_ );
	}
	uintSeqLen gc_bases = min(block->start_pos_, gc_range_);
	uintSeqLen n_count(0);
	uintSeqLen gc = reference.GCContentAbsolut( n_count, block->sequence_id_, block->start_pos_-gc_bases, block->start_pos_);
	uintPercent gc_percent;

	cur_var = block->first_variant_id_;

	for( uintSeqLen pos = 0; pos < tmp_errors_forward_.size(); ++pos ){
		ref_base = at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos);
		gc_percent = SafePercent(gc,gc_bases-n_count);

		// Enter values
		bool valid = tmp_errors_forward_.at(pos).coverage_sufficient_ && !IsN(ref_base) && !IsVariantPosition( cur_var, block->sequence_id_, block->start_pos_ + pos, reference, false );
		if(valid){
			++dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_forward_)][tmp_errors_forward_.at(pos).base_];
			++dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[gc_percent][tmp_errors_forward_.at(pos).base_];
			++gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_forward_)][gc_percent];

			++error_rates_by_distance_.at(ref_base).at(tmp_errors_forward_.at(pos).base_)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_forward_)][tmp_errors_forward_.at(pos).rate_];
			++error_rates_by_gc_.at(ref_base).at(tmp_errors_forward_.at(pos).base_)[gc_percent][tmp_errors_forward_.at(pos).rate_];
			++gc_by_distance_er_.at(ref_base).at(tmp_errors_forward_.at(pos).base_)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_forward_)][gc_percent];
		}

		// Set values for next iteration
		prev_ref_base = ref_base;
		SetDominantLastX( dom_ref_base5, seq_content_reference_last5, ref_base, 5, reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos);
		UpdateDistances(distance_to_start_of_error_region_forward_, (valid?tmp_errors_forward_.at(pos).rate_:0));
		auto new_pos = block->start_pos_ + pos;
		if(gc_bases<gc_range_){
			++gc_bases;
			if( reference.GC( block->sequence_id_, new_pos ) ){
				++gc;
			}
			else if( reference.N( block->sequence_id_, new_pos ) ){
				++n_count;
			}
		}
		else{
			reference.UpdateGC( gc, n_count, block->sequence_id_, new_pos - gc_bases, new_pos);
		}
	}

	// Remove already handled elements from tmp_errors_forward_
	tmp_errors_forward_.clear(); // In forward direction we can always handle all

	// Prepare distance_to_start_of_error_region_forward_ for next block
	if(block->next_block_){
		if( (*(block->next_block_)).sequence_id_ == block->sequence_id_){
			// Reference sequence changes between this block and the next or we are at the end of the file: Error region will be reseted
			distance_to_start_of_error_region_forward_ = 0;
		}
		else{
			distance_to_start_of_error_region_forward_ += (*(block->next_block_)).start_pos_ - block->start_pos_ - block->coverage_.size(); // Take the uncovered region between blocks into account
			if( reset_distance_ < distance_to_start_of_error_region_forward_ ){
				// Error region ends
				distance_to_start_of_error_region_forward_ = 0;
			}
		}
	}

	// Enter error_rates_ and dominant_errors_ in reverse direction
	uintSeqLen last_pos = tmp_errors_reverse_.size();
	if(block->next_block_ && (*(block->next_block_)).sequence_id_ == block->sequence_id_){
		// Everything that is potentially in an error region starting in the next block cannot be handled yet
		last_pos -= reset_distance_+1;
	}
	else{
		// Reference sequence changes between this block and the next or we are at the end of the file: We can handle all positions
		--last_pos; // Shift it from end pos to last valid pos
	}

	if(last_pos){
		ConstDna5StringReverseComplement reversed_seq(reference.ReferenceSequence(block->sequence_id_));
		uintSeqLen ref_pos = reference.SequenceLength(block->sequence_id_) + tmp_errors_reverse_.size() - block->start_pos_ - block->coverage_.size() - last_pos - 1; // tmp_errors_reverse_.size() - block->coverage_.size() is the shift to the left of block->start_pos_ due to the non-handleable reverse errors taken from the previous block
		seq_content_reference_last5.fill(0);
		if(ref_pos){
			prev_ref_base = at(reversed_seq, ref_pos - 1);
			GetDominantLastX( dom_ref_base5, seq_content_reference_last5, 5, reversed_seq, ref_pos );
		}
		else{
			prev_ref_base = 4;
			dom_ref_base5 = 4;
		}
		gc_bases = min(ref_pos, gc_range_);
		n_count = 0;
		gc = reference.GCContentAbsolut( n_count, block->sequence_id_, reference.SequenceLength(block->sequence_id_)-ref_pos, reference.SequenceLength(block->sequence_id_)-ref_pos+gc_bases);

		uintSeqLen distance_to_start_of_error_region_reverse(0);

		PrepareVariantPositionCheck( cur_var, block->sequence_id_, block->start_pos_ + block->coverage_.size() - tmp_errors_reverse_.size() + last_pos, reference, true );

		for( uintSeqLen pos = last_pos+1; pos--; ){
			ref_base = at(reversed_seq, ref_pos);
			gc_percent = SafePercent(gc,gc_bases-n_count);

			// Enter values
			bool valid = tmp_errors_reverse_.at(pos).coverage_sufficient_ && !IsN(ref_base) && !IsVariantPosition( cur_var, block->sequence_id_, block->start_pos_ + block->coverage_.size() - tmp_errors_reverse_.size() + pos, reference, true );
			if(valid){
				++dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_reverse)][tmp_errors_reverse_.at(pos).base_];
				++dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[gc_percent][tmp_errors_reverse_.at(pos).base_];
				++gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base5)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_reverse)][gc_percent];

				++error_rates_by_distance_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_reverse)][tmp_errors_reverse_.at(pos).rate_];
				++error_rates_by_gc_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_)[gc_percent][tmp_errors_reverse_.at(pos).rate_];
				++gc_by_distance_er_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_)[TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_reverse)][gc_percent];
			}

			// Set values for next iteration
			prev_ref_base = ref_base;
			SetDominantLastX( dom_ref_base5, seq_content_reference_last5, ref_base, 5, reversed_seq, ref_pos);
			UpdateDistances(distance_to_start_of_error_region_reverse, (valid?tmp_errors_reverse_.at(pos).rate_:0));
			auto new_pos = reference.SequenceLength(block->sequence_id_)-ref_pos-1;
			if(gc_bases<gc_range_){
				++gc_bases;
				if( reference.GC( block->sequence_id_, new_pos ) ){
					++gc;
				}
				else if( reference.N( block->sequence_id_, new_pos ) ){
					++n_count;
				}
			}
			else{
				reference.UpdateGC( gc, n_count, block->sequence_id_, new_pos + gc_bases, new_pos);
			}
			ref_pos++;
		}

		// Remove already handled elements from tmp_errors_reverse_
		tmp_errors_reverse_.erase(tmp_errors_reverse_.begin(), tmp_errors_reverse_.begin()+(last_pos+1));
	}
}

CoverageStats::CoverageBlock *CoverageStats::RemoveBlock(CoverageBlock *block){
	reusable_blocks_.push_back(block);
	return block->next_block_;
}

CoverageStats::CoverageBlock *CoverageStats::FindBlock(uintRefSeqId ref_seq_id, uintSeqLen ref_pos){
	CoverageBlock *block = last_block_;
	while(block->sequence_id_ > ref_seq_id){
		block = block->previous_block_;
	}
	while(block->start_pos_ > ref_pos){
		block = block->previous_block_;
	}

	return block;
}

bool CoverageStats::EnsureSpace(uintRefSeqId ref_seq_id, uintSeqLen start_pos, uintSeqLen end_pos, FullRecord *record, Reference &reference, uintReadLen minimum_sequence_length){
	if(0 == record->num_pointers_to_element_){
		printErr << "Accessing removed FullRecord: seqid(" << record->record_.rID << ") pos(" << record->record_.beginPos << ") name(" << record->record_.qName << ")" << std::endl;
		throw std::out_of_range( "Accessing removed FullRecord" );
		return false;
	}

	// Make sure the variation for the given reference sequence is already loaded
	if( reference.VariantPositionsLoaded() && !reference.VariantPositionsLoadedForSequence(ref_seq_id) ){
		lock_guard<mutex> lock(variant_loading_mutex_);
		if(!reference.ReadVariantPositions(ref_seq_id+1)){ // End is specified, so to get the current one +1 is needed
			return false;
		}
	}

	SetToMin(end_pos, reference.SequenceLength(ref_seq_id));
	if( last_block_ ){
		// As reads are ordered by position, we can jump over not covered regions
		if((*last_block_).sequence_id_ < ref_seq_id){
			// Account for the not covered region at the end of the last reference sequence
			zero_coverage_region_ += reference.SequenceLength((*last_block_).sequence_id_) - (*last_block_).start_pos_ - (*last_block_).coverage_.size();

			// Account for the not covered reference sequences in between
			for( auto seq_id = ref_seq_id; --seq_id > (*last_block_).sequence_id_ ; ){
				if(reference.SequenceLength(seq_id) > minimum_sequence_length ){
					zero_coverage_region_ += reference.SequenceLength(seq_id);
				}
				else{
					zero_coverage_region_ += 2*Reference::kMinDistToRefSeqEnds; // This is what we remove later from every sequence, so make sure we get to zero and not below
				}
			}

			// Account for the not covered region in the current reference sequence before the first position
			zero_coverage_region_ += start_pos;

			// Create first block of next reference sequence
			last_block_ = CreateBlock(ref_seq_id, start_pos);
		}
		else if((*last_block_).start_pos_+kBlockSize < start_pos){
			// Account for the empty coverage region
			zero_coverage_region_ += start_pos - (*last_block_).start_pos_ - kBlockSize;

			// Create first block after the empty coverage region
			auto first_var = (*last_block_).first_variant_id_;
			last_block_ = CreateBlock(ref_seq_id, start_pos);

			// Update first_variant_id_ for the new block
			if( reference.VariantPositionsLoaded() ){
				while( first_var < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(first_var) < (*last_block_).start_pos_ ){
					++first_var;
				}
				(*last_block_).first_variant_id_ = first_var;
			}
		}
		else if((*last_block_).start_pos_ > start_pos){
			// Add read to blocks before the last
			auto block = (*last_block_).previous_block_;
			while( block->start_pos_ > start_pos ){
				block->reads_.push_back(record);
				++(record->num_pointers_to_element_);
				block = block->previous_block_;
			};
			block->reads_.push_back(record);
			++(record->num_pointers_to_element_);
		}
	}
	else{
		// Initialize first block
		// Account for the not covered reference sequences before the first position
		for( auto seq_id = ref_seq_id; seq_id-- ; ){
			if(reference.SequenceLength(seq_id) > minimum_sequence_length ){
				zero_coverage_region_ += reference.SequenceLength(seq_id);
			}
			else{
				zero_coverage_region_ += 2*Reference::kMinDistToRefSeqEnds; // This is what we remove later from every sequence, so make sure we get to zero and not below
			}
		}

		// Account for the not covered region in the current reference sequence before the first position
		zero_coverage_region_ += start_pos;

		// Create new block
		first_block_ = new CoverageBlock(ref_seq_id, start_pos, NULL);
		first_block_->coverage_.resize( kBlockSize );
		first_block_->reads_.reserve( 2*kBlockSize );

		last_block_ = first_block_;
	}

	// Create blocks until end_pos is reached
	while((*last_block_).start_pos_+kBlockSize < end_pos){
		(*last_block_).reads_.push_back(record); // Add record to reads_ of all blocks it is part of (except last one, which is added after the loop)
		++(record->num_pointers_to_element_);
		last_block_ = CreateBlock((*last_block_).sequence_id_, (*last_block_).start_pos_+kBlockSize);
	}
	(*last_block_).reads_.push_back(record); // Add to reads_ for the last block it reaches into
	++(record->num_pointers_to_element_);

	// Ensure that last block is not reaching over the end of the reference sequence
	if( reference.SequenceLength((*last_block_).sequence_id_)-(*last_block_).start_pos_ < kBlockSize){
		(*last_block_).coverage_.resize( reference.SequenceLength((*last_block_).sequence_id_)-(*last_block_).start_pos_ );
	}

	return true;
}

void CoverageStats::AddFragment( uintRefSeqId ref_seq_id, uintSeqLen ref_pos, CoverageBlock *&block ){
	if(block){
		while( block->start_pos_+kBlockSize <= ref_pos || block->sequence_id_ < ref_seq_id ){
			block = block->next_block_;
		}
	}
	else{
		block = FindBlock(ref_seq_id, ref_pos);
	}
	++(block->unprocessed_fragments_);
}

void CoverageStats::RemoveFragment( uintRefSeqId ref_seq_id, uintSeqLen ref_pos, CoverageBlock *&block, uintFragCount &processed_fragments ){
	if(block){
		if( block->start_pos_+kBlockSize <= ref_pos || block->sequence_id_ < ref_seq_id ){
			auto old_block = block;
			do{
				block = block->next_block_;
			} while( block->start_pos_+kBlockSize <= ref_pos || block->sequence_id_ < ref_seq_id );

			// Subtract processed_fragments from old_block only after finding the new block, as it potentially allows other threads to remove old_block
			old_block->unprocessed_fragments_ -= processed_fragments;
			processed_fragments = 0;
		}
		else if(block->start_pos_ > ref_pos){ // Here we don't need to compare to sequence_id as pairs with reads on different contigs are not added to unprocessed_fragments_
			// Here we don't need to care if we set unprocessed_fragments_ to zero before finding the new block, as it is before this block and must have at least one unprocessed_fragments_
			block->unprocessed_fragments_ -= processed_fragments;
			processed_fragments = 0;

			do{
				block = block->previous_block_;
			} while( block->start_pos_ > ref_pos );
		}
	}
	else{
		block = FindBlock(ref_seq_id, ref_pos);
	}

	++processed_fragments; // Store processed fragments instead of subtracting them directly so that last fragments from this batch can be removed after the CleanUp step, so at least one block remains during CleanUp
}

reseq::uintRefSeqId CoverageStats::CleanUp(uintSeqLen &still_needed_position, Reference &reference, QualityStats &qualities, ErrorStats &errors, uintQual phred_quality_offset, CoverageBlock *cov_block, uintFragCount &processed_fragments){
	CoverageBlock *until_block;
	uintRefSeqId still_needed_reference_sequence(0);
	still_needed_position = 0;
	{
		lock_guard<mutex> lock(clean_up_mutex_);

		until_block = first_block_;
		while(!until_block->unprocessed_fragments_){ // The unprocessed fragments from the current read haven't been removed yet, so end of list cannot be reached
			CountBlock(until_block, reference);
			until_block = until_block->next_block_;
		}
		until_block = until_block->previous_block_;

		if(until_block){
			first_block_ = until_block->next_block_;
			still_needed_reference_sequence = first_block_->sequence_id_;
			still_needed_position = first_block_->start_pos_;
			first_block_->previous_block_ = NULL;
			until_block->next_block_ = NULL;
		}
	}

	// Subtract processed_fragments from cov_block after until_block is determined, as it potentially allows to remove cov_block
	cov_block->unprocessed_fragments_ -= processed_fragments;
	processed_fragments = 0;

	if(until_block){
		while(until_block->previous_block_){
			for( auto rec : until_block->reads_){
				EvalRead(rec, until_block, reference, qualities, errors, phred_quality_offset);
			}
			until_block = until_block->previous_block_;
		}

		for( auto rec : until_block->reads_){
			EvalRead(rec, until_block, reference, qualities, errors, phred_quality_offset);
		}

		lock_guard<mutex> lock(reuse_mutex_);
		do{
			until_block = RemoveBlock(until_block);
		} while(until_block);
	}

	reference.ClearVariants(still_needed_reference_sequence);

	return still_needed_reference_sequence;
}

bool CoverageStats::PreLoadVariants(Reference &reference){
	if(reference.VariantPositionsLoaded() && !reference.VariantPositionsCompletelyLoaded() && !reference.VariantPositionsLoadedForSequence((*last_block_).sequence_id_+2)){
		if( variant_loading_mutex_.try_lock() ){
			if( !reference.ReadVariantPositions((*last_block_).sequence_id_+2) ){
				variant_loading_mutex_.unlock();
				return false;
			}

			variant_loading_mutex_.unlock();
		}
	}

	return true;
}

bool CoverageStats::Finalize(const Reference &reference, QualityStats &qualities, ErrorStats &errors, uintQual phred_quality_offset, uintReadLen minimum_sequence_length){
	// Remove all reusable_blocks_ as they are not needed anymore
	for( auto block : reusable_blocks_ ){
		delete block;
	}
	reusable_blocks_.clear();

	if(first_block_){
		if( last_block_ ){
			// Update not covered regions after last block
			zero_coverage_region_ += reference.SequenceLength((*last_block_).sequence_id_)-(*last_block_).start_pos_-(*last_block_).coverage_.size();

			for( auto seq_id = (*last_block_).sequence_id_+1; seq_id < reference.NumberSequences(); ++seq_id ){
				if( reference.SequenceLength(seq_id) > minimum_sequence_length ){
					zero_coverage_region_ += reference.SequenceLength(seq_id);
				}
				else{
					zero_coverage_region_ += 2*Reference::kMinDistToRefSeqEnds; // This is what we remove later from every sequence, so make sure we get to zero and not below
				}
			}

			// Update coverage of remaining blocks and delete them
			while(first_block_->next_block_){
				CountBlock(first_block_, reference);
				for( auto rec : first_block_->reads_){
					EvalRead(rec, first_block_, reference, qualities, errors, phred_quality_offset);
				}

				first_block_ = first_block_->next_block_;
				delete first_block_->previous_block_;
			}

			CountBlock(first_block_, reference);
			for( auto rec : first_block_->reads_){
				EvalRead(rec, first_block_, reference, qualities, errors, phred_quality_offset);
			}
			delete first_block_;
			first_block_ = NULL;
			last_block_ = NULL;
		}
		else{
			printErr << "first_block_ set, but last_block_ not. Something dodgy is going on here" << std::endl;
			return false;
		}
	}
	else{
		if( last_block_ ){
			printErr << "last_block_ set, but first_block_ not. Something dodgy is going on here" << std::endl;
			return false;
		}
		else{
			// There have been reads in the file (or Finalize wouldn't have been called), but none made the criteria to add coverage
			SetAllZero(reference);
		}
	}

	UpdateZeroCoverageRegion(zero_coverage_region_-2*Reference::kMinDistToRefSeqEnds*reference.NumberSequences()); // Subtract excluded region due to minimum distance to reference sequence ends for accepting fragments
	tmp_errors_forward_.shrink_to_fit();
	tmp_errors_reverse_.shrink_to_fit();

	return true;
}

void CoverageStats::Shrink(){
	coverage_.Shrink();
	for( int strand=2; strand--; ){
		coverage_stranded_.at(strand).Shrink();
		coverage_stranded_percent_.at(strand).Shrink();
		coverage_stranded_percent_min_cov_10_.at(strand).Shrink();
		coverage_stranded_percent_min_cov_20_.at(strand).Shrink();
	}
	error_coverage_.Shrink();
	error_coverage_percent_.Shrink();
	error_coverage_percent_min_cov_10_.Shrink();
	error_coverage_percent_min_cov_20_.Shrink();
	ShrinkVect(error_coverage_percent_stranded_);
	ShrinkVect(error_coverage_percent_stranded_min_strand_cov_10_);
	ShrinkVect(error_coverage_percent_stranded_min_strand_cov_20_);
}

void CoverageStats::PreparePlotting(){
	error_rates_by_distance_sum_.Clear();
	error_rates_by_gc_sum_.Clear();
	for( int ref_base=4; ref_base--; ){
		for( int dom_error=5; dom_error--; ){
			error_rates_by_distance_sum_ += error_rates_by_distance_.at(ref_base).at(dom_error);
			error_rates_by_gc_sum_ += error_rates_by_gc_.at(ref_base).at(dom_error);
		}
	}
}
