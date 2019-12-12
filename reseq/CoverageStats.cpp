#include "CoverageStats.h"
using reseq::CoverageStats;

//include <algorithm>
using std::min;
using std::max;
//include <array>
using std::array;
//include <cmath>
using std::exp;
using std::pow;
using std::round;
#include <limits>
using std::numeric_limits;
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
using reseq::utilities::DominantBase;
using reseq::utilities::IsN;
using reseq::utilities::Percent;
using reseq::utilities::SafePercent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;
using reseq::utilities::TransformDistanceToStartOfErrorRegion;

void CoverageStats::EvalRead( FullRecord *record, CoverageStats::CoverageBlock *coverage_block, const Reference &reference, QualityStats &qualities, ErrorStats &errors, uintQual phred_quality_offset){
	// Fragments that did not made the filters have to_ref_pos_=0
	if(record->to_ref_pos_){
		ConstIupacStringReverseComplement reversed_seq(record->record_.seq);
		ReversedConstCharString reversed_qual(record->record_.qual);
		ReversedConstCigarString rev_cigar(record->record_.cigar);

		intSeqShift coverage_pos; // Can be before the start of the block, so negative values need to be allowed
		uintSeqLen ref_pos;
		if( hasFlagRC(record->record_) ){
			coverage_pos = record->to_ref_pos_-coverage_block->start_pos_-1;
			ref_pos = record->to_ref_pos_-1;
		}
		else{
			coverage_pos = record->from_ref_pos_-coverage_block->start_pos_;
			ref_pos = record->from_ref_pos_;
		}

		uintReadLen read_pos(0), not_considered_ref_bases(0);
		Dna ref_base;
		Dna5 base, dom_error;
		uintQual qual(0), last_qual(1);
		uintReadLen num_errors = 0;
		uintPercent error_rate;
		uintReadLenCalc error_rate_sum(0);
		for( const auto &cigar_element : (hasFlagRC(record->record_)?rev_cigar:record->record_.cigar) ){
			switch( cigar_element.operation ){
			case 'M':
			case '=':
			case 'X':
			case 'S': // Treat soft-clipping as match, so that bwa and bowtie2 behave the same
				for(auto i=cigar_element.count; i--; ){
					if(record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_){
						// We are currently not comparing the adapter to the reference
						if( hasFlagRC(record->record_) ){
							base = at(reversed_seq, read_pos);
							qual = at(reversed_qual, read_pos)-phred_quality_offset;
							ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), ref_pos));
						}
						else{
							base = at(record->record_.seq, read_pos);
							qual = at(record->record_.qual, read_pos)-phred_quality_offset;
							ref_base = at(reference.ReferenceSequence(record->record_.rID), ref_pos);
						}

						if( !CoveragePosValid(coverage_pos, coverage_block) ){
							++(not_considered_ref_bases);
						}
						else{
							if(0 > coverage_pos){
								error_rate = coverage_block->previous_coverage_.at(-1*coverage_pos-1).error_rate_.at(hasFlagRC(record->record_));
								dom_error = coverage_block->previous_coverage_.at(-1*coverage_pos-1).dom_error_.at(hasFlagRC(record->record_));
							}
							else{
								error_rate = coverage_block->coverage_.at(coverage_pos).error_rate_.at(hasFlagRC(record->record_));
								dom_error = coverage_block->coverage_.at(coverage_pos).dom_error_.at(hasFlagRC(record->record_));
							}
							qualities.AddRefBase(hasFlagLast(record->record_), ref_base, dom_error, record->tile_id_, qual, error_rate, last_qual, record->sequence_quality_, read_pos);
							errors.AddBase(hasFlagLast(record->record_), ref_base, dom_error, record->tile_id_, base, qual, read_pos, num_errors, error_rate);

							error_rate_sum += error_rate;

							if( ref_base != base ){
								++num_errors;
							}
						}

						last_qual = qual;

						++read_pos;

						if( hasFlagRC(record->record_) ){
							--coverage_pos;
							--ref_pos;
						}
						else{
							++coverage_pos;
							++ref_pos;
						}
					}
				}
				break;
			case 'N':
			case 'D':
				for(auto i=cigar_element.count; i--; ){
					if(record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_){
						// We are currently not comparing the adapter to the reference
						if( hasFlagRC(record->record_) ){
							ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), ref_pos));
						}
						else{
							ref_base = at(reference.ReferenceSequence(record->record_.rID), ref_pos);
						}

						if( !CoveragePosValid(coverage_pos, coverage_block) ){
							++(not_considered_ref_bases);
						}
						else{
							if(0 > coverage_pos){
								error_rate_sum += coverage_block->previous_coverage_.at(-1*coverage_pos-1).error_rate_.at(hasFlagRC(record->record_));
							}
							else{
								error_rate_sum += coverage_block->coverage_.at(coverage_pos).error_rate_.at(hasFlagRC(record->record_));
							}
							++num_errors;
						}

						if( hasFlagRC(record->record_) ){
							--coverage_pos;
							--ref_pos;
						}
						else{
							++coverage_pos;
							++ref_pos;
						}
					}
				}
				break;
			case 'I':
				if(record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_){
					if( hasFlagRC(record->record_) ){
						last_qual = at(reversed_qual, read_pos+cigar_element.count-1)-phred_quality_offset; // Directly fill last_qual as qual is not needed here
						ref_base = Complement::Dna5(at(reference.ReferenceSequence(record->record_.rID), ref_pos));
					}
					else{
						last_qual = at(record->record_.qual, read_pos+cigar_element.count-1)-phred_quality_offset; // Directly fill last_qual as qual is not needed here
						ref_base = at(reference.ReferenceSequence(record->record_.rID), ref_pos);
					}

					if( CoveragePosValid(coverage_pos, coverage_block) ){
						num_errors += cigar_element.count;
					}

					read_pos += cigar_element.count;
				}
				break;
			}
		}

		if(not_considered_ref_bases < record->to_ref_pos_-record->from_ref_pos_){ // At least one base was considered
			qualities.AddRefRead(hasFlagLast(record->record_), record->tile_id_, record->reference_gc_, record->sequence_quality_, utilities::Divide(error_rate_sum, static_cast<uintReadLenCalc>(record->to_ref_pos_-record->from_ref_pos_-not_considered_ref_bases)), record->fragment_length_);
		}
	}

	delete record;
}

inline void CoverageStats::ApplyZeroCoverageRegion(){
	tmp_coverage_.at( 0 ) += zero_coverage_region_;
	tmp_coverage_.at( 0 ) -= excluded_bases_;
	for( int strand=2; strand--; ){
		tmp_coverage_stranded_.at(strand).at(0) += zero_coverage_region_;
		tmp_coverage_stranded_.at(strand).at(0) -= excluded_bases_;
	}
}

reseq::uintCovCount CoverageStats::GetNormalErrorCount(double error_rate, uintCovCount coverage){
	// binom: (1-p)^(n-k) * p^k * prod(j=1..k)[(n+1-j)/j]
	double pmf = pow(1-error_rate, coverage);
	uintNucCount normal_error_count(1);
	double probability = 1-pmf; // Probability to have at least normal_error_count errors at position

	while(probability > kRejectingBinomialP){
		pmf *= error_rate/(1-error_rate) * static_cast<double>(coverage+1-normal_error_count)/normal_error_count;
		probability -= pmf;
		++normal_error_count;
	}

	return normal_error_count;
}

pair<double, double> CoverageStats::GetMaxNonSystematicError(CoverageBlock *block, const Reference &reference){
	uintNucCount total_bases(0), correct_bases(0);
	uintCovCount bases_forward, bases_reverse, min_cov(numeric_limits<uintCovCount>::max()), max_cov(0);
	for( uintSeqLen pos = 0; pos < block->coverage_.size(); ++pos ){
		if( block->coverage_.at(pos).valid_ ){
			bases_forward = 0;
			bases_reverse = 0;
			for( int i=5; i--; ){
				bases_forward += block->coverage_.at(pos).coverage_forward_.at(i);
				bases_reverse += block->coverage_.at(pos).coverage_reverse_.at(i);
			}

			total_bases += bases_forward + bases_reverse;
			auto ref_base = at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos);
			correct_bases += block->coverage_.at(pos).coverage_forward_.at(ref_base) + block->coverage_.at(pos).coverage_reverse_.at(Complement::Dna5(ref_base));

			SetToMin(min_cov, bases_forward);
			SetToMin(min_cov, bases_reverse);
			SetToMax(max_cov, bases_forward);
			SetToMax(max_cov, bases_reverse);
		}
	}

	if(3 > max_cov || correct_bases == total_bases){
		return {numeric_limits<uintCovCount>::max(), 0.0}; // Nothing is systematic for super low coverage or if we only have errors
	}

	// Linear approximation for errors over coverage (using the two points at 0.25 and 0.75 of the coverage interval)
	double error_rate = static_cast<double>(total_bases-correct_bases)/total_bases;
	++block_error_rate_[static_cast<uintPercent>(round(error_rate*100))];
	error_rate /= 3.0; // Systematic is only for one base

	uintCovCount lower_coverage = max(static_cast<uintCovCount>(2), min_cov+(max_cov-min_cov)/4); // GetNormalErrorCount calculation doesn't make sense for coverage lower than 2
	auto lower_errors = GetNormalErrorCount(error_rate, lower_coverage);
	++block_non_systematic_error_rate_[Percent(min(lower_errors, lower_coverage), lower_coverage)];

	uintCovCount upper_coverage = max(min_cov+1, max_cov-(max_cov-min_cov)/4); // Prevent min_cov == max_cov
	auto upper_errors = GetNormalErrorCount(error_rate, upper_coverage);
	++block_non_systematic_error_rate_[Percent(min(upper_errors, upper_coverage), upper_coverage)];

	pair<double, double> non_systematic_error;
	non_systematic_error.second = static_cast<double>(upper_errors-lower_errors) / (upper_coverage-lower_coverage); // slope
	non_systematic_error.first = static_cast<double>(lower_errors) - non_systematic_error.second*lower_coverage; // y-interrcept

	return non_systematic_error;
}

void CoverageStats::UpdateCoverageAtSinglePosition(CoveragePosition &coverage, Dna5 ref_base, pair<double, double> non_systematic_error){
	Dna5 reverse_base = Complement::Dna5(ref_base);

	// Stored error rate for the simulation is the rate of the dominant error, real error rate is calculated further below for plotting stats
	uintCovCount errors(1); // Enforce at least 2 errors to call it a systematic error
	uintPercent error_rate;
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
	if(coverage_forward && coverage.valid_){
		if(errors <= GetMaxNonSystematicError(non_systematic_error, coverage_forward)){
			error_rate = 0;
			dom_error = 4;
		}
		else{
			error_rate = Percent(errors,coverage_forward);
		}
		tmp_errors_forward_.emplace_back(dom_error, error_rate, coverage_threshold_ <= coverage_forward);
		coverage.dom_error_.at(0) = dom_error;
		coverage.error_rate_.at(0) = error_rate;
	}
	else{
		tmp_errors_forward_.emplace_back(4, 0, false);
		coverage.dom_error_.at(0) = 4;
		coverage.error_rate_.at(0) = 0;
	}

	errors = 1; // Enforce at least 2 errors to call it a systematic error
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
	if(coverage_reverse && coverage.valid_){
		if(errors <= GetMaxNonSystematicError(non_systematic_error, coverage_reverse)){
			error_rate = 0;
			dom_error = 4;
		}
		else{
			error_rate = Percent(errors,coverage_reverse);
		}
		tmp_errors_reverse_.emplace_back(dom_error, error_rate, coverage_threshold_ <= coverage_reverse);
		coverage.dom_error_.at(1) = dom_error;
		coverage.error_rate_.at(1) = error_rate;
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
	++tmp_coverage_.at( min(cov_sum, kMaxCoverage) );
	++tmp_coverage_stranded_.at(0).at( min(coverage_forward, kMaxCoverage) );
	++tmp_coverage_stranded_.at(1).at( min(coverage_reverse, kMaxCoverage) );
	if(cov_sum){
		++tmp_coverage_stranded_percent_.at(0).at( Percent(coverage_forward, cov_sum) );
		++tmp_coverage_stranded_percent_.at(1).at( Percent(coverage_reverse, cov_sum) );

		if(coverage.valid_){
			++tmp_error_coverage_.at( error_sum );
			++tmp_error_coverage_percent_.at( Percent(error_sum, cov_sum) );
			if(coverage_forward && coverage_reverse){
				++tmp_error_coverage_percent_stranded_.at( Percent(errors_forward, coverage_forward) ).at( Percent(errors_reverse, coverage_reverse) );
			}
		}

		if( 10 <= cov_sum ){
			++tmp_coverage_stranded_percent_min_cov_10_.at(0).at( Percent(coverage_forward, cov_sum) );
			++tmp_coverage_stranded_percent_min_cov_10_.at(1).at( Percent(coverage_reverse, cov_sum) );
			if(coverage.valid_){
				++tmp_error_coverage_percent_min_cov_10_.at( Percent(error_sum, cov_sum) );
			}

			if( 20 <= cov_sum ){
				++tmp_coverage_stranded_percent_min_cov_20_.at(0).at( Percent(coverage_forward, cov_sum) );
				++tmp_coverage_stranded_percent_min_cov_20_.at(1).at( Percent(coverage_reverse, cov_sum) );
				if(coverage.valid_){
					++tmp_error_coverage_percent_min_cov_20_.at( Percent(error_sum, cov_sum) );

					if(10 <=coverage_forward && 10 <= coverage_reverse){
						++tmp_error_coverage_percent_stranded_min_strand_cov_10_.at( Percent(errors_forward, coverage_forward) ).at( Percent(errors_reverse, coverage_reverse) );
						if(20 <=coverage_forward && 20 <= coverage_reverse){
							++tmp_error_coverage_percent_stranded_min_strand_cov_20_.at( Percent(errors_forward, coverage_forward) ).at( Percent(errors_reverse, coverage_reverse) );
						}
					}
				}
			}
		}
	}
}

void CoverageStats::UpdateDistances(uintSeqLen &distance_to_start_of_error_region, uintPercent &start_rate, uintPercent error_rate) const{
	if( distance_to_start_of_error_region ){
		if(start_rate < error_rate){
			distance_to_start_of_error_region = 0;
			start_rate = error_rate;
		}
		else if( ++distance_to_start_of_error_region >= reset_distance_ ){
			distance_to_start_of_error_region = 0;
			start_rate = 0;
		}
	}
	else{
		if(error_rate){ // If we assign an error rate it is systematic
			distance_to_start_of_error_region = 1;
			start_rate = error_rate;
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
			new_block->previous_coverage_.clear();
			new_block->reads_.clear();
		}
		else{
			reuse_mutex_.unlock();
			new_block = new CoverageBlock(seq_id, start_pos, last_block_);
			new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
		}
	}
	else{
		new_block = new CoverageBlock(seq_id, start_pos, last_block_);
		new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
	}

	new_block->reads_.reserve( (*last_block_).reads_.capacity() );
	new_block->coverage_.resize( kBlockSize );
	new_block->first_variant_id_ = (*last_block_).first_variant_id_;
	(*last_block_).next_block_ = new_block;

	return new_block;
}

void CoverageStats::UpdateFirstVariant(CoverageBlock &block, const Reference &reference){
	if( reference.VariantPositionsLoaded() ){
		while( block.first_variant_id_ < reference.VariantPositions(block.sequence_id_).size() && reference.VariantPositions(block.sequence_id_).at(block.first_variant_id_) < (*last_block_).start_pos_ ){
			++block.first_variant_id_;
		}
	}
}

bool CoverageStats::IsVariantPosition( intVariantId &cur_var, uintRefSeqId ref_seq_id, uintSeqLen pos, const Reference &reference ) const{
	if( reference.VariantPositionsLoaded() ){
		while(cur_var < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(cur_var) < pos){
			++cur_var;
		}

		return cur_var < reference.VariantPositions(ref_seq_id).size() && reference.VariantPositions(ref_seq_id).at(cur_var) == pos;
	}

	return false;
}

void CoverageStats::InitBlock(CoverageBlock &block, const Reference &reference){
	// Ensure that the block is not reaching over the end of the reference sequence
	if( reference.SequenceLength(block.sequence_id_) < block.start_pos_+kBlockSize){
		block.coverage_.resize( reference.SequenceLength(block.sequence_id_)-block.start_pos_ );
	}

	// Check for invalid positions
	auto cur_var = block.first_variant_id_;
	for(auto pos = block.start_pos_; pos < block.start_pos_+block.coverage_.size(); ++pos){
		if( IsN( at(reference.ReferenceSequence(block.sequence_id_), pos) ) || IsVariantPosition( cur_var, block.sequence_id_, pos, reference  ) ){
			block.coverage_.at(pos-block.start_pos_).valid_ = false;
		}
	}
}

void CoverageStats::ProcessBlock(CoverageBlock *block, const Reference &reference){
	auto non_systematic_error = GetMaxNonSystematicError(block, reference);
	for( uintSeqLen pos = 0; pos < block->coverage_.size(); ++pos ){
		UpdateCoverageAtSinglePosition( block->coverage_.at(pos), at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos), non_systematic_error );
	}

	// Save information that is still needed in next_block_
	if(block->next_block_ && block->start_pos_+kBlockSize == (*block->next_block_).start_pos_ && block->sequence_id_ == (*block->next_block_).sequence_id_){
		for(uintReadLen pos_before = 1; pos_before < maximum_read_length_on_reference_+1; ++pos_before){
			(*block->next_block_).previous_coverage_.emplace_back( block->coverage_.at(kBlockSize-pos_before) );
		}
	}
}

void CoverageStats::CountBlock(CoverageBlock *block, const Reference &reference){
	// Enter error_rates_ and dominant_errors_ in forward direction
	Dna5 ref_base, prev_ref_base(4);
	DominantBase dom_ref_base;
	if(block->start_pos_){
		prev_ref_base = at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ - 1);
		dom_ref_base.Set( reference.ReferenceSequence(block->sequence_id_), block->start_pos_ );
	}
	uintSeqLen gc_bases = min(block->start_pos_, gc_range_);
	uintSeqLen n_count(0);
	uintSeqLen gc = reference.GCContentAbsolut( n_count, block->sequence_id_, block->start_pos_-gc_bases, block->start_pos_);
	uintPercent gc_percent;

	for( uintSeqLen pos = 0; pos < tmp_errors_forward_.size(); ++pos ){
		ref_base = at(reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos);
		gc_percent = SafePercent(gc,gc_bases-n_count);

		// Enter values
		if(tmp_errors_forward_.at(pos).valid_and_coverage_sufficient_){
			auto distance = TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_forward_);

			++tmp_dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(tmp_errors_forward_.at(pos).base_);
			++tmp_dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(gc_percent).at(tmp_errors_forward_.at(pos).base_);
			++tmp_gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(gc_percent);
			++tmp_dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(start_rate_forward_).at(tmp_errors_forward_.at(pos).base_);
			++tmp_start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(start_rate_forward_);
			++tmp_start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(gc_percent).at(start_rate_forward_);

			++tmp_error_rates_by_distance_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(distance).at(tmp_errors_forward_.at(pos).rate_);
			++tmp_error_rates_by_gc_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(gc_percent).at(tmp_errors_forward_.at(pos).rate_);
			++tmp_gc_by_distance_er_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(distance).at(gc_percent);
			++tmp_error_rates_by_start_rates_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(start_rate_forward_).at(tmp_errors_forward_.at(pos).rate_);
			++tmp_start_rates_by_distance_er_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(distance).at(start_rate_forward_);
			++tmp_start_rates_by_gc_er_.at(ref_base).at(tmp_errors_forward_.at(pos).base_).at(gc_percent).at(start_rate_forward_);
		}

		// Set values for next iteration
		prev_ref_base = ref_base;
		dom_ref_base.Update( ref_base, reference.ReferenceSequence(block->sequence_id_), block->start_pos_ + pos);
		UpdateDistances(distance_to_start_of_error_region_forward_, start_rate_forward_, (tmp_errors_forward_.at(pos).valid_and_coverage_sufficient_?tmp_errors_forward_.at(pos).rate_:0));
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
			start_rate_forward_ = 0;
		}
		else{
			distance_to_start_of_error_region_forward_ += (*(block->next_block_)).start_pos_ - block->start_pos_ - block->coverage_.size(); // Take the uncovered region between blocks into account
			if( reset_distance_ <= distance_to_start_of_error_region_forward_ ){
				// Error region ends
				distance_to_start_of_error_region_forward_ = 0;
				start_rate_forward_ = 0;
			}
		}
	}

	// Enter error_rates_ and dominant_errors_ in reverse direction
	intSeqShift last_pos = tmp_errors_reverse_.size(); // Might get negative when we subtract
	if(block->next_block_ && (*(block->next_block_)).sequence_id_ == block->sequence_id_){
		// Everything that is potentially in an error region starting in the next block cannot be handled yet
		last_pos -= reset_distance_+1;
	}
	else{
		// Reference sequence changes between this block and the next or we are at the end of the file: We can handle all positions
		--last_pos; // Shift it from end pos to last valid pos
	}

	if(0 <= last_pos){
		ConstDna5StringReverseComplement reversed_seq(reference.ReferenceSequence(block->sequence_id_));
		uintSeqLen ref_pos = reference.SequenceLength(block->sequence_id_) + tmp_errors_reverse_.size() - block->start_pos_ - block->coverage_.size() - last_pos - 1; // tmp_errors_reverse_.size() - block->coverage_.size() is the shift to the left of block->start_pos_ due to the non-handleable reverse errors taken from the previous block
		dom_ref_base.Clear();
		if(ref_pos){
			prev_ref_base = at(reversed_seq, ref_pos - 1);
			dom_ref_base.Set( reversed_seq, ref_pos );
		}
		else{
			prev_ref_base = 4;
		}
		gc_bases = min(ref_pos, gc_range_);
		n_count = 0;
		gc = reference.GCContentAbsolut( n_count, block->sequence_id_, reference.SequenceLength(block->sequence_id_)-ref_pos, reference.SequenceLength(block->sequence_id_)-ref_pos+gc_bases);

		// Initialize distance and start_rate
		uintSeqLen distance_to_start_of_error_region_reverse(0);
		uintPercent start_rate_reverse(0);
		for( uintSeqLen pos = tmp_errors_reverse_.size(); pos-- > last_pos; ){
			UpdateDistances(distance_to_start_of_error_region_reverse, start_rate_reverse, (tmp_errors_reverse_.at(pos).valid_and_coverage_sufficient_?tmp_errors_reverse_.at(pos).rate_:0));
		}

		for( uintSeqLen pos = last_pos+1; pos--; ){
			ref_base = at(reversed_seq, ref_pos);
			gc_percent = SafePercent(gc,gc_bases-n_count);

			// Enter values
			if(tmp_errors_reverse_.at(pos).valid_and_coverage_sufficient_){
				auto distance = TransformDistanceToStartOfErrorRegion(distance_to_start_of_error_region_reverse);

				++tmp_dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(tmp_errors_reverse_.at(pos).base_);
				++tmp_dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(gc_percent).at(tmp_errors_reverse_.at(pos).base_);
				++tmp_gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(gc_percent);
				++tmp_dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(start_rate_reverse).at(tmp_errors_reverse_.at(pos).base_);
				++tmp_start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(distance).at(start_rate_reverse);
				++tmp_start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_ref_base.Get()).at(gc_percent).at(start_rate_reverse);

				++tmp_error_rates_by_distance_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(distance).at(tmp_errors_reverse_.at(pos).rate_);
				++tmp_error_rates_by_gc_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(gc_percent).at(tmp_errors_reverse_.at(pos).rate_);
				++tmp_gc_by_distance_er_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(distance).at(gc_percent);
				++tmp_error_rates_by_start_rates_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(start_rate_reverse).at(tmp_errors_reverse_.at(pos).rate_);
				++tmp_start_rates_by_distance_er_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(distance).at(start_rate_reverse);
				++tmp_start_rates_by_gc_er_.at(ref_base).at(tmp_errors_reverse_.at(pos).base_).at(gc_percent).at(start_rate_reverse);
			}

			// Set values for next iteration
			prev_ref_base = ref_base;
			dom_ref_base.Update( ref_base, reversed_seq, ref_pos);
			UpdateDistances(distance_to_start_of_error_region_reverse, start_rate_reverse, (tmp_errors_reverse_.at(pos).valid_and_coverage_sufficient_?tmp_errors_reverse_.at(pos).rate_:0));
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

void CoverageStats::Prepare(uintCovCount average_coverage, uintReadLen average_read_length, uintReadLen maximum_read_length_on_reference){
	reusable_blocks_.reserve(100);
	tmp_errors_forward_.reserve(2*kBlockSize);
	tmp_errors_reverse_.reserve(2*kBlockSize);
	distance_to_start_of_error_region_forward_ = 0;
	start_rate_forward_ = 0;

	if(coverage_threshold_ > average_coverage/4 ){
		if( 1 < average_coverage/4 ){
			coverage_threshold_ = average_coverage/4; // If half the average coverage is lower than the threshold, so we do not just take the statistics from repeat regions (the coverage per strand is already half the total coverage, so in total we have a quarter)
		}
		else{
			coverage_threshold_ = 1; // Do not set anything without having coverage
		}
	}
	gc_range_ = average_read_length/2;
	reset_distance_ = average_read_length;

	maximum_read_length_on_reference_ = maximum_read_length_on_reference;

	auto max_error_dist = TransformDistanceToStartOfErrorRegion(reset_distance_-1)+1;
	for( auto ref_base=4; ref_base--; ){
		for( auto prev_ref_base=5; prev_ref_base--; ){
			for( auto dom_base=4; dom_base--; ){
				SetDimensions( tmp_dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_base), max_error_dist, 5 );
				SetDimensions( tmp_dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_base), 101, 5 );
				SetDimensions( tmp_gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base), max_error_dist, 101 );
				SetDimensions( tmp_dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_base), 101, 5 );
				SetDimensions( tmp_start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base), max_error_dist, 101 );
				SetDimensions( tmp_start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_base), 101, 101 );
			}
		}

		for( auto dom_error=5; dom_error--; ){
			SetDimensions( tmp_error_rates_by_distance_.at(ref_base).at(dom_error), max_error_dist, 101 );
			SetDimensions( tmp_error_rates_by_gc_.at(ref_base).at(dom_error), 101, 101 );
			SetDimensions( tmp_gc_by_distance_er_.at(ref_base).at(dom_error), max_error_dist, 101 );
			SetDimensions( tmp_error_rates_by_start_rates_.at(ref_base).at(dom_error), 101, 101 );
			SetDimensions( tmp_start_rates_by_distance_er_.at(ref_base).at(dom_error), max_error_dist, 101 );
			SetDimensions( tmp_start_rates_by_gc_er_.at(ref_base).at(dom_error), 101, 101 );
		}
	}

	// Setting the size, those aren't atomics
	block_error_rate_[100];
	block_non_systematic_error_rate_[100];

	tmp_coverage_.resize(kMaxCoverage+1);
	for( auto strand=2; strand--; ){
		tmp_coverage_stranded_.at(strand).resize(kMaxCoverage+1);
		tmp_coverage_stranded_percent_.at(strand).resize(101);
		tmp_coverage_stranded_percent_min_cov_10_.at(strand).resize(101);
		tmp_coverage_stranded_percent_min_cov_20_.at(strand).resize(101);
	}
	tmp_error_coverage_.resize(kMaxCoverage+1);
	tmp_error_coverage_percent_.resize(101);
	tmp_error_coverage_percent_min_cov_10_.resize(101);
	tmp_error_coverage_percent_min_cov_20_.resize(101);
	SetDimensions( tmp_error_coverage_percent_stranded_, 101, 101 );
	SetDimensions( tmp_error_coverage_percent_stranded_min_strand_cov_10_, 101, 101 );
	SetDimensions( tmp_error_coverage_percent_stranded_min_strand_cov_20_, 101, 101 );
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

bool CoverageStats::EnsureSpace(uintRefSeqId ref_seq_id, uintSeqLen start_pos, uintSeqLen end_pos, FullRecord *record, Reference &reference){
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
			CountEmptyEndOfSequence(reference);

			// Account for the not covered reference sequences in between
			for( auto seq_id = ref_seq_id; --seq_id > (*last_block_).sequence_id_ ; ){
				zero_coverage_region_ += reference.SequenceLength(seq_id);
				excluded_bases_ += reference.SumExcludedBases(seq_id);
				num_exclusion_regions_ += reference.NumExcludedRegions(seq_id);
			}

			// Account for the not covered region in the current reference sequence before the first position
			zero_coverage_region_ += start_pos;
			excluded_bases_ += reference.SumExcludedBases(ref_seq_id);
			num_exclusion_regions_ += reference.NumExcludedRegions(ref_seq_id);

			// Create first block of next reference sequence
			last_block_ = CreateBlock(ref_seq_id, start_pos);
			InitBlock(*last_block_, reference);
		}
		else if((*last_block_).start_pos_+kBlockSize < start_pos){
			// Account for the empty coverage region
			zero_coverage_region_ += start_pos - (*last_block_).start_pos_ - kBlockSize;

			// Create first block after the empty coverage region
			last_block_ = CreateBlock(ref_seq_id, start_pos);
			UpdateFirstVariant(*last_block_, reference);
			InitBlock(*last_block_, reference);
		}
	}
	else{
		// Initialize first block
		// Account for the not covered reference sequences before the first position
		for( auto seq_id = ref_seq_id; seq_id-- ; ){
			zero_coverage_region_ += reference.SequenceLength(seq_id);
			excluded_bases_ += reference.SumExcludedBases(seq_id);
			num_exclusion_regions_ += reference.NumExcludedRegions(seq_id);
		}

		// Account for the not covered region in the current reference sequence before the first position
		zero_coverage_region_ += start_pos;
		excluded_bases_ += reference.SumExcludedBases(ref_seq_id);
		num_exclusion_regions_ += reference.NumExcludedRegions(ref_seq_id);

		// Create new block
		first_block_ = new CoverageBlock(ref_seq_id, start_pos, NULL);
		first_block_->coverage_.resize( kBlockSize );
		first_block_->previous_coverage_.reserve(maximum_read_length_on_reference_);
		first_block_->reads_.reserve( 2*kBlockSize );

		last_block_ = first_block_;
		InitBlock(*last_block_, reference);
	}

	// Create blocks until end_pos is reached
	while((*last_block_).start_pos_+kBlockSize < end_pos){
		last_block_ = CreateBlock((*last_block_).sequence_id_, (*last_block_).start_pos_+kBlockSize);
		UpdateFirstVariant(*last_block_, reference);
		InitBlock(*last_block_, reference);
	}

	// Add read to last block it potentially overlaps (most of the time it is last_block_, but it is not guaranteed so use FindBlock)
	auto reg_block = FindBlock(ref_seq_id, end_pos);
	reg_block->reads_.push_back(record);

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
			ProcessBlock(until_block, reference);
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

bool CoverageStats::Finalize(const Reference &reference, QualityStats &qualities, ErrorStats &errors, uintQual phred_quality_offset, mutex &print_mutex){
	// Remove all reusable_blocks_ as they are not needed anymore
	for( auto block : reusable_blocks_ ){
		delete block;
	}
	auto final_num_blocks = reusable_blocks_.size();
	reusable_blocks_.clear();

	if(first_block_){
		if( last_block_ ){
			// Update not covered regions after last block
			CountEmptyEndOfSequence(reference);

			for( auto seq_id = (*last_block_).sequence_id_+1; seq_id < reference.NumberSequences(); ++seq_id ){
				zero_coverage_region_ += reference.SequenceLength(seq_id);
				excluded_bases_ += reference.SumExcludedBases(seq_id);
				num_exclusion_regions_ += reference.NumExcludedRegions(seq_id);
			}

			// Update coverage of remaining blocks and delete them
			while(first_block_->next_block_){
				++final_num_blocks;
				ProcessBlock(first_block_, reference);
				CountBlock(first_block_, reference);
				for( auto rec : first_block_->reads_){
					EvalRead(rec, first_block_, reference, qualities, errors, phred_quality_offset);
				}

				first_block_ = first_block_->next_block_;
				delete first_block_->previous_block_;
			}

			++final_num_blocks;
			ProcessBlock(first_block_, reference);
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

	{
		lock_guard<mutex> lock(print_mutex);
		printInfo << "Needed to create " << final_num_blocks << " coverage blocks." << std::endl;
		auto total_size = reference.TotalSize();
		printInfo << "Excluded " << excluded_bases_ << " of " << total_size << " bases [" << static_cast<uintPercentPrint>(Percent(excluded_bases_, total_size)) << "%] in the reference due to N's distributed over " << num_exclusion_regions_ << " regions." << std::endl;
	}

	ApplyZeroCoverageRegion();
	tmp_errors_forward_.shrink_to_fit();
	tmp_errors_reverse_.shrink_to_fit();

	for( auto ref_base=4; ref_base--; ){
		for( auto prev_ref_base=5; prev_ref_base--; ){
			for( auto dom_base=4; dom_base--; ){
				dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_base) );
				dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_base) );
				gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
				dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_base) );
				start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
				start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_base).Acquire( tmp_start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
			}
		}

		for( auto dom_error=5; dom_error--; ){
			error_rates_by_distance_.at(ref_base).at(dom_error).Acquire( tmp_error_rates_by_distance_.at(ref_base).at(dom_error) );
			error_rates_by_gc_.at(ref_base).at(dom_error).Acquire( tmp_error_rates_by_gc_.at(ref_base).at(dom_error) );
			gc_by_distance_er_.at(ref_base).at(dom_error).Acquire( tmp_gc_by_distance_er_.at(ref_base).at(dom_error) );
			error_rates_by_start_rates_.at(ref_base).at(dom_error).Acquire( tmp_error_rates_by_start_rates_.at(ref_base).at(dom_error) );
			start_rates_by_distance_er_.at(ref_base).at(dom_error).Acquire( tmp_start_rates_by_distance_er_.at(ref_base).at(dom_error) );
			start_rates_by_gc_er_.at(ref_base).at(dom_error).Acquire( tmp_start_rates_by_gc_er_.at(ref_base).at(dom_error) );
		}
	}

	coverage_.Acquire( tmp_coverage_ );
	for( auto strand=2; strand--; ){
		coverage_stranded_.at(strand).Acquire( tmp_coverage_stranded_.at(strand) );
		coverage_stranded_percent_.at(strand).Acquire( tmp_coverage_stranded_percent_.at(strand) );
		coverage_stranded_percent_min_cov_10_.at(strand).Acquire( tmp_coverage_stranded_percent_min_cov_10_.at(strand) );
		coverage_stranded_percent_min_cov_20_.at(strand).Acquire( tmp_coverage_stranded_percent_min_cov_20_.at(strand) );
	}
	error_coverage_.Acquire( tmp_error_coverage_ );
	error_coverage_percent_.Acquire( tmp_error_coverage_percent_ );
	error_coverage_percent_min_cov_10_.Acquire( tmp_error_coverage_percent_min_cov_10_ );
	error_coverage_percent_min_cov_20_.Acquire( tmp_error_coverage_percent_min_cov_20_ );
	error_coverage_percent_stranded_.Acquire( tmp_error_coverage_percent_stranded_ );
	error_coverage_percent_stranded_min_strand_cov_10_.Acquire( tmp_error_coverage_percent_stranded_min_strand_cov_10_ );
	error_coverage_percent_stranded_min_strand_cov_20_.Acquire( tmp_error_coverage_percent_stranded_min_strand_cov_20_ );

	return true;
}

void CoverageStats::Shrink(){
	for( auto ref_base=4; ref_base--; ){
		for( auto prev_ref_base=5; prev_ref_base--; ){
			for( auto dom_base=4; dom_base--; ){
				ShrinkVect( dominant_errors_by_distance_.at(ref_base).at(prev_ref_base).at(dom_base) );
				ShrinkVect( dominant_errors_by_gc_.at(ref_base).at(prev_ref_base).at(dom_base) );
				ShrinkVect( gc_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
				ShrinkVect( dominant_errors_by_start_rates_.at(ref_base).at(prev_ref_base).at(dom_base) );
				ShrinkVect( start_rates_by_distance_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
				ShrinkVect( start_rates_by_gc_de_.at(ref_base).at(prev_ref_base).at(dom_base) );
			}
		}

		for( auto dom_error=5; dom_error--; ){
			ShrinkVect( error_rates_by_distance_.at(ref_base).at(dom_error) );
			ShrinkVect( error_rates_by_gc_.at(ref_base).at(dom_error) );
			ShrinkVect( gc_by_distance_er_.at(ref_base).at(dom_error) );
			ShrinkVect( error_rates_by_start_rates_.at(ref_base).at(dom_error) );
			ShrinkVect( start_rates_by_distance_er_.at(ref_base).at(dom_error) );
			ShrinkVect( start_rates_by_gc_er_.at(ref_base).at(dom_error) );
		}
	}

	block_error_rate_.Shrink();
	block_non_systematic_error_rate_.Shrink();

	ShrinkVect( error_coverage_percent_stranded_ );
	ShrinkVect( error_coverage_percent_stranded_min_strand_cov_10_ );
	ShrinkVect( error_coverage_percent_stranded_min_strand_cov_20_ );
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
