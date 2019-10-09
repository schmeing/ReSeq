#include "QualityStats.h"
using reseq::QualityStats;
using reseq::SeqQualityStats;

#include "reportingUtils.hpp"
//include "utilities.hpp"
using reseq::utilities::Divide;
using reseq::utilities::Percent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;

void QualityStats::SplitPairedSequenceQuality(){
	// Split sequence_quality_mean_paired_per_tile_ by template_segment
	for( uintTempSeq template_segment=2; template_segment--; ){
		sequence_quality_mean_per_tile_.at(template_segment).Clear();
	}

	for( auto tile_id = sequence_quality_mean_paired_per_tile_.size(); tile_id--; ){
		for( auto sq_first = sequence_quality_mean_paired_per_tile_.at(tile_id).to(); sq_first-- > sequence_quality_mean_paired_per_tile_.at(tile_id).from();){
			for( auto sq_second = sequence_quality_mean_paired_per_tile_.at(tile_id).at(sq_first).to(); sq_second-- > sequence_quality_mean_paired_per_tile_.at(tile_id).at(sq_first).from();){
				sequence_quality_mean_per_tile_.at(0)[tile_id][sq_first] += sequence_quality_mean_paired_per_tile_.at(tile_id).at(sq_first).at(sq_second);
				sequence_quality_mean_per_tile_.at(1)[tile_id][sq_second] += sequence_quality_mean_paired_per_tile_.at(tile_id).at(sq_first).at(sq_second);

				// Sum up tiles
				sequence_quality_mean_paired_[sq_first][sq_second] += sequence_quality_mean_paired_per_tile_.at(tile_id).at(sq_first).at(sq_second);
			}
		}
	}
	ShrinkVect(sequence_quality_mean_paired_);

	for( uintTempSeq template_segment=2; template_segment--; ){
		ShrinkVect(sequence_quality_mean_per_tile_.at(template_segment));
	}
}

void QualityStats::SumTiles(){
	Vect<SeqQualityStats<uintFragCount>> sequence_quality_tile_sum;

	for( uintTempSeq template_segment=2; template_segment--; ){
		// From Reference
		base_quality_stats_reference_.at(template_segment).Clear();
		for( auto ref_base = base_quality_stats_per_tile_per_error_reference_.at(template_segment).size(); ref_base--; ){
			for( auto dom_error = base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).size(); dom_error--; ){
				for( auto tile_id = base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).from(); tile_id < base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).to(); ++tile_id ){
					base_quality_stats_reference_.at(template_segment) += base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error)[tile_id];
				}
			}
		}
		ShrinkVect(base_quality_stats_reference_.at(template_segment));

		sequence_quality_tile_sum.Clear();
		for( auto tile_id = sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment).size(); tile_id--; ){
			sequence_quality_tile_sum += sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment).at(tile_id);
		}
		average_sequence_quality_for_gc_.at(template_segment).Clear();
		for( uintPercent gc_percentage = sequence_quality_tile_sum.to(); gc_percentage-- > sequence_quality_tile_sum.from(); ){
			sequence_quality_tile_sum.at(gc_percentage).Calculate(false);
			average_sequence_quality_for_gc_.at(template_segment)[gc_percentage] = sequence_quality_tile_sum.at(gc_percentage).mean_;
		}
		average_sequence_quality_for_gc_.at(template_segment).Shrink();

		// From raw reads
		sequence_quality_mean_.at(template_segment).Clear();
		sequence_quality_mean_.at(template_segment).SetOffset(2);
		for( auto tile_id = sequence_quality_mean_per_tile_.at(template_segment).size(); tile_id--; ){
			sequence_quality_mean_.at(template_segment) += sequence_quality_mean_per_tile_.at(template_segment).at(tile_id);
		}
		sequence_quality_mean_.at(template_segment).Shrink();

		base_quality_stats_.at(template_segment).Clear();
		base_quality_for_sequence_.at(template_segment).Clear();
		base_quality_for_preceding_quality_.at(template_segment).Clear();
		for( auto called_base = 5; called_base--; ){
			sequence_quality_tile_sum.Clear();
			for( auto tile_id = base_quality_stats_per_tile_.at(template_segment).at(called_base).size(); tile_id--; ){
				base_quality_stats_.at(template_segment) += base_quality_stats_per_tile_.at(template_segment).at(called_base)[tile_id];

				base_quality_for_sequence_.at(template_segment) += base_quality_for_sequence_per_tile_.at(template_segment).at(called_base)[tile_id];
				base_quality_for_preceding_quality_.at(template_segment) += base_quality_for_preceding_quality_per_tile_.at(template_segment).at(called_base)[tile_id];

				sequence_quality_tile_sum += sequence_quality_for_base_per_tile_.at(template_segment).at(called_base)[tile_id];
			}

			average_sequence_quality_for_base_.at(template_segment).at(called_base).Clear();
			for( uintPercent nuc_percentage = sequence_quality_tile_sum.to(); nuc_percentage-- > sequence_quality_tile_sum.from(); ){
				sequence_quality_tile_sum.at(nuc_percentage).Calculate(false);
				average_sequence_quality_for_base_.at(template_segment).at(called_base)[nuc_percentage] = sequence_quality_tile_sum.at(nuc_percentage).mean_;
			}
			average_sequence_quality_for_base_.at(template_segment).at(called_base).Shrink();
		}
		ShrinkVect(base_quality_stats_.at(template_segment));
		ShrinkVect(base_quality_for_sequence_.at(template_segment));
		ShrinkVect(base_quality_for_preceding_quality_.at(template_segment));
	}
}

void QualityStats::CalculateQualityStats(){
	// Calculate quality means, etc.
	SeqQualityStats<uintNucCount> *stats;
	char mean_difference;
	SeqQualityStats<uintNucCount> base_quality_stats_tile_sum;
	for( uintTempSeq template_segment=2; template_segment--; ){
		for( auto pos = base_quality_stats_reference_.at(template_segment).size(); pos--; ){
			// base_quality_stats_reference_
			stats = &base_quality_stats_reference_.at(template_segment).at(pos);
			stats->Calculate();

			base_quality_mean_reference_.at(template_segment)[pos] = stats->mean_;
			base_quality_minimum_reference_.at(template_segment)[pos] = stats->minimum_;
			base_quality_first_quartile_reference_.at(template_segment)[pos] = stats->first_quartile_;
			base_quality_median_reference_.at(template_segment)[pos] = stats->median_;
			base_quality_third_quartile_reference_.at(template_segment)[pos] = stats->third_quartile_;
			base_quality_maximum_reference_.at(template_segment)[pos] = stats->maximum_;
		}

		for( auto pos = base_quality_stats_.at(template_segment).size(); pos--; ){
			// base_quality_stats_
			stats = &base_quality_stats_.at(template_segment).at(pos);
			stats->Calculate();

			base_quality_mean_.at(template_segment)[pos] = stats->mean_;
			base_quality_minimum_.at(template_segment)[pos] = stats->minimum_;
			base_quality_first_quartile_.at(template_segment)[pos] = stats->first_quartile_;
			base_quality_median_.at(template_segment)[pos] = stats->median_;
			base_quality_third_quartile_.at(template_segment)[pos] = stats->third_quartile_;
			base_quality_maximum_.at(template_segment)[pos] = stats->maximum_;

			// tile_quality_mean_difference_
			for( auto tile_id = base_quality_stats_per_tile_.at(template_segment).at(0).size(); tile_id--; ){
				base_quality_stats_tile_sum.Clear();
				for( auto base = base_quality_stats_per_tile_.at(template_segment).size(); base--; ){
					base_quality_stats_tile_sum += base_quality_stats_per_tile_.at(template_segment).at(base)[tile_id][pos];
				}
				base_quality_stats_tile_sum.Calculate();

				mean_difference = base_quality_stats_tile_sum.mean_ - stats->mean_;
				if( mean_difference ){ // Do not put in empty bins that are later removed again by shrinking
					tile_quality_mean_difference_.at(template_segment)[tile_id][pos] = mean_difference;
				}
			}
		}

		for( auto pos = base_quality_stats_per_strand_.at(template_segment).size(); pos--; ){
			// base_quality_stats_per_strand_
			stats = &base_quality_stats_per_strand_.at(template_segment).at(pos); // template_segment really is strand here
			stats->Calculate(false);

			base_quality_mean_per_strand_.at(template_segment)[pos] = stats->mean_; // template_segment really is strand here
		}

		// Fill mean_sequence_quality_mean_by_fragment_length_
		uintSeqLen min_frag_length(sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(0).from()), max_frag_length(0);
		for(uintTileId tile_id=sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).to(); tile_id--; ){
			SetToMin(min_frag_length, sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id).from());
			SetToMax(max_frag_length, sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id).to());
		}
		mean_sequence_quality_mean_by_fragment_length_.at(template_segment).Clear();
		uintFragCount counts(0), sum(0), last_counts(0), last_sum(0);
		uintSeqLen from, to(max_frag_length), last_to(max_frag_length);
		for(uintSeqLen frag_length=max_frag_length; frag_length-- > min_frag_length; ){
			from = frag_length;
			for(uintTileId tile_id=sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).to(); tile_id--; ){
				for(uintQual qual=sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id)[frag_length].from(); qual<sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id)[frag_length].to(); ++qual){
					counts += sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id)[frag_length][qual];
					sum += qual*sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).at(tile_id)[frag_length][qual];
				}
			}

			if(100 <= counts){ // Bin fragment lengths together to take the mean over at least 100 sequences
				for(uintSeqLen len = to; len-- > from;){
					mean_sequence_quality_mean_by_fragment_length_.at(template_segment)[len] = Divide(sum, counts);
				}
				last_counts = counts;
				counts = 0;
				last_sum = sum;
				sum = 0;
				last_to = to;
				to = frag_length;
			}
		}

		if(counts){ // Write out last fragment bin if it did not reach 100 counts
			for(uintSeqLen len = last_to; len-- > from;){
				mean_sequence_quality_mean_by_fragment_length_.at(template_segment)[len] = Divide(sum+last_sum, counts+last_counts);
			}
		}
	}
}

void QualityStats::Prepare(uintTileId num_tiles, uintQual size_qual, uintReadLen size_pos, uintSeqLen maximum_fragment_length){
	// Resize vectors to necessary size
	for( auto template_segment=2; template_segment--; ){
		for( auto ref_base=4; ref_base--; ){
			for( auto dom_error=5; dom_error--; ){
				SetDimensions( tmp_base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error), num_tiles, size_pos, size_qual );
				SetDimensions( tmp_error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error), num_tiles, size_pos, 101 );
				SetDimensions( tmp_base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error), num_tiles, 101, size_qual );
			}
			SetDimensions( tmp_base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, 101, size_qual );
			SetDimensions( tmp_preceding_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, size_pos, size_qual );
			SetDimensions( tmp_base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, 101, size_qual );
			SetDimensions( tmp_sequence_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base), num_tiles, size_pos, size_qual );
		}

		SetDimensions( tmp_sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment), num_tiles, 101, size_qual );
		SetDimensions( tmp_sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment), num_tiles, 101, size_qual );
		SetDimensions( tmp_sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment), num_tiles, maximum_fragment_length+1, size_qual );
		SetDimensions( tmp_mean_error_rate_for_gc_per_tile_reference_.at(template_segment), num_tiles, 101, 101 );
		SetDimensions( tmp_mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment), num_tiles, maximum_fragment_length+1, 101 );
		SetDimensions( tmp_gc_for_fragment_length_per_tile_reference_.at(template_segment), num_tiles, maximum_fragment_length+1, 101 );

		for( auto base=5; base--; ){
			SetDimensions( tmp_base_quality_for_sequence_per_tile_.at(template_segment).at(base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_base_quality_for_preceding_quality_per_tile_.at(template_segment).at(base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_base_quality_stats_per_tile_.at(template_segment).at(base), num_tiles, size_pos, size_qual );
			SetDimensions( tmp_preceding_quality_for_sequence_per_tile_.at(template_segment).at(base), num_tiles, size_qual, size_qual );
			SetDimensions( tmp_preceding_quality_for_position_per_tile_.at(template_segment).at(base), num_tiles, size_pos, size_qual );
			SetDimensions( tmp_sequence_quality_for_position_per_tile_.at(template_segment).at(base), num_tiles, size_pos, size_qual );

			SetDimensions( tmp_sequence_quality_for_base_per_tile_.at(template_segment).at(base), num_tiles, 101, size_qual );

			tmp_nucleotide_quality_.at(template_segment).at(base).resize(size_qual);
		}

		SetDimensions( tmp_base_quality_stats_per_strand_.at(template_segment), size_pos, size_qual );

		SetDimensions( tmp_sequence_quality_mean_for_gc_per_tile_.at(template_segment), num_tiles, 101, size_qual );
		tmp_sequence_quality_probability_mean_.at(template_segment).resize(size_qual);
		tmp_sequence_quality_minimum_.at(template_segment).resize(size_qual);
		tmp_sequence_quality_first_quartile_.at(template_segment).resize(size_qual);
		tmp_sequence_quality_median_.at(template_segment).resize(size_qual);
		tmp_sequence_quality_third_quartile_.at(template_segment).resize(size_qual);
		tmp_sequence_quality_maximum_.at(template_segment).resize(size_qual);
		SetDimensions( tmp_sequence_quality_content_.at(template_segment), size_qual, size_pos );
	}

	SetDimensions( tmp_sequence_quality_mean_paired_per_tile_, num_tiles, size_qual, size_qual );

	SetDimensions( tmp_homoquality_distribution_, size_qual, size_pos );
}

void QualityStats::Finalize(uintFragCount total_number_reads){
	// Copy vectors to final ones
	for( auto template_segment=2; template_segment--; ){
		for( auto ref_base=4; ref_base--; ){
			for( auto dom_error=5; dom_error--; ){
				base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error) );
				error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error) );
				base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error) );
			}
			base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(ref_base) );
			preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base) );
			preceding_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_preceding_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base) );
			base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base) );
			preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base) );
			sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base) );
			sequence_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base).Acquire( tmp_sequence_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base) );
		}

		sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment).Acquire( tmp_sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment) );
		sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment).Acquire( tmp_sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment) );
		sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment).Acquire( tmp_sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment) );
		mean_error_rate_for_gc_per_tile_reference_.at(template_segment).Acquire( tmp_mean_error_rate_for_gc_per_tile_reference_.at(template_segment) );
		mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment).Acquire( tmp_mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment) );
		gc_for_fragment_length_per_tile_reference_.at(template_segment).Acquire( tmp_gc_for_fragment_length_per_tile_reference_.at(template_segment) );

		for( auto base=5; base--; ){
			base_quality_for_sequence_per_tile_.at(template_segment).at(base).Acquire( tmp_base_quality_for_sequence_per_tile_.at(template_segment).at(base) );
			base_quality_for_preceding_quality_per_tile_.at(template_segment).at(base).Acquire( tmp_base_quality_for_preceding_quality_per_tile_.at(template_segment).at(base) );
			base_quality_stats_per_tile_.at(template_segment).at(base).Acquire( tmp_base_quality_stats_per_tile_.at(template_segment).at(base) );
			preceding_quality_for_sequence_per_tile_.at(template_segment).at(base).Acquire( tmp_preceding_quality_for_sequence_per_tile_.at(template_segment).at(base) );
			preceding_quality_for_position_per_tile_.at(template_segment).at(base).Acquire( tmp_preceding_quality_for_position_per_tile_.at(template_segment).at(base) );
			sequence_quality_for_position_per_tile_.at(template_segment).at(base).Acquire( tmp_sequence_quality_for_position_per_tile_.at(template_segment).at(base) );

			sequence_quality_for_base_per_tile_.at(template_segment).at(base).Acquire( tmp_sequence_quality_for_base_per_tile_.at(template_segment).at(base) );

			nucleotide_quality_.at(template_segment).at(base).Acquire( tmp_nucleotide_quality_.at(template_segment).at(base) );
		}

		base_quality_stats_per_strand_.at(template_segment).Acquire( tmp_base_quality_stats_per_strand_.at(template_segment) );

		sequence_quality_mean_for_gc_per_tile_.at(template_segment).Acquire( tmp_sequence_quality_mean_for_gc_per_tile_.at(template_segment) );
		sequence_quality_probability_mean_.at(template_segment).Acquire( tmp_sequence_quality_probability_mean_.at(template_segment) );
		sequence_quality_minimum_.at(template_segment).Acquire( tmp_sequence_quality_minimum_.at(template_segment) );
		sequence_quality_first_quartile_.at(template_segment).Acquire( tmp_sequence_quality_first_quartile_.at(template_segment) );
		sequence_quality_median_.at(template_segment).Acquire( tmp_sequence_quality_median_.at(template_segment) );
		sequence_quality_third_quartile_.at(template_segment).Acquire( tmp_sequence_quality_third_quartile_.at(template_segment) );
		sequence_quality_maximum_.at(template_segment).Acquire( tmp_sequence_quality_maximum_.at(template_segment) );
		sequence_quality_content_.at(template_segment).Acquire( tmp_sequence_quality_content_.at(template_segment) );
	}

	sequence_quality_mean_paired_per_tile_.Acquire( tmp_sequence_quality_mean_paired_per_tile_ );

	homoquality_distribution_.Acquire( tmp_homoquality_distribution_);

	// Update sequence qualities that did not appear in some reads, so the zero appearance values are correct
	for( uintTempSeq template_segment=2; template_segment--; ){
		for( auto qual = sequence_quality_content_.at(template_segment).from(); qual < sequence_quality_content_.at(template_segment).to(); ++qual){
			if(sequence_quality_content_.at(template_segment).at(qual).size()){
				sequence_quality_content_.at(template_segment).at(qual)[0] += total_number_reads - SumVect(sequence_quality_content_.at(template_segment).at(qual));
			}
		}
	}
}

void QualityStats::Shrink(){
	for( uintTempSeq template_segment=2; template_segment--; ){
		for( auto ref_base = 4; ref_base--; ){
			for( auto dom_error = 5; dom_error--; ){
				ShrinkVect(base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error));
			}
			ShrinkVect(base_quality_for_preceding_quality_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(preceding_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(preceding_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(base_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(preceding_quality_for_sequence_quality_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(sequence_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(sequence_quality_for_position_per_tile_reference_.at(template_segment).at(ref_base));
		}

		ShrinkVect(sequence_quality_mean_for_gc_per_tile_reference_.at(template_segment));
		ShrinkVect(sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(template_segment));
		ShrinkVect(sequence_quality_mean_for_fragment_length_per_tile_reference_.at(template_segment));
		ShrinkVect(mean_error_rate_for_gc_per_tile_reference_.at(template_segment));
		ShrinkVect(mean_error_rate_for_fragment_length_per_tile_reference_.at(template_segment));
		ShrinkVect(gc_for_fragment_length_per_tile_reference_.at(template_segment));

		for(auto called_base = 5; called_base--; ){
			ShrinkVect(base_quality_for_sequence_per_tile_.at(template_segment).at(called_base));
			ShrinkVect(base_quality_for_preceding_quality_per_tile_.at(template_segment).at(called_base));
			ShrinkVect(preceding_quality_for_sequence_per_tile_.at(template_segment).at(called_base));
			ShrinkVect(preceding_quality_for_position_per_tile_.at(template_segment).at(called_base));
			ShrinkVect(sequence_quality_for_position_per_tile_.at(template_segment).at(called_base));

			ShrinkVect(base_quality_stats_per_tile_.at(template_segment).at(called_base));

			ShrinkVect(sequence_quality_for_base_per_tile_.at(template_segment).at(called_base));

			nucleotide_quality_.at(template_segment).at(called_base).Shrink();
		}

		ShrinkVect(base_quality_stats_per_strand_.at(template_segment));

		ShrinkVect(sequence_quality_mean_for_gc_per_tile_.at(template_segment));
		ShrinkVect(sequence_quality_probability_mean_.at(template_segment));
		sequence_quality_minimum_.at(template_segment).Shrink();
		sequence_quality_first_quartile_.at(template_segment).Shrink();
		sequence_quality_median_.at(template_segment).Shrink();
		sequence_quality_third_quartile_.at(template_segment).Shrink();
		sequence_quality_maximum_.at(template_segment).Shrink();
		ShrinkVect(sequence_quality_content_.at(template_segment));
	}

	ShrinkVect(sequence_quality_mean_paired_per_tile_);

	ShrinkVect(homoquality_distribution_);
}

void QualityStats::PrepareEstimation(){
	for( uintTempSeq template_segment=2; template_segment--; ){
		for( auto ref_base = base_quality_stats_per_tile_per_error_reference_.at(template_segment).size(); ref_base--; ){
			base_quality_stats_per_tile_reference_.at(template_segment).at(ref_base).Clear();
			error_rate_for_position_per_tile_reference_.at(template_segment).at(ref_base).Clear();
			base_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base).Clear();

			for( auto dom_error = base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).size(); dom_error--; ){
				base_quality_stats_per_tile_reference_.at(template_segment).at(ref_base) += base_quality_stats_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error);
				error_rate_for_position_per_tile_reference_.at(template_segment).at(ref_base) += error_rate_for_position_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error);
				base_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base) += base_quality_for_error_rate_per_tile_per_error_reference_.at(template_segment).at(ref_base).at(dom_error);
			}

			ShrinkVect(base_quality_stats_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(error_rate_for_position_per_tile_reference_.at(template_segment).at(ref_base));
			ShrinkVect(base_quality_for_error_rate_per_tile_reference_.at(template_segment).at(ref_base));
		}
	}
}

void QualityStats::PreparePlotting(){
	SplitPairedSequenceQuality();
	SumTiles();
	CalculateQualityStats();
}

void QualityStats::PrepareTesting(){
	if(sequence_quality_mean_for_gc_per_tile_reference_.at(0).size()){
		// Only if there are accepted reads, we do something
		PrepareEstimation();
		PreparePlotting();

		for( uintTempSeq template_segment=2; template_segment--; ){
			for( uintBaseCall nucleotide=5; nucleotide--; ){
				nucleotide_quality_.at(template_segment).at(nucleotide).Calculate();
			}
		}
	}
}
