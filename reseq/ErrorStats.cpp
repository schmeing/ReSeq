#include "ErrorStats.h"
using reseq::ErrorStats;

#include <algorithm>
using std::max;
using std::min;

#include "reportingUtils.hpp"

//include "utlities.hpp"
using reseq::utilities::getConst;

void ErrorStats::Prepare(uintTileId num_tiles, uintQual size_qual, uintReadLen size_pos, uintReadLen size_indel){
	for( auto template_segment=2; template_segment--; ){
		for( auto ref_base=4; ref_base--; ){
			for( auto dom_error=5; dom_error--; ){
				SetDimensions( tmp_called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, 5, size_qual );
				SetDimensions( tmp_called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, 5, size_pos );
				SetDimensions( tmp_called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, 5, size_pos );
				SetDimensions( tmp_called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, 5, 101 );

				SetDimensions( tmp_error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, size_pos, size_qual );
				SetDimensions( tmp_error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, size_pos, size_pos );
				SetDimensions( tmp_error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error), num_tiles, size_pos, 101 );
			}

			for( auto prev_base=5; prev_base--; ){
				for( auto prev_called_base=6; prev_called_base--; ){
					tmp_called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(prev_base).at(prev_called_base).resize(size_qual);
				}
			}
		}

		tmp_errors_per_read_.at(template_segment).resize(size_pos);
	}

	for( auto indel_type=2; indel_type--; ){
		for( auto last_call=6; last_call--; ){
			SetDimensions( tmp_indel_by_indel_pos_.at(indel_type).at(last_call), size_indel, 7 );
			SetDimensions( tmp_indel_by_position_.at(indel_type).at(last_call), size_pos, 7 );
			SetDimensions( tmp_indel_by_gc_.at(indel_type).at(last_call), 101, 7 );
			SetDimensions( tmp_indel_pos_by_position_.at(indel_type).at(last_call), size_pos, size_indel );
			SetDimensions( tmp_indel_pos_by_gc_.at(indel_type).at(last_call), 101, size_indel );
			SetDimensions( tmp_gc_by_position_.at(indel_type).at(last_call), size_pos, 101 );
		}
	}
}

void ErrorStats::Finalize(){
	for( auto template_segment=2; template_segment--; ){
		for( auto ref_base=4; ref_base--; ){
			for( auto dom_error=5; dom_error--; ){
				called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
				called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
				called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
				called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error) );

				error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
				error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
				error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error).Acquire( tmp_error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error) );
			}

			for( auto prev_base=5; prev_base--; ){
				for( auto prev_called_base=6; prev_called_base--; ){
					called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(prev_base).at(prev_called_base).Acquire( tmp_called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(prev_base).at(prev_called_base) );
				}
			}
		}

		errors_per_read_.at(template_segment).Acquire( tmp_errors_per_read_.at(template_segment) );
	}

	for( auto indel_type=2; indel_type--; ){
		for( auto last_call=6; last_call--; ){
			indel_by_indel_pos_.at(indel_type).at(last_call).Acquire( tmp_indel_by_indel_pos_.at(indel_type).at(last_call) );
			indel_by_position_.at(indel_type).at(last_call).Acquire( tmp_indel_by_position_.at(indel_type).at(last_call) );
			indel_by_gc_.at(indel_type).at(last_call).Acquire( tmp_indel_by_gc_.at(indel_type).at(last_call) );
			indel_pos_by_position_.at(indel_type).at(last_call).Acquire( tmp_indel_pos_by_position_.at(indel_type).at(last_call) );
			indel_pos_by_gc_.at(indel_type).at(last_call).Acquire( tmp_indel_pos_by_gc_.at(indel_type).at(last_call) );
			gc_by_position_.at(indel_type).at(last_call).Acquire( tmp_gc_by_position_.at(indel_type).at(last_call) );
		}
	}
}

void ErrorStats::Shrink(){
	for( uintTempSeq template_segment=2; template_segment--; ){
		for( auto ref_base = 4; ref_base--; ){
			for( auto dom_error = 5; dom_error--; ){
				ShrinkVect(called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(called_bases_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(called_bases_by_error_num_per_tile_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(called_bases_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error));

				ShrinkVect(error_num_by_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(error_num_by_position_per_tile_.at(template_segment).at(ref_base).at(dom_error));
				ShrinkVect(error_num_by_error_rate_per_tile_.at(template_segment).at(ref_base).at(dom_error));
			}

			for(auto called_base = 5; called_base--; ){
				for(auto prev_base = 6; prev_base--; ){ // Previously called base here, therefore size of 6 and not 5 like before
					ShrinkVect(called_bases_by_base_quality_per_previous_called_base_.at(template_segment).at(ref_base).at(called_base).at(prev_base));
				}
			}
		}
	}

	for( auto indel_type=2; indel_type--; ){
		for( auto last_call=6; last_call--; ){
			ShrinkVect( indel_by_indel_pos_.at(indel_type).at(last_call) );
			ShrinkVect( indel_by_position_.at(indel_type).at(last_call) );
			ShrinkVect( indel_by_gc_.at(indel_type).at(last_call) );
			ShrinkVect( indel_pos_by_position_.at(indel_type).at(last_call) );
			ShrinkVect( indel_pos_by_gc_.at(indel_type).at(last_call) );
			ShrinkVect( gc_by_position_.at(indel_type).at(last_call) );
		}
	}
}

void ErrorStats::PreparePlotting(){
	for( uintTempSeq template_segment=2; template_segment--; ){
		for(auto ref_base = called_bases_by_base_quality_.at(template_segment).size(); ref_base--; ){
			for(auto called_base = called_bases_by_base_quality_.at(template_segment).at(ref_base).size(); called_base--; ){
				// Clear vectors
				called_bases_by_base_quality_.at(template_segment).at(ref_base).at(called_base).Clear();
				called_bases_by_position_.at(template_segment).at(ref_base).at(called_base).Clear();

				// Fill vectors
				for(auto dom_error = called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).size(); dom_error--; ){
					for(auto tile_id = called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).to(); tile_id-- > called_bases_by_base_quality_per_tile_.at(template_segment).at(ref_base).at(dom_error).from(); ){
						called_bases_by_base_quality_.at(template_segment).at(ref_base).at(called_base) += getConst(called_bases_by_base_quality_per_tile_).at(template_segment).at(ref_base).at(dom_error).at(tile_id)[called_base]; // called_base might not exist in all tiles for all dom_errors (make const to not change object)
						called_bases_by_position_.at(template_segment).at(ref_base).at(called_base) += getConst(called_bases_by_position_per_tile_).at(template_segment).at(ref_base).at(dom_error).at(tile_id)[called_base]; // called_base might not exist in all tiles for all dom_errors (make const to not change object)
					}
				}

				// Shrink vectors
				ShrinkVect(called_bases_by_base_quality_.at(template_segment).at(ref_base).at(called_base));
				ShrinkVect(called_bases_by_position_.at(template_segment).at(ref_base).at(called_base));
			}
		}
	}

	// Prepare indel vectors
	// Clear vectors
	for( uintInDelType type=2; type--; ){
		indel_error_by_length_.at(type).Clear();
		indel_error_by_position_.at(type).Clear();
		indel_error_by_gc_.at(type).Clear();
	}

	for( uintInDelType prev_indel=2; prev_indel--; ){
		// Fill vectors
		for( uintBaseCall prev_call=6; prev_call--; ){
			for(auto indel_pos=indel_by_indel_pos_.at(prev_indel).at(prev_call).to(); indel_pos-- > indel_by_indel_pos_.at(prev_indel).at(prev_call).from(); ){ // indel_pos=0 does not exist for prev_indel=1(deletion) as no indel is treated as prev_indel=0(insertion)
				// Deletions
				indel_error_by_length_.at(1)[indel_pos+1] += getConst(indel_by_indel_pos_).at(prev_indel).at(prev_call).at(indel_pos)[1]; // There might be no deletions (make const to not change object)

				// Insertions
				for( uintBaseCall ins_nuc=max(static_cast<size_t>(2), indel_by_indel_pos_.at(prev_indel).at(prev_call).at(indel_pos).from()); ins_nuc<indel_by_indel_pos_.at(prev_indel).at(prev_call).at(indel_pos).to(); ++ins_nuc){
					indel_error_by_length_.at(0)[indel_pos+1] += indel_by_indel_pos_.at(prev_indel).at(prev_call).at(indel_pos).at(ins_nuc);
				}
			}

			for(auto read_pos=indel_by_position_.at(prev_indel).at(prev_call).to(); read_pos-- > indel_by_position_.at(prev_indel).at(prev_call).from(); ){
				// Deletions
				indel_error_by_position_.at(1)[read_pos] += getConst(indel_by_position_).at(prev_indel).at(prev_call).at(read_pos)[1]; // There might be no deletions (make const to not change object)
				// Insertions
				for( uintBaseCall ins_nuc=max(static_cast<size_t>(2), indel_by_position_.at(prev_indel).at(prev_call).at(read_pos).from()); ins_nuc<indel_by_position_.at(prev_indel).at(prev_call).at(read_pos).to(); ++ins_nuc){
					indel_error_by_position_.at(0)[read_pos] += indel_by_position_.at(prev_indel).at(prev_call).at(read_pos).at(ins_nuc);
				}
			}

			for(auto gc=indel_by_gc_.at(prev_indel).at(prev_call).to(); gc-- > indel_by_gc_.at(prev_indel).at(prev_call).from(); ){
				// Deletions
				indel_error_by_gc_.at(1)[gc] += getConst(indel_by_gc_).at(prev_indel).at(prev_call).at(gc)[1]; // There might be no deletions (make const to not change object)
				// Insertions
				for( uintBaseCall ins_nuc=max(static_cast<size_t>(2), indel_by_gc_.at(prev_indel).at(prev_call).at(gc).from()); ins_nuc<indel_by_gc_.at(prev_indel).at(prev_call).at(gc).to(); ++ins_nuc){
					indel_error_by_gc_.at(0)[gc] += indel_by_gc_.at(prev_indel).at(prev_call).at(gc).at(ins_nuc);
				}
			}

		}
	}

	for( uintInDelType type=2; type--; ){
		// Transform vector with indel bases by length up to that base into vector with indels by total indel length
		// Start from the back and remove the counts of longer length from the ones with shorter length
		if(indel_error_by_length_.at(type).size()){
			uintNucCount excess_counts(0);
			for(auto length=indel_error_by_length_.at(type).to(); --length; ){
				indel_error_by_length_.at(type).at(length) -= excess_counts;
				excess_counts += indel_error_by_length_.at(type).at(length);
			}
		}

		// Shrink vectors
		indel_error_by_length_.at(type).Shrink();
		indel_error_by_position_.at(type).Shrink();
		indel_error_by_gc_.at(type).Shrink();
	}
}

void ErrorStats::PrepareSimulation(){
	max_len_deletion_ = 0;
	for( auto &vect : indel_by_indel_pos_.at(1) ){
		if( vect.to() > max_len_deletion_ ){
			max_len_deletion_ = vect.to();
		}
	}
}
