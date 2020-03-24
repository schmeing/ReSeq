#include "Simulator.h"
using reseq::Simulator;

#include <algorithm>
using std::max;
using std::min;
#include <array>
using std::array;
//include <atomic>
using std::atomic;
#include <cmath>
using std::round;
#include <exception>
using std::exception;
#include <iomanip>
using std::setprecision;
#include <iostream>
#include <limits>
using std::numeric_limits;
#include <map>
using std::map;
//include <mutex>
using std::mutex;
using std::lock_guard;
//include <random>
using std::discrete_distribution;
using std::mt19937_64;
//include <set>
using std::set;
#include <sstream>
using std::stringstream;
#include <stdexcept>
using std::runtime_error;
//include <string>
using std::string;
#include <thread>
using std::thread;
//include <utility>
using std::pair;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/seq_io.h>
using seqan::appendValue;
using seqan::CharString;
using seqan::Dna;
using seqan::DnaString;
using seqan::Dna5;
using seqan::Dna5String;
using seqan::Exact;
using seqan::Exception;
using seqan::resizeSpace;
using seqan::StringSet;
//include <seqan/bam_io.h>
using seqan::atEnd;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamHeader; // <seqan/bam_io/bam_header_record.h>
using seqan::CigarElement;
using seqan::FormattedFileContext;
using seqan::readRecord; // <>
using seqan::readHeader; // <>

//include "utilities.hpp"
using reseq::utilities::ConstDna5StringReverseComplement;
using reseq::utilities::CigarString;
using reseq::utilities::Complement;
using reseq::utilities::at;
using reseq::utilities::CreateDir;
using reseq::utilities::DeleteFile;
using reseq::utilities::Divide;
using reseq::utilities::DominantBase;
using reseq::utilities::IsGC;
using reseq::utilities::Percent;
using reseq::utilities::ReverseComplementorDna;
using reseq::utilities::SafePercent;
using reseq::utilities::SetToMax;
using reseq::utilities::TransformDistanceToStartOfErrorRegion;
using reseq::utilities::TrueRandom;

double Simulator::CoveragePropLostFromAdapters(const DataStats &stats){
	uintRefLenCalc adapter_bases(0), total_bases(0);
	for(uintTempSeq seg=2; seg--; ){
		for(uintSeqLen frag_len=stats.ReadLengthsByFragmentLength(seg).from(); frag_len < stats.ReadLengthsByFragmentLength(seg).to(); ++frag_len){
			for(uintReadLen read_len=stats.ReadLengthsByFragmentLength(seg).at(frag_len).from(); read_len < stats.ReadLengthsByFragmentLength(seg).at(frag_len).to(); ++read_len){
				total_bases += stats.ReadLengthsByFragmentLength(seg).at(frag_len).at(read_len) * read_len;
				if(frag_len < read_len){
					adapter_bases += (stats.ReadLengthsByFragmentLength(seg).at(frag_len).at(read_len)-stats.NonMappedReadLengthsByFragmentLength(seg)[frag_len][read_len]) * (read_len-frag_len);
					adapter_bases += stats.NonMappedReadLengthsByFragmentLength(seg)[frag_len][read_len] * read_len;
				}
			}
		}
	}
	return static_cast<double>(adapter_bases)/total_bases; // Compensate the part that is made up of adapters and therefore does not produce reference coverage
}

reseq::uintFragCount Simulator::CoverageToNumberPairs(double coverage, uintRefLenCalc total_ref_size, double average_read_length, double adapter_part){
	return round(coverage * total_ref_size / average_read_length / 2 / (1-adapter_part));
}

double Simulator::NumberPairsToCoverage(uintFragCount total_pairs, uintRefLenCalc total_ref_size, double average_read_length, double adapter_part){
	return static_cast<double>(total_pairs) / total_ref_size * average_read_length * 2 * (1-adapter_part);
}

bool Simulator::Flush(){
	if( flush_mutex_.try_lock() ){
		array<StringSet<CharString> *, 2> old_output_ids;
		array<StringSet<Dna5String> *, 2> old_output_seqs;
		array<StringSet<CharString> *, 2> old_output_quals;

		output_mutex_.lock();
		for( int i=2; i--; ){
			old_output_ids.at(i) = output_ids_.at(i);
			output_ids_.at(i) = new StringSet<CharString>;
			reserve(*output_ids_.at(i), kBatchSize, Exact());

			old_output_seqs.at(i) = output_seqs_.at(i);
			output_seqs_.at(i) = new StringSet<Dna5String>;
			reserve(*output_seqs_.at(i), kBatchSize, Exact());

			old_output_quals.at(i) = output_quals_.at(i);
			output_quals_.at(i) = new StringSet<CharString>;
			reserve(*output_quals_.at(i), kBatchSize, Exact());
		}
		output_mutex_.unlock();

		for( int i=2; i--; ){
			try{
				writeRecords(dest_.at(i), *old_output_ids.at(i), *old_output_seqs.at(i), *old_output_quals.at(i));
			}
			catch(const Exception &e){
				print_mutex_.lock();
				printErr << "Could not write records " << written_records_+1 << " to " << (written_records_+length(*old_output_ids.at(i))) << ": " << e.what() << std::endl;
				print_mutex_.unlock();
				return false;
			}
		}

		written_records_ += length(*old_output_ids.at(0));
		print_mutex_.lock();
		printInfo << "Generated " << written_records_ << " read pairs." << std::endl;
		print_mutex_.unlock();

		for( int i=2; i--; ){
			delete old_output_ids.at(i);
			delete old_output_seqs.at(i);
			delete old_output_quals.at(i);
		}

		flush_mutex_.unlock();
	}

	return true;
}

bool Simulator::Output(const SimPair &sim_reads){
	output_mutex_.lock();
	for(uintTempSeq template_segment=2; template_segment--; ){
		appendValue(*output_ids_.at(template_segment), sim_reads.id_.at(template_segment));
		appendValue(*output_seqs_.at(template_segment), sim_reads.seq_.at(template_segment));
		appendValue(*output_quals_.at(template_segment), sim_reads.qual_.at(template_segment));
	}

	if(length(*output_ids_.at(0)) >= kBatchSize){
		output_mutex_.unlock();
		return this->Flush();
	}
	else{
		output_mutex_.unlock();
		return true;
	}
}

inline void Simulator::IncrementBlockPos(uintSeqLen &block_pos, const SimBlock* &block, intVariantId &cur_var){
	if( block->sys_errors_.size() <= ++block_pos ){
		block = block->next_block_;
		block_pos = 0;
		cur_var = 0;
	}
}

void Simulator::GetSysErrorFromBlock(Dna5 &dom_error, uintPercent &error_rate, const SimBlock* &block, uintSeqLen &block_pos, intVariantId &cur_var, uintSeqLen &var_pos, uintAlleleId allele){
	bool no_variant(true);

	if(var_pos){
		no_variant = false;

		dom_error = block->sys_errors_.at(block_pos).first;
		error_rate = block->sys_errors_.at(block_pos).second;

		if(++var_pos >= block->err_variants_.at(cur_var).var_errors_.size()){
			var_pos = 0;
			IncrementBlockPos(block_pos, block, ++cur_var);
		}
	}
	else{
		while(cur_var < block->err_variants_.size() && block->err_variants_.at(cur_var).position_ <= block_pos){
			if(block->err_variants_.at(cur_var).InAllele(allele)){
				if(0 == block->err_variants_.at(cur_var).var_errors_.size()){
					// Deletion
					IncrementBlockPos(block_pos, block, ++cur_var);
				}
				else{
					no_variant = false;

					dom_error = block->err_variants_.at(cur_var).var_errors_.at(0).first;
					error_rate = block->err_variants_.at(cur_var).var_errors_.at(0).second;

					if(1 == block->err_variants_.at(cur_var).var_errors_.size()){
						// Substitution
						IncrementBlockPos(block_pos, block, ++cur_var);
						++cur_var;
					}
					else{
						// Insertion
						var_pos = 1;
					}

					break;
				}
			}
			else{
				++cur_var;
			}
		}
	}

	if(no_variant){
		dom_error = block->sys_errors_.at(block_pos).first;
		error_rate = block->sys_errors_.at(block_pos).second;

		IncrementBlockPos(block_pos, block, cur_var);
	}
}

bool Simulator::FillReadPart(
		SimPair &sim_reads,
		uintTempSeq template_segment,
		uintTileId tile_id,
		const DnaString &org_sequence,
		uintReadLen org_pos,
		char base_cigar_element,
		const SimBlock* block,
		uintSeqLen block_pos,
		uintAdapterId adapter_id,
		ReadFillParameter &par,
		pair<intVariantId, uintSeqLen> first_variant,
		uintAlleleId allele,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen){
	uintReadLen cigar_element_length(0);
	char cigar_element(base_cigar_element);

	ErrorStats::InDelDef indel;
	double prob_sum;
	vector<double > prob_indel, prob_quality, prob_base_call;
	Dna5 dom_error;

	auto var_pos = first_variant.second;
	auto cur_var = first_variant.first;

	while(par.read_pos_ < par.read_length_ && org_pos < length(org_sequence)){
		// Still in the original fragment
		indel = static_cast<ErrorStats::InDelDef>(estimates.InDels(par.previous_indel_type_, par.base_call_).Draw(prob_indel, prob_sum, {par.indel_pos_, par.read_pos_, par.gc_seq_}, rdist.ZeroToOne(rgen)));
		if( 0.0 == prob_sum ){
			indel = ErrorStats::kNoInDel;
		}
		
		if(ErrorStats::kNoInDel == indel){
			// Load systematic error
			if( block ){
				GetSysErrorFromBlock(dom_error, par.error_rate_, block, block_pos, cur_var, var_pos, allele);
			}
			else{
				dom_error = adapter_sys_error_.at(template_segment).at(adapter_id).at(org_pos).first;
				par.error_rate_ = adapter_sys_error_.at(template_segment).at(adapter_id).at(org_pos).second;
			}

			// Roll quality
			par.qual_ = estimates.Quality(template_segment, tile_id, at(org_sequence, org_pos)).Draw(prob_quality, prob_sum, {par.seq_qual_, par.qual_, par.read_pos_, par.error_rate_}, rdist.ZeroToOne(rgen));
			if( 0.0 == prob_sum ){
				if(par.read_pos_){
					par.qual_ = at(sim_reads.qual_.at(template_segment), par.read_pos_-1) - stats.PhredQualityOffset(); // Take last quality again
				}
				else{
					par.qual_ = estimates.Quality(template_segment, tile_id, at(org_sequence, org_pos)).MaxValue();
				}
			}
			at(sim_reads.qual_.at(template_segment), par.read_pos_) = par.qual_ + stats.PhredQualityOffset();

			// Roll base call
			par.base_call_ = estimates.BaseCall(template_segment, tile_id, at(org_sequence, org_pos), dom_error).Draw(prob_base_call, prob_sum, {par.qual_, par.read_pos_, par.num_errors_, par.error_rate_}, rdist.ZeroToOne(rgen));
			if( 0.0 == prob_sum ){
				par.base_call_ = at(org_sequence, org_pos);
			}
			at(sim_reads.seq_.at(template_segment), par.read_pos_) = par.base_call_;

			// Update cigar
			if(base_cigar_element == cigar_element){
				++cigar_element_length;
			}
			else{
				append( sim_reads.cigar_.at(template_segment), CigarElement<>(cigar_element, cigar_element_length ) );
				cigar_element = base_cigar_element;
				cigar_element_length = 1;
				par.indel_pos_ = 0;
				par.previous_indel_type_ = 0;
			}

			// Update num_errors_
			if(par.base_call_ != at(org_sequence, org_pos)){
				++par.num_errors_;
			}

			// Update position
			++par.read_pos_;
			++org_pos;
		}
		else if(ErrorStats::kDeletion == indel){
			// Load systematic error
			if( block ){
				par.error_rate_ = block->sys_errors_.at(block_pos).second;

				if( var_pos && ++var_pos >= block->err_variants_.at(cur_var).var_errors_.size() ){
					var_pos = 0;
				}
				if( 0 == var_pos && block->sys_errors_.size() <= ++block_pos ){
					block = block->next_block_;
					block_pos = 0;
					cur_var = 0;
				}
			}
			else{
				par.error_rate_ = adapter_sys_error_.at(template_segment).at(adapter_id).at(org_pos).second;
			}

			// Update cigar
			if('D' == cigar_element){
				++cigar_element_length;
				++par.indel_pos_;
			}
			else{
				append( sim_reads.cigar_.at(template_segment), CigarElement<>(cigar_element, cigar_element_length ) );
				cigar_element = 'D';
				cigar_element_length = 1;
				par.indel_pos_ = 1;
				par.previous_indel_type_ = 1;
			}

			// Update num_errors_
			++par.num_errors_;

			// Update position
			++org_pos;
		}
		else{
			// Insertion
			// Roll quality
			at(sim_reads.qual_.at(template_segment), par.read_pos_) = stats.PhredQualityOffset() + estimates.Quality(template_segment, tile_id, at(org_sequence, org_pos)).Draw(prob_quality, prob_sum, {par.seq_qual_, par.qual_, par.read_pos_, par.error_rate_}, rdist.ZeroToOne(rgen));
			if( 0.0 == prob_sum ){
				at(sim_reads.qual_.at(template_segment), par.read_pos_) = stats.PhredQualityOffset() + par.qual_; // Take last quality again
			}
			
			// Take base call from indel
			at(sim_reads.seq_.at(template_segment), par.read_pos_) = indel-2; // This base call is not stored in par.base_call, because the normal calls have more predictive power for indel calls

			// Update cigar
			if('I' == cigar_element){
				++cigar_element_length;
				++par.indel_pos_;
			}
			else{
				append( sim_reads.cigar_.at(template_segment), CigarElement<>(cigar_element, cigar_element_length ) );
				cigar_element = 'I';
				cigar_element_length = 1;
				par.indel_pos_ = 1;
				par.previous_indel_type_ = 0;
			}

			// Update position
			++par.num_errors_;
			++par.read_pos_;
		}
	}

	// Append final cigar element
	if(cigar_element_length){
		append( sim_reads.cigar_.at(template_segment), CigarElement<>(cigar_element, cigar_element_length ) );
	}

	return true;
}

bool Simulator::FillRead(
		SimPair &sim_reads,
		uintTempSeq template_segment,
		uintTileId tile_id,
		uintSeqLen fragment_length,
		const SimBlock* start_block,
		uintSeqLen block_start_pos,
		pair<intVariantId, uintSeqLen> first_variant,
		uintAlleleId allele,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen){
	ReadFillParameter par;

	par.read_length_ = rdist.ReadLength( template_segment, fragment_length, rgen );

	resize( sim_reads.seq_.at(template_segment), par.read_length_ );
	resize( sim_reads.qual_.at(template_segment), par.read_length_ );

	clear(sim_reads.cigar_.at(template_segment));

	uintAdapterId adapter_id(0);

	// Get gc content and mean error rate of read
	uintReadLen seq_length(min(par.read_length_, static_cast<uintReadLen>(length(sim_reads.org_seq_.at(template_segment))))), read_pos( seq_length );
	uintReadLenCalc mean_error_rate(0);
	if( seq_length ){
		uintSeqLen block_pos(block_start_pos);
		const SimBlock *block(start_block);

		Dna5 dom_error;
		uintPercent error_rate;
		auto var_pos = first_variant.second;
		auto cur_var = first_variant.first;

		for( ; read_pos--; ){
			// gc
			if( IsGC(sim_reads.org_seq_.at(template_segment), read_pos) ){
				++par.gc_seq_;
			}

			// error rate
			GetSysErrorFromBlock(dom_error, error_rate, block, block_pos, cur_var, var_pos, allele);
			mean_error_rate += error_rate;
		}

		par.gc_seq_ = Percent(par.gc_seq_, seq_length);
		mean_error_rate = Divide(mean_error_rate, static_cast<uintReadLenCalc>(seq_length));
	}
	else{
		// Adapter only read
		adapter_id = rdist.Adapter(template_segment,rgen);
		read_pos = length(stats.Adapters().Sequence(template_segment, adapter_id));

		for( ; read_pos--; ){
			// gc
			if( IsGC(stats.Adapters().Sequence(template_segment, adapter_id), read_pos) ){
				++par.gc_seq_;
			}

			// error rate
			mean_error_rate += adapter_sys_error_.at(template_segment).at(adapter_id).at(read_pos).second;
		}

		par.gc_seq_ = Percent(par.gc_seq_, static_cast<uintReadLen>(length(stats.Adapters().Sequence(template_segment, adapter_id))));
		mean_error_rate = Divide(mean_error_rate, static_cast<uintReadLenCalc>(length(stats.Adapters().Sequence(template_segment, adapter_id))));
	}

	// Roll sequence quality
	vector<double> prob_seq_qual;
	double prob_sum;
	par.seq_qual_ = estimates.SequenceQuality(template_segment, tile_id).Draw(prob_seq_qual, prob_sum, {par.gc_seq_, mean_error_rate, fragment_length/QualityStats::kSqFragmentLengthBinSize}, rdist.ZeroToOne(rgen));
	if( 0.0 == prob_sum ){
		par.seq_qual_ = estimates.SequenceQuality(template_segment, tile_id).MostLikely(); // Take the one with the highest general probability, as parameter specific probabilities do not exist
	}
	
	// Reference part of the sequence
	if( !FillReadPart(sim_reads, template_segment, tile_id, sim_reads.org_seq_.at(template_segment), 0, 'M', start_block, block_start_pos, 0, par, first_variant, allele, stats, estimates, rdist, rgen) ){
		return false;
	}

	if(par.read_pos_ < par.read_length_){
		// Adapter needs to be added as we reached the fragment end
		if(0 == adapter_id){
			adapter_id = rdist.Adapter(template_segment,rgen);
		}

		uintReadLen adapter_pos(0);
		if(0 == par.read_pos_){
			// Read starts with the adapter: Potentially remove some bases at the beginning of the adapter
			discrete_distribution<uint64_t> rdist_cut( stats.Adapters().StartCut(template_segment, adapter_id).begin(), stats.Adapters().StartCut(template_segment, adapter_id).end() );
			adapter_pos = rdist_cut(rgen) + stats.Adapters().StartCut(template_segment, adapter_id).from();
		}

		// Adapter part of the sequence
		if( !FillReadPart(sim_reads, template_segment, tile_id, stats.Adapters().Sequence(template_segment, adapter_id), adapter_pos, 'S', NULL, 0, adapter_id, par, {0,0}, 0, stats, estimates, rdist, rgen) ){
			return false;
		}

		// Rubbish part of the sequence
		if(par.read_pos_ < par.read_length_){
			// The sequence is longer than the adapter
			append( sim_reads.cigar_.at(template_segment), CigarElement<>('H', par.read_length_-par.read_pos_) );

			// Adapter does not go until end of read: Add poly-A tail
			double prob_sum;
			vector<double > prob_quality;
			auto tail_length = rdist.PolyATailLength( rgen );
			for(uintReadLen pos_tail=0; pos_tail < tail_length && par.read_pos_ < par.read_length_; ++pos_tail){
				// There is no clear way to set the qualities here, so just do something
				par.qual_ = estimates.Quality(template_segment, tile_id, 0).Draw(prob_quality, prob_sum, {par.seq_qual_, par.qual_, par.read_pos_, par.error_rate_}, rdist.ZeroToOne(rgen));
				if( 0.0 == prob_sum ){
					par.qual_ = at(sim_reads.qual_.at(template_segment), par.read_pos_-1) - stats.PhredQualityOffset(); // Take last quality again
				}
				at(sim_reads.qual_.at(template_segment), par.read_pos_) = par.qual_ + stats.PhredQualityOffset();

				// Set called base to A
				at(sim_reads.seq_.at(template_segment), par.read_pos_++) = 0;
			}

			// After the tail add random bases
			while(par.read_pos_ < par.read_length_){
				// There is no clear way to set the qualities here, so just do something
				par.qual_ = estimates.Quality(template_segment, tile_id, 0).Draw(prob_quality, prob_sum, {par.seq_qual_, par.qual_, par.read_pos_, par.error_rate_}, rdist.ZeroToOne(rgen));
				if( 0.0 == prob_sum ){
					par.qual_ = at(sim_reads.qual_.at(template_segment), par.read_pos_-1) - stats.PhredQualityOffset(); // Take last quality again
				}
				at(sim_reads.qual_.at(template_segment), par.read_pos_) = par.qual_ + stats.PhredQualityOffset();

				// Set random called base to A
				at(sim_reads.seq_.at(template_segment), par.read_pos_++) = rdist.OverrunBases(rgen);
			}
		}
	}

	return true;
}

void Simulator::CreateReadId(
		CharString &id,
		uintRefSeqBin block_number,
		uintFragCount read_number,
		uintSeqLen start_pos,
		uintSeqLen end_pos,
		uintAlleleId allele,
		uintTile tile,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		const CigarString &cigar){
	// start_pos should start at 1 (not at 0 as the array)[0 means adapter only], end_pos should be the last included position (not the first not included position like usual)
	stringstream readid_stream;
	readid_stream << record_base_identifier_ << block_number << '_' << read_number;

	if(start_pos && 1 < ref.NumAlleles()){
		readid_stream << "_allele" << allele;
	}

	readid_stream << ':' << start_pos << ':';
	if(start_pos){
		readid_stream << ref.ReferenceIdFirstPart(ref_seq_id);
	}
	else{
		readid_stream << "Adapter";
	}
	readid_stream << ':' << end_pos << ':' << tile << ":1337:1337 ";

	for( const auto &cigar_element : cigar ){
		readid_stream << cigar_element.count << cigar_element.operation;
	}

	id = readid_stream.str().c_str();
}

bool Simulator::CreateReads(
		const Reference &ref,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		SimPair &sim_reads,
		mt19937_64 &rgen,
		uintFragCount counts,
		bool strand,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		uintSeqLen fragment_length,
		uintFragCount &read_number,
		const SimBlock *start_block,
		uintSeqLen start_position_forward,
		uintSeqLen end_position_forward,
		pair<intVariantId, uintSeqLen> start_variant,
		pair<intVariantId, uintSeqLen> end_variant
		){
		// Find block corresponding to fragment end
		const SimBlock *end_block = start_block;
		if(start_block){
			while(end_block->start_pos_+end_block->sys_errors_.size() < end_position_forward){
				end_block = end_block->next_block_;
			}
			end_block = end_block->partner_block_; // Second read is reversed, so we need the reversed block
		}

		// Define values for this reference position
		uintSeqLen print_start_position, print_end_position;
		uintRefSeqBin block_id;
		array<const SimBlock *, 2> block;
		array<uintSeqLen, 2> block_start_pos;
		array<pair<intVariantId, uintSeqLen>, 2> variant = {{{0,0},{0,0}}};
	if(fragment_length){
		block_id = start_block->id_;

		if(strand){
			print_start_position = end_position_forward;
			print_end_position = start_position_forward+1; // +1 because counting starts from 1 there
		}
		else{
			print_start_position = start_position_forward+1; // +1 because counting starts from 1 there
			print_end_position = end_position_forward;
		}

		block.at(strand) = start_block;
		block.at(!strand) = end_block; // First read is reversed and at fragment end for reverse strand

		block_start_pos.at(strand) = start_position_forward - start_block->start_pos_;
		block_start_pos.at(!strand) = end_block->start_pos_ + end_block->sys_errors_.size() - end_position_forward;

		if(ref.VariantsLoaded()){
			variant.at(strand) = {start_variant.first-start_block->first_variant_id_, start_variant.second};
			variant.at(!strand) = {end_block->first_variant_id_-end_variant.first, (end_variant.second?length(ref.Variants(ref_seq_id).at(end_variant.first).var_seq_)-end_variant.second:0)};
		}
	}
	else{
		// Adapter only pair
		block_id = 0;
		print_start_position = 0;
		print_end_position = 0;
		block.fill(NULL);
		block_start_pos.fill(0);
	}

	// Simulate
	for( auto n = counts; n--; ){
		++read_number;

		auto tile_id = rdist.TileId( rgen );

		for(uintTempSeq template_segment=2; template_segment--; ){
			if(!FillRead(sim_reads, template_segment, tile_id, fragment_length, block.at(template_segment), block_start_pos.at(template_segment), variant.at(template_segment), allele, stats, estimates, rdist, rgen)){
				return false;
			}

			CreateReadId(sim_reads.id_.at(template_segment), block_id, read_number, print_start_position, print_end_position, allele, stats.Tiles().Tiles().at(tile_id), ref_seq_id, ref, sim_reads.cigar_.at(template_segment));
		}

		this->Output( sim_reads );
	}

	return true;
}

void Simulator::ResetSystematicErrorCounters(const Reference &ref){
	sys_last_base_ = 4;
	sys_dom_base_.Clear();

	if(ref.VariantsLoaded()){
		sys_last_var_pos_per_allele_.clear();
		sys_last_var_pos_per_allele_.resize(ref.NumAlleles(), numeric_limits<uintSeqLen>::max());
		sys_last_base_per_allele_.clear();
		sys_last_base_per_allele_.resize(ref.NumAlleles(), 4);

		for(auto allele = ref.NumAlleles(); allele--; ){
			sys_dom_base_per_allele_.at(allele).Clear();
		}
	}

	sys_gc_ = 0;
	sys_gc_bases_ = 0;

	distance_to_start_of_error_region_ = 0;
	start_error_rate_ = 0;
}

bool Simulator::LoadSysErrorRecord(uintRefSeqId ref_id, const Reference &ref){
	try{
		readRecord(sys_id_, sys_dom_err_, sys_err_perc_, sys_file_);
	}
	catch(const Exception &e){
		printErr << "Could not read systematic error profile for reference sequence '" << ref.ReferenceId(ref_id) << "': " << e.what() << std::endl;
		simulation_error_ = true;
		return false;
	}

	if(length(sys_dom_err_) != ref.SequenceLength(ref_id)){
		printErr << "Systematic error profile '" << sys_id_ << "' (length " << length(sys_dom_err_) << ") does not match reference sequence '" << ref.ReferenceId(ref_id) << "' (length " << ref.SequenceLength(ref_id) << "). Wrong file or order incorrect?"<< std::endl;
		simulation_error_ = true;
		return false;
	}

	return true;
}

void Simulator::SetSystematicErrorVariantsReverse(uintSeqLen &start_dist_error_region, uintPercent &start_rate, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	if(ref.VariantsLoaded()){
		uintAlleleId chosen_allele;
		vector<pair<Dna5,uintPercent>> tmp_sys_errors;
		Dna5 base, last_base;
		uintReadLen gc(0), gc_bases(0);

		auto var_id=block.first_variant_id_;
		for(; 0 <= var_id && ref.Variants(ref_seq_id).at(var_id).position_ >= end_pos; --var_id){
			auto &var( ref.Variants(ref_seq_id).at(var_id) );
			chosen_allele = var.FirstAllele();
			tmp_sys_errors.clear();
			uintSeqLen rev_pos = ref.SequenceLength(ref_seq_id)-var.position_-1;
			auto block_pos = block.start_pos_+block.sys_errors_.size()-var.position_-1;

			// Get last_base
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele) == var.position_+1 ){
				last_base = sys_last_base_per_allele_.at(chosen_allele);
			}
			else{
				if(var.position_+1 < ref.SequenceLength(ref_seq_id)){
					last_base = Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), var.position_+1));
				}
				else{
					last_base = 4;
				}
			}

			// Get dom_base
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele) <= var.position_+DominantBase::kLastX ){
				// Variant within DominantBase::kLastX of last variant
				for( auto pos = sys_last_var_pos_per_allele_.at(chosen_allele)-1; pos > var.position_; --pos ){
					sys_dom_base_per_allele_.at(chosen_allele).Update( Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos)) );
				}
			}
			else{
				sys_dom_base_per_allele_.at(chosen_allele).Clear();
				if(rev_pos){
					// If we have bases before the variant add them
					sys_dom_base_per_allele_.at(chosen_allele).Set( ConstDna5StringReverseComplement(ref.ReferenceSequence(ref_seq_id)), rev_pos-1 );
				}
			}

			// Get start_dist_error_region (Variants cannot start an error region)
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < block_pos; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = block.start_pos_+block.sys_errors_.size()-ref.Variants(ref_seq_id).at(var_id+1).position_-1; pos < block_pos; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}

			// Get GC (Variants' GC is not taken into account)
			if(ref.Variants(ref_seq_id).size() > var_id+1 && var_id != block.first_variant_id_ && ref.Variants(ref_seq_id).at(var_id+1).position_ < var.position_ + sys_gc_range_/3 ){ // Updating is ~3 times as much work as recounting
				for(auto pos = ref.SequenceLength(ref_seq_id)-ref.Variants(ref_seq_id).at(var_id+1).position_-1; pos < rev_pos; ++pos){
					UpdateGC(gc, gc_bases, pos, ConstDna5StringReverseComplement(ref.ReferenceSequence(ref_seq_id)));
				}
			}
			else{
				gc_bases = min(rev_pos, static_cast<uintSeqLen>(sys_gc_range_));
				gc = ref.GCContentAbsolut(ref_seq_id, var.position_+1, var.position_+gc_bases+1);
			}

			// Draw systematic errors
			if( 0 < length(var.var_seq_)){
				tmp_sys_errors.reserve(length(var.var_seq_));
				for(uintSeqLen vpos=length(var.var_seq_); vpos--; ){
					base = Complement::Dna(at(var.var_seq_, vpos));
					sys_dom_base_per_allele_.at(chosen_allele).Update(base);
					DrawSystematicError(tmp_sys_errors, base, last_base, sys_dom_base_per_allele_.at(chosen_allele).Get(), SafePercent(gc, gc_bases), TransformDistanceToStartOfErrorRegion(start_dist_error_region), start_rate, estimates);
					last_base = base;
				}
			}

			// Insert variant
			block.err_variants_.emplace_back(block_pos, tmp_sys_errors, var.allele_);

			// Update information for next variants
			auto ref_allele = ref.NumAlleles(); // For this allele we already used the reference, so in case we have another allele identical to the reference up to the current variance we can copy the values
			if(sys_last_var_pos_per_allele_.at(chosen_allele) >= ref.SequenceLength(ref_seq_id) || sys_last_var_pos_per_allele_.at(chosen_allele)+length(var.var_seq_) > var.position_+DominantBase::kLastX){
				ref_allele = chosen_allele;
			}
			for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
				if(var.InAllele(allele)){
					if( allele != chosen_allele){ // Dom base for chosen allele was already updated during the drawing process
						// Update dom base
						if(sys_last_var_pos_per_allele_.at(allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(allele)+length(var.var_seq_) <= var.position_+DominantBase::kLastX){
							for( auto pos = sys_last_var_pos_per_allele_.at(allele)-1; pos > var.position_; --pos ){
								sys_dom_base_per_allele_.at(allele).Update( Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos)) );
							}

							for(uintSeqLen vpos=length(var.var_seq_); vpos--; ){
								sys_dom_base_per_allele_.at(allele).Update( Complement::Dna(at(var.var_seq_, vpos)) );
							}
						}
						else{
							if(ref_allele < ref.NumAlleles()){
								// We can copy from an identical allele
								sys_dom_base_per_allele_.at(allele) = sys_dom_base_per_allele_.at(ref_allele);
							}
							else{
								ref_allele = allele;

								sys_dom_base_per_allele_.at(allele).Clear();
								if(rev_pos){
									// If we have bases before the variant add them
									sys_dom_base_per_allele_.at(allele).Set( ConstDna5StringReverseComplement(ref.ReferenceSequence(ref_seq_id)), rev_pos-1 );
								}

								for(uintSeqLen vpos=length(var.var_seq_); vpos--; ){
									sys_dom_base_per_allele_.at(allele).Update( Complement::Dna(at(var.var_seq_, vpos)) );
								}
							}
						}
					}

					// Update last base
					sys_last_var_pos_per_allele_.at(allele) = var.position_;
					sys_last_base_per_allele_.at(allele) = last_base;
				}
			}
		}

		// Update start_dist_error_region for rest of block as we don't do this while drawing the systematic errors
		if(sys_from_file_){
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = block.start_pos_+block.sys_errors_.size()-ref.Variants(ref_seq_id).at(var_id+1).position_-1; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
		}
	}
}

bool Simulator::CreateUnit(uintRefSeqId ref_id, uintRefSeqBin first_block_id, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&first_reverse_block, SimUnit *&unit){
	unit = new SimUnit(ref_id);
	// at the end of the function block will be filled with the reverse block at the beginning of the reference sequence

	// Create all blocks for the reverse strand of the unit
	auto block = new SimBlock(first_block_id, 0, NULL, block_seed_gen_());
	first_reverse_block = block;
	auto old_block( block );
	while( block->start_pos_ + kBlockSize < ref.SequenceLength(ref_id) ){
		block = new SimBlock(block->id_+1, block->start_pos_ + kBlockSize, NULL, block_seed_gen_());
		block->next_block_ = old_block; // Reverse direction
		old_block->partner_block_ = block; // partner_block_ is previous block for reverse strand
		old_block = block;
	}

	// Load reverse record for next sequence from systematic error file if given
	if(sys_from_file_){
		if( !LoadSysErrorRecord(ref_id, ref) ){
			return false; // Stop as ReadSystematicErrors would crash
		}
	}

	// Additional stuff in case a variant file has been specified
	if(ref.VariantsLoaded()){
		// Make sure the variation for this unit is already loaded
		if( !ref.VariantsLoadedForSequence(ref_id) ){
			lock_guard<mutex> lock(var_read_mutex_);
			if(!ref.ReadVariants(ref_id+1)){ // End is specified, so to get the current one +1
				return false;
			}
		}

		// Assign last variant (first in reverse direction) to each reverse block
		intVariantId first_var = 0;
		auto var_block = first_reverse_block;
		while(var_block->partner_block_){
			while(first_var < ref.Variants(unit->ref_seq_id_).size() && ref.Variants(unit->ref_seq_id_).at(first_var).position_ < var_block->start_pos_ + kBlockSize){
				++first_var;
			}

			var_block->first_variant_id_ = first_var-1;
			var_block = var_block->partner_block_;
		}
		var_block->first_variant_id_ = ref.Variants(unit->ref_seq_id_).size()-1;
	}

	if(ref.MethylationLoaded()){
		if( !ref.MethylationLoadedForSequence(ref_id) ){
			lock_guard<mutex> lock(methylation_read_mutex_);
			if(!ref.ReadMethylation(ref_id+1)){ // End is specified, so to get the current one +1
				return false;
			}
		}

		// Assign last methylation region (first in reverse direction) to each reverse block
		intVariantId first_meth = 0;
		auto meth_block = first_reverse_block;
		while(meth_block->partner_block_){
			while(first_meth < ref.UnmethylatedRegions(unit->ref_seq_id_).size() && ref.UnmethylatedRegions(unit->ref_seq_id_).at(first_meth).first < meth_block->start_pos_ + kBlockSize){
				++first_meth;
			}

			meth_block->first_methylation_id_ = first_meth-1;
			meth_block = meth_block->partner_block_;
		}
		meth_block->first_methylation_id_ = ref.UnmethylatedRegions(unit->ref_seq_id_).size()-1;
	}

	// Set sytematic error for all reverse blocks of the unit
	ConstDna5StringReverseComplement reversed_ref( ref.ReferenceSequence(ref_id) );
	ResetSystematicErrorCounters(ref);
	uintSeqLen start_pos(0), end_pos(ref.SequenceLength(unit->ref_seq_id_) - block->start_pos_);
	while(block){
		block->sys_errors_.reserve(end_pos - start_pos);
		if(sys_from_file_){
			ReadSystematicErrors(block->sys_errors_, start_pos, end_pos);
			SetSystematicErrorVariantsReverse(distance_to_start_of_error_region_, start_error_rate_, *block, unit->ref_seq_id_, block->start_pos_, ref, stats, estimates);
		}
		else{
			auto tmp_distance_to_start_of_error_region = distance_to_start_of_error_region_;
			auto tmp_start_error_rate = start_error_rate_;
			SetSystematicErrors(block->sys_errors_, reversed_ref, start_pos, end_pos, stats, estimates);
			SetSystematicErrorVariantsReverse(tmp_distance_to_start_of_error_region, tmp_start_error_rate, *block, unit->ref_seq_id_, block->start_pos_, ref, stats, estimates);
		}

		start_pos = end_pos;
		end_pos += kBlockSize;
		block = block->next_block_;
	}

	ResetSystematicErrorCounters(ref); // Reset them again for the creation of the forward strand that will follow block by block

	// Load forward record for next sequence from systematic error file if given
	if(sys_from_file_){
		LoadSysErrorRecord(ref_id, ref);
	}

	return true;
}

void Simulator::SetSystematicErrorVariantsForward(uintSeqLen &start_dist_error_region, uintPercent &start_rate, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	if(ref.VariantsLoaded()){
		uintAlleleId chosen_allele;
		vector<pair<Dna5,uintPercent>> tmp_sys_errors;
		Dna5 base, last_base;
		uintReadLen gc(0), gc_bases(0);

		auto var_id=block.first_variant_id_;
		for(; var_id < ref.Variants(ref_seq_id).size() && ref.Variants(ref_seq_id).at(var_id).position_ < end_pos; ++var_id){
			auto &var( ref.Variants(ref_seq_id).at(var_id) );
			chosen_allele = var.FirstAllele();
			tmp_sys_errors.clear();

			// Get last_base
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele)+1 == var.position_ ){
				last_base = sys_last_base_per_allele_.at(chosen_allele);
			}
			else{
				if(var.position_){
					last_base = at(ref.ReferenceSequence(ref_seq_id), var.position_-1);
				}
				else{
					last_base = 4;
				}
			}

			// Get dom_base
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele)+DominantBase::kLastX >= var.position_ ){
				for( auto pos = sys_last_var_pos_per_allele_.at(chosen_allele)+1; pos < var.position_; ++pos ){
					sys_dom_base_per_allele_.at(chosen_allele).Update( at(ref.ReferenceSequence(ref_seq_id), pos) );
				}
			}
			else{
				sys_dom_base_per_allele_.at(chosen_allele).Clear();
				if(var.position_){
					// If we have bases before the variant add them
					sys_dom_base_per_allele_.at(chosen_allele).Set( ref.ReferenceSequence(ref_seq_id), var.position_-1 );
				}
			}

			// Get start_dist_error_region (Variants cannot start an error region)
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < var.position_-block.start_pos_; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = ref.Variants(ref_seq_id).at(var_id-1).position_-block.start_pos_; pos < var.position_-block.start_pos_; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}

			// Get GC (Variants' GC is not taken into account)
			if(var_id && var_id != block.first_variant_id_ && ref.Variants(ref_seq_id).at(var_id-1).position_ + sys_gc_range_/3 > var.position_){ // Updating is ~3 times as much work as recounting
				for(auto pos = ref.Variants(ref_seq_id).at(var_id-1).position_; pos < var.position_; ++pos){
					UpdateGC(gc, gc_bases, pos, ref.ReferenceSequence(ref_seq_id));
				}
			}
			else{
				gc_bases = min(var.position_, static_cast<uintSeqLen>(sys_gc_range_));
				gc = ref.GCContentAbsolut(ref_seq_id, var.position_-gc_bases, var.position_);
			}

			// Draw systematic errors
			if( 0 < length(var.var_seq_)){
				tmp_sys_errors.reserve(length(var.var_seq_));
				for(uintSeqLen vpos=0; vpos < length(var.var_seq_); ++vpos){
					base = at(var.var_seq_, vpos);
					sys_dom_base_per_allele_.at(chosen_allele).Update(base);
					DrawSystematicError(tmp_sys_errors, base, last_base, sys_dom_base_per_allele_.at(chosen_allele).Get(), SafePercent(gc, gc_bases), TransformDistanceToStartOfErrorRegion(start_dist_error_region), start_rate, estimates);
					last_base = base;
				}
			}

			// Insert variant
			block.err_variants_.emplace_back(var.position_-block.start_pos_, tmp_sys_errors, var.allele_);

			// Update information for next variants
			auto ref_allele = ref.NumAlleles(); // For this allele we already used the reference, so in case we have another allele identical to the reference up to the current variance we can copy the values
			if(sys_last_var_pos_per_allele_.at(chosen_allele) >= ref.SequenceLength(ref_seq_id) || sys_last_var_pos_per_allele_.at(chosen_allele)+DominantBase::kLastX < var.position_+length(var.var_seq_)){
				ref_allele = chosen_allele;
			}
			for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
				if(var.InAllele(allele)){
					if( allele != chosen_allele ){ // Dom base for chosen allele was already updated during the drawing process
						// Update dom base
						if(sys_last_var_pos_per_allele_.at(allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(allele)+DominantBase::kLastX >= var.position_+length(var.var_seq_)){
							for( auto pos = sys_last_var_pos_per_allele_.at(allele)+1; pos < var.position_; ++pos ){
								sys_dom_base_per_allele_.at(allele).Update( at(ref.ReferenceSequence(ref_seq_id), pos) );
							}

							for(uintSeqLen vpos=0; vpos < length(var.var_seq_); ++vpos){
								sys_dom_base_per_allele_.at(allele).Update( at(var.var_seq_, vpos) );
							}
						}
						else{
							if(ref_allele < ref.NumAlleles()){
								// We can copy from an identical allele
								sys_dom_base_per_allele_.at(allele) = sys_dom_base_per_allele_.at(ref_allele);
							}
							else{
								ref_allele = allele;

								sys_dom_base_per_allele_.at(allele).Clear();
								if(var.position_){
									sys_dom_base_per_allele_.at(allele).Set( ref.ReferenceSequence(ref_seq_id), var.position_-1 );
								}

								for(uintSeqLen vpos=0; vpos < length(var.var_seq_); ++vpos){
									sys_dom_base_per_allele_.at(allele).Update( at(var.var_seq_, vpos) );
								}
							}
						}
					}

					// Update last base
					sys_last_var_pos_per_allele_.at(allele) = var.position_;
					sys_last_base_per_allele_.at(allele) = last_base;
				}
			}
		}

		// Update start_dist_error_region for rest of block as we don't do this while drawing the systematic errors
		if(sys_from_file_){
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = ref.Variants(ref_seq_id).at(var_id-1).position_-block.start_pos_; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, start_rate, block.sys_errors_.at(pos).second);
				}
			}
		}
	}
}

bool Simulator::CreateBlock(Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	SimBlock *block;
	SimUnit *unit;

	// Create new block
	if(!last_unit_){
		// Create first unit
		if( ref.NumberSequences() ){
			// Take first reference sequence that is long enough
			uintRefSeqId ref_id=0;
			while(ref.NumberSequences() > ref_id && ref.SequenceLength(ref_id) < stats.FragmentDistribution().InsertLengths().to() ){
				++ref_id;
			}
			if(ref.NumberSequences() > ref_id){
				if(!CreateUnit(ref_id, 0, ref, stats, estimates, block, unit)){
					return false;
				}
				block = new SimBlock(1, 0, block, block_seed_gen_());

				first_unit_ = unit;
				last_unit_ = unit;
				current_unit_ = unit;
				unit->first_block_ = block;
				unit->last_block_ = block;
				current_block_ = block;
			}
			else{
				// There is no reference sequence that is long enough
				printErr << "All reference sequences are too short for simulating. They should have at least " << stats.FragmentDistribution().InsertLengths().to() << " bases" << std::endl;
				return false;
			}
		}
	}
	else{
		if( last_unit_->last_block_->start_pos_ + kBlockSize >= ref.SequenceLength(last_unit_->ref_seq_id_) ){
			// Create next unit as the current one is completed: Go to the next reference sequence that is long enough if exists
			uintRefSeqId ref_id=last_unit_->ref_seq_id_+1;
			while(ref.NumberSequences() > ref_id && ref.SequenceLength(ref_id) < stats.FragmentDistribution().InsertLengths().to() ){
				++ref_id;
			}
			if( ref.NumberSequences() > ref_id ){
				if(!CreateUnit(ref_id, last_unit_->last_block_->id_+1, ref, stats, estimates, block, unit)){
					return false;
				}
				block = new SimBlock(last_unit_->last_block_->id_+1, 0, block, block_seed_gen_());

				last_unit_->next_unit_ = unit;
				last_unit_ = unit;
				unit->first_block_ = block;
				unit->last_block_ = block;
			}
			else{
				// No new reference sequence anymore: Simulation is complete
				return false;
			}
		}
		else{
			// The new block can be created in the current unit
			unit = last_unit_;

			block = new SimBlock(unit->last_block_->id_+1, unit->last_block_->start_pos_ + kBlockSize, unit->last_block_->partner_block_->partner_block_, block_seed_gen_());
			if(ref.VariantsLoaded()){
				// Find first variant that is in the new block (or one of the upcoming ones in case this block does not have variation), which is one after the last variant(reverse first) of the previous block
				block->first_variant_id_ = unit->last_block_->partner_block_->first_variant_id_+1;
			}
			if(ref.MethylationLoaded()){
				// Find first methylation region that is in the new block (or one of the upcoming ones in case this block does not have methylation)
				block->first_methylation_id_ = unit->last_block_->partner_block_->first_methylation_id_;
				if(0 > block->first_methylation_id_ || ref.UnmethylatedRegions(unit->ref_seq_id_).at(block->first_methylation_id_).second <= block->start_pos_ ){
					++block->first_methylation_id_;
				}

			}
			unit->last_block_->next_block_ = block;
			unit->last_block_ = block;
		}
	}

	if(!simulation_error_){
		// Set systematic errors for new block
		auto end_pos = min(block->start_pos_ + kBlockSize, ref.SequenceLength(unit->ref_seq_id_));
		block->sys_errors_.reserve(end_pos - block->start_pos_);
		if(sys_from_file_){
			ReadSystematicErrors(block->sys_errors_, block->start_pos_, end_pos);
			SetSystematicErrorVariantsForward(distance_to_start_of_error_region_, start_error_rate_, *block, unit->ref_seq_id_, end_pos, ref, stats, estimates);
		}
		else{
			auto tmp_distance_to_start_of_error_region = distance_to_start_of_error_region_;
			auto tmp_start_error_rate = start_error_rate_;
			SetSystematicErrors(block->sys_errors_, ref.ReferenceSequence(unit->ref_seq_id_), block->start_pos_, end_pos, stats, estimates);
			SetSystematicErrorVariantsForward(tmp_distance_to_start_of_error_region, tmp_start_error_rate, *block, unit->ref_seq_id_, end_pos, ref, stats, estimates);
		}

		// Delete old blocks and units that are not needed anymore
		SimBlock *del_block(first_unit_->first_block_);
		SimUnit *del_unit(first_unit_);
		while( del_block->finished_ ){
			if(first_unit_->first_block_->next_block_){
				first_unit_->first_block_ = first_unit_->first_block_->next_block_;
			}
			else if(first_unit_->next_unit_ && first_unit_->next_unit_->first_block_){
				first_unit_ = first_unit_->next_unit_;
				ref.ClearVariants(del_unit->ref_seq_id_ + 1);
				ref.ClearMethylation(del_unit->ref_seq_id_ + 1);
				delete del_unit;
				del_unit = first_unit_;
			}
			else{
				printErr << "Ran out of simulation blocks, but simulation is not complete.";
				simulation_error_ = true;
				return false;
			}

			delete del_block->partner_block_;
			delete del_block;
			del_block = first_unit_->first_block_;
		}

		return true;
	}
	else{
		return false;
	}
}

bool Simulator::GetNextBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&block, SimUnit *&unit ){
	// See if we can already read in more variants, so we are always one reference sequence ahead of the simulation
	if(!ref.VariantsCompletelyLoaded() && current_unit_ && !ref.VariantsLoadedForSequence(current_unit_->ref_seq_id_+2)){
		if( var_read_mutex_.try_lock() ){
			if( !ref.ReadVariants(current_unit_->ref_seq_id_+2) ){
				var_read_mutex_.unlock();
				return false;
			}

			var_read_mutex_.unlock();
		}
	}

	if(!ref.MethylationCompletelyLoaded() && current_unit_ && !ref.MethylationLoadedForSequence(current_unit_->ref_seq_id_+2)){
		if( methylation_read_mutex_.try_lock() ){
			if( !ref.ReadMethylation(current_unit_->ref_seq_id_+2) ){
				methylation_read_mutex_.unlock();
				return false;
			}

			methylation_read_mutex_.unlock();
		}
	}

	lock_guard<mutex> lock(block_creation_mutex_);

	// Create a new block to keep the buffer for long fragments
	if( !CreateBlock(ref, stats, estimates) ){
		return false;
	}

	if(!current_unit_){
		// No new reference sequence anymore: Simulation is complete
		return false;
	}

	// Set block and unit to simulate from
	unit = current_unit_;
	block = current_block_;

	// Find the next block for next time
	if( current_block_->next_block_ ){
		// Still blocks available in the current unit
		current_block_ = current_block_->next_block_;
	}
	else{
		// No new blocks in current unit: Get next unit if exists
		current_unit_ = current_unit_->next_unit_;
		if(current_unit_){
			current_block_ = current_unit_->first_block_;
		}
	}

	return true;
}

inline bool Simulator::VariantInsideCurrentFragment( intVariantId cur_var_id, uintSeqLen cur_end_position, intSeqShift end_pos_shift, uintRefSeqId ref_seq_id, const Reference &ref ) const{
	 // variant with that id exists && ( is relevant to end surrounding || is relevant to start surrounding) // if insert length is very short start surrounding can be longer than end surrounding
	return ref.Variants(ref_seq_id).size() > cur_var_id && ( ref.Variants(ref_seq_id).at(cur_var_id).position_ < cur_end_position + end_pos_shift );
}

void Simulator::HandleGCModAndEndPosShiftForNewVariants(
		VariantBiasVarModifiers &bias_mod,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_end_position) const{
	while( VariantInsideCurrentFragment(bias_mod.unhandled_variant_id_.at(allele), cur_end_position, bias_mod.end_pos_shift_.at(allele), ref_seq_id, ref ) && 0==bias_mod.unhandled_bases_in_variant_.at(allele) ){
		auto &var( ref.Variants(ref_seq_id).at(bias_mod.unhandled_variant_id_.at(allele)) );

		if(!var.InAllele(allele)){
			++bias_mod.unhandled_variant_id_.at(allele);
		}
		else{
			// ---GC modification---
			// Add GC from variant
			for(uintSeqLen pos=0; pos < length(var.var_seq_) && var.position_ + pos < cur_end_position + bias_mod.end_pos_shift_.at(allele); ++pos){
				if( IsGC(var.var_seq_, pos) ){
					++bias_mod.gc_mod_.at(allele);
				}
			}

			// Remove GC from reference
			if( ref.GC(ref_seq_id, var.position_) ){
				--bias_mod.gc_mod_.at(allele);
			}

			// ---End pos shift---
			if(0 == length(var.var_seq_)){
				// Deletion
				++bias_mod.end_pos_shift_.at(allele);
				++bias_mod.unhandled_variant_id_.at(allele);
			}
			else if(1 == length(var.var_seq_)){
				++bias_mod.unhandled_variant_id_.at(allele);
			}
			else{
				// Insertion
				if(var.position_ + length(var.var_seq_) <= cur_end_position + bias_mod.end_pos_shift_.at(allele)){
					// Completely in current fragment
					bias_mod.end_pos_shift_.at(allele) -= length(var.var_seq_) - 1;
					++bias_mod.unhandled_variant_id_.at(allele);
				}
				else{
					// Reaches out of fragment
					bias_mod.unhandled_bases_in_variant_.at(allele) = var.position_ + length(var.var_seq_) - (cur_end_position + bias_mod.end_pos_shift_.at(allele));
					bias_mod.end_pos_shift_.at(allele) -= length(var.var_seq_) - bias_mod.unhandled_bases_in_variant_.at(allele) - 1;
					// Variant is not handled completely so don't increase unhandled_variant_id
				}
			}
		}
	}
}



void Simulator::HandleSurroundingVariantsBeforeCenter(
		Surrounding &mod_surrounding,
		uintSeqLen center_position,
		intSeqShift initial_pos_shift,
		intVariantId center_var,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		bool reverse) const{
	auto cur_var = center_var;
	auto pos_shift = initial_pos_shift;
	intSeqShift sur_pos;
	while( 0 <= --cur_var ){ // + Break statements at beginning
		if(reverse){
			sur_pos = static_cast<intSeqShift>(center_position) - ref.Variants(ref_seq_id).at(cur_var).position_ + pos_shift;
		}
		else{
			sur_pos = static_cast<intSeqShift>(ref.Variants(ref_seq_id).at(cur_var).position_) - center_position + pos_shift;
		}

		if( 0 > sur_pos || Surrounding::Length() <= sur_pos ){
			break;
		}

		if(ref.Variants(ref_seq_id).at(cur_var).InAllele(allele)){
			if(0 == length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
				// Deletion
				if(reverse){
					auto new_base_pos = static_cast<intSeqShift>(center_position) + pos_shift - static_cast<intSeqShift>(Surrounding::Length());
					if(0 > new_base_pos){
						mod_surrounding.DeleteSurroundingBaseShiftingOnRightSide(sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), ref.SequenceLength(ref_seq_id) + new_base_pos)));
					}
					else{
						mod_surrounding.DeleteSurroundingBaseShiftingOnRightSide(sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), new_base_pos)));
					}
					--pos_shift;
				}
				else{
					if(++pos_shift > center_position){
						mod_surrounding.DeleteSurroundingBaseShiftingOnLeftSide(sur_pos, at(ref.ReferenceSequence(ref_seq_id), ref.SequenceLength(ref_seq_id) + center_position - pos_shift));
					}
					else{
						mod_surrounding.DeleteSurroundingBaseShiftingOnLeftSide(sur_pos, at(ref.ReferenceSequence(ref_seq_id), center_position - pos_shift));
					}
				}
			}
			else{
				// Base modification (Insertions may also include a base modification and we skip checking first if we exchange it with the same)
				if(reverse){
					mod_surrounding.ChangeSurroundingBase(sur_pos, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)));
				}
				else{
					mod_surrounding.ChangeSurroundingBase(sur_pos, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0));
				}

				if(1 < length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
					// Insertion (First base is the already handled)
					if(reverse){
						mod_surrounding.InsertSurroundingBasesShiftingOnRightSide(sur_pos, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1)));
						pos_shift += length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
					}
					else{
						mod_surrounding.InsertSurroundingBasesShiftingOnLeftSide(sur_pos, suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1));
						pos_shift -= length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
					}
				}
			}
		}
	}
}

void Simulator::HandleSurroundingVariantsAfterCenter(
		Surrounding &mod_surrounding,
		uintSeqLen center_position,
		intSeqShift initial_pos_shift,
		intVariantId center_var,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		bool reverse) const{
	auto pos_shift = initial_pos_shift;
	intSeqShift sur_pos;
	for(auto cur_var = center_var; cur_var < ref.Variants(ref_seq_id).size(); ++cur_var){ // + Break statements at beginning
		if(reverse){
			sur_pos = static_cast<intSeqShift>(center_position) - ref.Variants(ref_seq_id).at(cur_var).position_ + pos_shift;
		}
		else{
			sur_pos = static_cast<intSeqShift>(ref.Variants(ref_seq_id).at(cur_var).position_) - center_position + pos_shift;
		}

		if( 0 > sur_pos || Surrounding::Length() <= sur_pos ){
			break;
		}

		if(ref.Variants(ref_seq_id).at(cur_var).InAllele(allele)){
			if(0 == length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
				// Deletion
				if(reverse){
					++pos_shift;
					mod_surrounding.DeleteSurroundingBaseShiftingOnLeftSide(sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), (center_position + pos_shift)%ref.SequenceLength(ref_seq_id))));
				}
				else{
					mod_surrounding.DeleteSurroundingBaseShiftingOnRightSide(sur_pos, at(ref.ReferenceSequence(ref_seq_id), (center_position - pos_shift + Surrounding::Length())%ref.SequenceLength(ref_seq_id)));
					--pos_shift;
				}
			}
			else{
				// Base modification (Insertions may also include a base modification and we skip checking first if we exchange it with the same)
				if(reverse){
					mod_surrounding.ChangeSurroundingBase(sur_pos, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)));
				}
				else{
					mod_surrounding.ChangeSurroundingBase(sur_pos, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0));
				}

				if(1 < length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
					// Insertion (First base is the already handled)
					if(reverse){
						if(sur_pos){
							mod_surrounding.InsertSurroundingBasesShiftingOnLeftSide(sur_pos-1, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1)));
							pos_shift -= length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
						}
					}
					else{
						if(sur_pos+1 < Surrounding::Length()){
							mod_surrounding.InsertSurroundingBasesShiftingOnRightSide(sur_pos+1, suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1));
							pos_shift += length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
						}
					}
				}
			}
		}
	}
}

void Simulator::VariantModStartSurrounding(
		VariantBiasVarModifiers &bias_mod,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		const Surrounding &surrounding_start) const{
	// Variants in the roll around are ignored for variant modifications of surrounding (but to fill in deletions still use the proper roll around base)
	// Copy reference surrounding into allele
	bias_mod.surrounding_start_.at(allele) = surrounding_start;

	if(ref.Variants(ref_seq_id).size()){
		// Handle variants before center (Shifting to/from the left)
		intSeqShift pos_shift = Surrounding::kStartPos;
		auto cur_var = bias_mod.first_variant_id_;
		if(bias_mod.start_variant_pos_){
			// Handle the base substitution (as the insertion might include it)
			bias_mod.surrounding_start_.at(allele).ChangeSurroundingBase(pos_shift, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0) );
			// Partial insertion (as the rest is shifted to the right)
			bias_mod.surrounding_start_.at(allele).InsertSurroundingBasesShiftingOnLeftSide(pos_shift, infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, bias_mod.start_variant_pos_+1));
			pos_shift -= bias_mod.start_variant_pos_;
		}

		HandleSurroundingVariantsBeforeCenter(bias_mod.surrounding_start_.at(allele), cur_start_position, pos_shift, cur_var, allele, ref_seq_id, ref, false);

		// Handle variants after center (Shifting to/from the right)
		pos_shift = Surrounding::kStartPos;
		cur_var = bias_mod.first_variant_id_;
		if(bias_mod.start_variant_pos_){
			if(length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) > bias_mod.start_variant_pos_+1){
				// Partial insertion (as the rest is shifted to the left)
				if(pos_shift+1 < Surrounding::Length()){
					bias_mod.surrounding_start_.at(allele).InsertSurroundingBasesShiftingOnRightSide(pos_shift+1, suffix(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, bias_mod.start_variant_pos_+1));
					pos_shift += length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) - bias_mod.start_variant_pos_ - 1;
				}
			}

			// Skip this insertion for the subsequent variant handling
			++cur_var;
		}

		HandleSurroundingVariantsAfterCenter(bias_mod.surrounding_start_.at(allele), cur_start_position, pos_shift, cur_var, allele, ref_seq_id, ref, false);
	}
}

void Simulator::PrepareBiasModForCurrentStartPos(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		uintSeqLen first_fragment_length,
		const Surrounding &surrounding_start) const{
	auto cur_end_position = cur_start_position + first_fragment_length - 1;

	if(ref.VariantsLoaded()){
		// Store the modifications from the reference sequence to select correct biases for counts generation
		bias_mod.unhandled_variant_id_.clear();
		bias_mod.unhandled_variant_id_.resize(ref.NumAlleles(), bias_mod.first_variant_id_);
		bias_mod.unhandled_bases_in_variant_.clear();
		bias_mod.unhandled_bases_in_variant_.resize(ref.NumAlleles(), 0);
		bias_mod.gc_mod_.clear();
		bias_mod.gc_mod_.resize(ref.NumAlleles(), 0);
		bias_mod.end_pos_shift_.clear();
		bias_mod.end_pos_shift_.resize(ref.NumAlleles(), 0);

		if(bias_mod.start_variant_pos_){
			// ---Handle gc changes and end pos shift due to insertion at the beginning---
			// Add GC from variant
			for(uintSeqLen pos=bias_mod.start_variant_pos_; pos < length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ + pos - bias_mod.start_variant_pos_ < cur_end_position; ++pos){
				if( IsGC(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, pos) ){
					++bias_mod.gc_mod_.at(0);
				}
			}

			// Remove GC from reference as we are starting from the non-reference part
			if( ref.GC(ref_seq_id, cur_start_position) ){
				--bias_mod.gc_mod_.at(0);
			}

			// End pos shift cannot be longer than insert length or variant length
			bias_mod.end_pos_shift_.at(0) = static_cast<intSeqShift>(1) - static_cast<intSeqShift>( min(cur_end_position-cur_start_position, static_cast<uintSeqLen>(length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_))-bias_mod.start_variant_pos_) );

			++bias_mod.unhandled_variant_id_.at(0); // Handled first_variant_id_

			for(uintAlleleId allele=1; allele < ref.NumAlleles(); ++allele){
				bias_mod.gc_mod_.at(allele) = bias_mod.gc_mod_.at(0);
				bias_mod.end_pos_shift_.at(allele) = bias_mod.end_pos_shift_.at(0);
				++bias_mod.unhandled_variant_id_.at(allele);
			}
		}

		// Check for every allele if variants cause changes to gc or surrounding
		for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
			if( !AlleleSkipped(bias_mod, allele, ref.Variants(ref_seq_id), cur_start_position) ){
				VariantModStartSurrounding(bias_mod, allele, ref_seq_id, ref, cur_start_position, surrounding_start); // Must be before HandleGCModAndEndPosShiftForNewVariants, so that bias_mod.unhandled_variant_id_ is untouched
				HandleGCModAndEndPosShiftForNewVariants(bias_mod, allele, ref_seq_id, ref, cur_end_position);
			}
		}
	}

	bias_mod.last_end_position_.clear();
	bias_mod.last_end_position_.resize(ref.NumAlleles(), cur_end_position);
	bias_mod.last_gc_ = 0;
	bias_mod.last_gc_end_ = cur_start_position;
}

void Simulator::VariantModEndSurrounding(
		VariantBiasVarModifiers &bias_mod,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen last_position) const{
	// Variants in the roll around are ignored for variant modifications of surrounding (but to fill in deletions still use the proper roll around base)
	// Copy reference surrounding into allele if it doesn't need to much shifting afterwards
	if(ref.VariantsLoaded() && ref.Variants(ref_seq_id).size()){
		// Handle variants after center (Shifting to/from the left)
		intSeqShift pos_shift = Surrounding::kStartPos;
		auto cur_var = bias_mod.unhandled_variant_id_.at(allele); // At the variants where pos == last_position if exists or after it otherwise
		if( bias_mod.unhandled_bases_in_variant_.at(allele) ){
			// Partial insertion (as the rest is shifted to the right)
			bias_mod.surrounding_end_.at(allele).InsertSurroundingBasesShiftingOnLeftSide(pos_shift, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele))));
			pos_shift -= bias_mod.unhandled_bases_in_variant_.at(allele);
			++cur_var; // This insertion is handled, do not handle at again in HandleSurroundingVariantsAfterCenter
		}
		else if( bias_mod.start_variant_pos_ && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == last_position && length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) > bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1){
			if(pos_shift){
				// Partial insertion (as the rest is shifted to the right)
				bias_mod.surrounding_end_.at(allele).InsertSurroundingBasesShiftingOnLeftSide(pos_shift-1, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1)));
				pos_shift -= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) - (bias_mod.start_variant_pos_ - bias_mod.end_pos_shift_.at(allele) + 1);
			}
		}

		HandleSurroundingVariantsAfterCenter(bias_mod.surrounding_end_.at(allele), last_position, pos_shift, cur_var, allele, ref_seq_id, ref, true);

		// Handle variants before center (Shifting to/from the right)
		pos_shift = Surrounding::kStartPos;
		cur_var = bias_mod.unhandled_variant_id_.at(allele); // At the variants where pos == last_position if exists or after it otherwise
		if( bias_mod.unhandled_bases_in_variant_.at(allele) ){
			if(pos_shift+1 < Surrounding::Length()){
				// Handle the base substitution (as the insertion might include it)
				bias_mod.surrounding_end_.at(allele).ChangeSurroundingBase(pos_shift+1, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)) );

				// Partial insertion (as the rest is shifted to the left)
				bias_mod.surrounding_end_.at(allele).InsertSurroundingBasesShiftingOnRightSide(pos_shift+1, ReverseComplementorDna(infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele))));
				pos_shift += length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele)-1;
			}
		}
		else if( bias_mod.start_variant_pos_ && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == last_position ){
			cur_var = bias_mod.first_variant_id_;
			// Handle the base substitution (as the insertion might include it)
			bias_mod.surrounding_end_.at(allele).ChangeSurroundingBase(pos_shift, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)) );

			// Partial insertion (as the rest is shifted to the left)
			bias_mod.surrounding_end_.at(allele).InsertSurroundingBasesShiftingOnRightSide(pos_shift, ReverseComplementorDna(infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1)));
			pos_shift += bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele);
		}

		HandleSurroundingVariantsBeforeCenter(bias_mod.surrounding_end_.at(allele), last_position, pos_shift, cur_var, allele, ref_seq_id, ref, true);
	}
}

void Simulator::UpdateBiasModForCurrentFragmentLength(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		uintSeqLen cur_end_position,
		uintSeqLen last_end_position,
		uintAlleleId allele ) const{
	if( ref.VariantsLoaded() && cur_end_position > bias_mod.last_end_position_.at(allele) ){
		bool need_new_variants = false;
		if(bias_mod.start_variant_pos_ && last_end_position+1-cur_start_position <= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_)-bias_mod.start_variant_pos_){
			// Still in start variant insertion
			auto &var( ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_) );

			auto stop_pos = bias_mod.start_variant_pos_ + cur_end_position-cur_start_position;
			if( stop_pos > length(var.var_seq_) ){
				stop_pos = length(var.var_seq_);
				need_new_variants = true;
			}

			auto start_pos = bias_mod.start_variant_pos_ + last_end_position+1-cur_start_position - 1;
			bias_mod.end_pos_shift_.at(allele) -= stop_pos - start_pos;
			for(auto pos = start_pos; pos < stop_pos; ++pos){
				if( IsGC(var.var_seq_, pos) ){
					++bias_mod.gc_mod_.at(allele);
				}
			}
		}
		else if(bias_mod.unhandled_bases_in_variant_.at(allele)){
			// Unhandled bases only occur for insertions and the reference base has already been covered
			auto &var( ref.Variants(ref_seq_id).at(bias_mod.unhandled_variant_id_.at(allele)) );

			// ---GC modification---
			// Add GC from variant
			auto start_pos = length(var.var_seq_) - bias_mod.unhandled_bases_in_variant_.at(allele);
			auto stop_pos = start_pos + cur_end_position-last_end_position;

			if(stop_pos > length(var.var_seq_)){
				stop_pos = length(var.var_seq_);
				need_new_variants = true;
			}

			bias_mod.end_pos_shift_.at(allele) -= stop_pos - start_pos;
			bias_mod.unhandled_bases_in_variant_.at(allele) -= stop_pos - start_pos;

			for(auto pos = start_pos; pos < stop_pos; ++pos){
				if( IsGC(var.var_seq_, pos) ){
					++bias_mod.gc_mod_.at(allele);
				}
			}

			if( 0 == bias_mod.unhandled_bases_in_variant_.at(allele) ){
				++bias_mod.unhandled_variant_id_.at(allele);
			}
		}
		else{
			need_new_variants = true;
		}

		if(need_new_variants){
			HandleGCModAndEndPosShiftForNewVariants(bias_mod, allele, ref_seq_id, ref, cur_end_position);
		}
	}
}

inline void Simulator::PrepareEndSurroundingsForCurrentFragmentLength(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_end_position,
		uintAlleleId allele ) const{
	auto corrected_pos = cur_end_position + bias_mod.end_pos_shift_.at(allele);
	if(corrected_pos < ref.SequenceLength(ref_seq_id)){
		// Last surrounding is too far away: create new one
		ref.ReverseSurrounding( bias_mod.surrounding_end_.at(allele), ref_seq_id, corrected_pos );
		VariantModEndSurrounding( bias_mod, allele, ref_seq_id, ref, corrected_pos );
	}
}

void Simulator::PrepareBiasModForCurrentFragmentLength(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		uintSeqLen fragment_length,
		uintAlleleId allele ) const{
	auto cur_end_position = cur_start_position + fragment_length - 1;

	UpdateBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, ref, cur_start_position, cur_end_position, bias_mod.last_end_position_.at(allele), allele );

	if( bias_mod.start_variant_pos_ && fragment_length <= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_)-bias_mod.start_variant_pos_ ){
		UpdateBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, ref, cur_start_position, cur_end_position+1, cur_end_position, allele );
		PrepareEndSurroundingsForCurrentFragmentLength( bias_mod, ref_seq_id, ref, cur_end_position, allele);
	}
	else{
		PrepareEndSurroundingsForCurrentFragmentLength( bias_mod, ref_seq_id, ref, cur_end_position, allele);
		UpdateBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, ref, cur_start_position, cur_end_position+1, cur_end_position, allele );
	}

	bias_mod.last_end_position_.at(allele) = cur_end_position+1;
}

reseq::uintPercent Simulator::GetGCPercent(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_end_position,
		uintSeqLen fragment_length,
		uintAlleleId allele ){
	if( cur_end_position < bias_mod.last_gc_end_ ){
		return Percent(bias_mod.last_gc_ - ref.GCContentAbsolut( ref_seq_id, cur_end_position, bias_mod.last_gc_end_ ) + bias_mod.gc_mod_.at(allele), fragment_length);
	}
	else{
		bias_mod.last_gc_ += ref.GCContentAbsolut( ref_seq_id, bias_mod.last_gc_end_, cur_end_position );
		bias_mod.last_gc_end_ = cur_end_position;
		return Percent(bias_mod.last_gc_ + bias_mod.gc_mod_.at(allele), fragment_length);
	}
}

void Simulator::CheckForInsertedBasesToStartFrom( VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position, const Reference &ref ) const{
	if(ref.VariantsLoaded() && bias_mod.first_variant_id_ < ref.Variants(ref_seq_id).size() && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == cur_start_position){
		if(bias_mod.start_variant_pos_){
			if(++bias_mod.start_variant_pos_ >= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_)){
				bias_mod.start_variant_pos_ = 0;
				++bias_mod.first_variant_id_; // We have handled the current variant in this while loop pass
			}
		}
		else{
			while(bias_mod.first_variant_id_ < ref.Variants(ref_seq_id).size() && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == cur_start_position && 2 > length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) ){
				++bias_mod.first_variant_id_; // All deletions and substitutions are handled with the first variant of the position
			}
		}

		if(0 == bias_mod.start_variant_pos_){
			if( bias_mod.first_variant_id_ < ref.Variants(ref_seq_id).size() && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == cur_start_position ){
				// We have an insertion to handle
				bias_mod.start_variant_pos_ = 1; // Starting from the first base is handled during the first pass of the while loop, so now we use the second, third, ... as start positions
			}
		}
	}
}

void Simulator::GetOrgSeq(
		SimPair &sim_reads,
		bool strand,
		uintAlleleId allele,
		uintSeqLen fragment_length,
		uintSeqLen cur_start_position,
		uintSeqLen cur_end_position,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		const DataStats &stats,
		const VariantBiasVarModifiers &bias_mod ) const{
	if(ref.VariantsLoaded()){
		// Forward
		ref.ReferenceSequence( sim_reads.org_seq_.at(strand), ref_seq_id, cur_start_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(strand).to()+stats.Errors().MaxLenDeletion())), false, ref.Variants(ref_seq_id), bias_mod.StartVariant(), allele );

		// Reverse
		ref.ReferenceSequence( sim_reads.org_seq_.at(!strand), ref_seq_id, cur_end_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(!strand).to()+stats.Errors().MaxLenDeletion())), true, ref.Variants(ref_seq_id), bias_mod.EndVariant(ref.Variants(ref_seq_id), cur_end_position, allele), allele );
	}
	else{
		// Forward
		ref.ReferenceSequence( sim_reads.org_seq_.at(strand), ref_seq_id, cur_start_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(strand).to()+stats.Errors().MaxLenDeletion())), false );

		// Reverse
		ref.ReferenceSequence( sim_reads.org_seq_.at(!strand), ref_seq_id, cur_end_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(!strand).to()+stats.Errors().MaxLenDeletion())), true );
	}
}

void Simulator::CTConversion(
		DnaString &read,
		const Reference &ref,
		uintRefSeqId seq_id,
		uintSeqLen start_pos,
		uintAlleleId allele,
		intVariantId cur_methylation_start,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen,
		bool reversed ) const{
	auto cur_meth = cur_methylation_start;
	uintReadLen read_pos = 0;
	auto ref_pos = start_pos;

	auto &conversion_rate = ref.Unmethylation(seq_id, allele);
	auto &regions = ref.UnmethylatedRegions(seq_id);

	if(reversed){
		// Bring cur_meth to last in region
		while( cur_meth < regions.size() && regions.at(cur_meth).first <= ref_pos ){
			++cur_meth;
		}
		--cur_meth; // Go from one behind the current region to the last in it

		// Jump over region without information
		if( cur_meth && regions.at(cur_meth).second <= ref_pos ){
			read_pos += ref_pos - (regions.at(cur_meth).second-1);
			ref_pos = regions.at(cur_meth).second-1;
		}

		while( cur_meth && read_pos < length(read) ){
			// Convert C->T in region
			while( ref_pos >= regions.at(cur_meth).first && read_pos < length(read) ){
				if( 1 == at(read, read_pos) ){
					if( rdist.ZeroToOne(rgen) < conversion_rate.at(cur_meth) ){
						at(read, read_pos) = 3;
					}
				}

				--ref_pos;
				++read_pos;
			}

			// Jump over region without information
			if( --cur_meth && regions.at(cur_meth).second <= ref_pos ){
				read_pos += ref_pos - (regions.at(cur_meth).second-1);
				ref_pos = regions.at(cur_meth).second-1;
			}
		}
	}
	else{
		// Jump over region without information
		if( cur_meth < regions.size() && regions.at(cur_meth).first > ref_pos ){
			read_pos += regions.at(cur_meth).first - ref_pos;
			ref_pos = regions.at(cur_meth).first;
		}

		while( cur_meth < regions.size() && read_pos < length(read) ){
			// Convert C->T in region
			while( ref_pos < regions.at(cur_meth).second && read_pos < length(read) ){
				if( 1 == at(read, read_pos) ){
					if( rdist.ZeroToOne(rgen) < conversion_rate.at(cur_meth) ){
						at(read, read_pos) = 3;
					}
				}

				++ref_pos;
				++read_pos;
			}

			// Jump over region without information
			if( ++cur_meth < regions.size() ){
				read_pos += regions.at(cur_meth).first - ref_pos;
				ref_pos = regions.at(cur_meth).first;
			}
		}
	}
}

void Simulator::CTConversion(
		DnaString &read,
		const Reference &ref,
		uintRefSeqId seq_id,
		uintSeqLen start_pos,
		uintAlleleId allele,
		intVariantId cur_methylation_start,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen,
		bool reversed,
		const vector<Reference::Variant> &variants,
		pair<intVariantId, uintSeqLen> first_variant // {id, posCurrentlyAt}
		) const{
	auto cur_meth = cur_methylation_start;
	uintReadLen read_pos = 0;
	auto ref_pos = start_pos;

	auto &conversion_rate = ref.Unmethylation(seq_id, allele);
	auto &regions = ref.UnmethylatedRegions(seq_id);

	auto cur_var = first_variant.first;
	uintSeqLen var_bases_left = 0;

	bool deletion;

	if(reversed){
		// Set var_bases_left if we start in a variant
		if(0 <= cur_var && variants.at(cur_var).position_ == ref_pos && 1 < length(variants.at(cur_var).var_seq_) && variants.at(cur_var).InAllele(allele)){
			// We cannot start at a deletions and we don't care about substitutions, so only handle insertions
			var_bases_left = length(variants.at(cur_var).var_seq_) - first_variant.second;
		}

		// Bring cur_meth to last in region
		while( cur_meth < regions.size() && regions.at(cur_meth).first <= ref_pos ){
			++cur_meth;
		}
		--cur_meth; // Go from one behind the current region to the last in it

		// Jump over region without information
		if( 0 <= cur_meth && regions.at(cur_meth).second <= ref_pos ){
			// Handle start variant
			if(var_bases_left){
				read_pos += var_bases_left;
				var_bases_left = 0;
				--ref_pos;
				--cur_var;
			}

			// Handle all variants before the next unmethylated region
			while( 0 <= cur_var && regions.at(cur_meth).second <= variants.at(cur_var).position_ && read_pos < length(read) ){
				if( variants.at(cur_var).InAllele(allele) ){
					read_pos += ref_pos - variants.at(cur_var).position_;
					read_pos += length(variants.at(cur_var).var_seq_);
					ref_pos = variants.at(cur_var).position_-1; // Go one before the variant, because we handle the variant here, and it is only at exactly one ref_pos, so cannot reach into unmethylated region
				}
				--cur_var;
			}

			// Handle from the last variant before the region to the beginning of the region
			read_pos += ref_pos - (regions.at(cur_meth).second-1);
			ref_pos = regions.at(cur_meth).second-1;
		}

		while( 0 <= cur_meth && read_pos < length(read) ){
			// Convert C->T in region
			while( ref_pos >= regions.at(cur_meth).first && read_pos < length(read) ){
				if( 0 == var_bases_left ){
					while( 0 <= cur_var && variants.at(cur_var).position_ == ref_pos && !variants.at(cur_var).InAllele(allele) ){
						--cur_var;
					}

					if( 0 <= cur_var && variants.at(cur_var).position_ == ref_pos ){
						if( 0 == length(variants.at(cur_var).var_seq_) ){
							deletion = true;
							--cur_var;
						}
						else{
							var_bases_left = length(variants.at(cur_var).var_seq_);
						}
					}
				}

				if( deletion ){
					deletion = false;
				}
				else{
					if( 1 == at(read, read_pos) ){
						if( rdist.ZeroToOne(rgen) < conversion_rate.at(cur_meth) ){
							at(read, read_pos) = 3;
						}
					}
				}

				if(var_bases_left){
					if(0 == --var_bases_left){
						--cur_var;
					}
				}
				if(0 == var_bases_left){
					--ref_pos;
				}
				++read_pos;
			}

			// Jump over region without information
			if( 0 <= --cur_meth && regions.at(cur_meth).second <= ref_pos ){
				// Handle all variants before the next unmethylated region
				while( 0 <= cur_var && regions.at(cur_meth).second <= variants.at(cur_var).position_ && read_pos < length(read) ){
					if( variants.at(cur_var).InAllele(allele) ){
						read_pos += ref_pos - variants.at(cur_var).position_;
						read_pos += length(variants.at(cur_var).var_seq_);
						ref_pos = variants.at(cur_var).position_-1; // Go one before the variant, because we handle the variant here, and it is only at exactly one ref_pos, so cannot reach into unmethylated region
					}
					--cur_var;
				}

				// Handle from the last variant before the region to the beginning of the region
				read_pos += ref_pos - (regions.at(cur_meth).second-1);
				ref_pos = regions.at(cur_meth).second-1;
			}
		}
	}
	else{
		// Set var_bases_left if we start in a variant
		if( cur_var < variants.size() && variants.at(cur_var).position_ == ref_pos && 1 < length(variants.at(cur_var).var_seq_) && variants.at(cur_var).InAllele(allele) ){
			// We cannot start at a deletions and we don't care about substitutions, so only handle insertions
			var_bases_left = length(variants.at(cur_var).var_seq_) - first_variant.second;
		}

		// Jump over region without information
		if( cur_meth < regions.size() && regions.at(cur_meth).first > ref_pos ){
			// Handle start variant
			if(var_bases_left){
				read_pos += var_bases_left;
				var_bases_left = 0;
				++ref_pos;
				++cur_var;
			}

			// Handle all variants before the next unmethylated region
			while( cur_var < variants.size() && regions.at(cur_meth).first > variants.at(cur_var).position_ && read_pos < length(read) ){
				if( variants.at(cur_var).InAllele(allele) ){
					read_pos += variants.at(cur_var).position_ - ref_pos;
					read_pos += length(variants.at(cur_var).var_seq_);
					ref_pos = variants.at(cur_var).position_+1; // Go one after the variant, because we handle the variant here, and it is only at exactly one ref_pos, so cannot reach into unmethylated region
				}
				++cur_var;
			}

			// Handle from the last variant before the region to the beginning of the region
			read_pos += regions.at(cur_meth).first - ref_pos;
			ref_pos = regions.at(cur_meth).first;
		}

		while( cur_meth < regions.size() && read_pos < length(read) ){
			// Convert C->T in region
			while( ref_pos < regions.at(cur_meth).second && read_pos < length(read) ){
				if(0 == var_bases_left){
					while( cur_var < variants.size() && variants.at(cur_var).position_ == ref_pos && !variants.at(cur_var).InAllele(allele) ){
						++cur_var;
					}

					if( cur_var < variants.size() && variants.at(cur_var).position_ == ref_pos ){
						if( 0 == length(variants.at(cur_var).var_seq_) ){
							deletion = true;
							++cur_var;
						}
						else{
							var_bases_left = length(variants.at(cur_var).var_seq_);
						}
					}
				}

				if( deletion ){
					deletion = false;
				}
				else{
					if( 1 == at(read, read_pos) ){
						if( rdist.ZeroToOne(rgen) < conversion_rate.at(cur_meth) ){
							at(read, read_pos) = 3;
						}
					}
				}

				if(var_bases_left){
					if(0 == --var_bases_left){
						++cur_var;
					}
				}
				if(0 == var_bases_left){
					++ref_pos;
				}
				++read_pos;
			}

			// Jump over region without information
			if( ++cur_meth < regions.size() && regions.at(cur_meth).first > ref_pos ){
				// Handle all variants before the next unmethylated region
				while( cur_var < variants.size() && regions.at(cur_meth).first > variants.at(cur_var).position_ && read_pos < length(read) ){
					if( variants.at(cur_var).InAllele(allele) ){
						read_pos += variants.at(cur_var).position_ - ref_pos;
						read_pos += length(variants.at(cur_var).var_seq_);
						ref_pos = variants.at(cur_var).position_+1; // Go one after the variant, because we handle the variant here, and it is only at exactly one ref_pos, so cannot reach into unmethylated region
					}
					++cur_var;
				}

				// Handle from the last variant before the region to the beginning of the region
				read_pos += regions.at(cur_meth).first - ref_pos;
				ref_pos = regions.at(cur_meth).first;
			}
		}
	}
}

void Simulator::CTConversion(
		SimPair &sim_reads,
		bool strand,
		uintAlleleId allele,
		uintSeqLen cur_start_position,
		uintSeqLen cur_end_position,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		const VariantBiasVarModifiers &bias_mod,
		intVariantId cur_methylation_start,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen ) const{
	if(ref.MethylationLoaded()){
		if(ref.VariantsLoaded()){
			// Forward
			CTConversion( sim_reads.org_seq_.at(strand), ref, ref_seq_id, cur_start_position, allele, cur_methylation_start, rdist, rgen, false, ref.Variants(ref_seq_id), bias_mod.StartVariant());

			// Reverse
			CTConversion( sim_reads.org_seq_.at(!strand), ref, ref_seq_id, cur_end_position, allele, cur_methylation_start, rdist, rgen, true, ref.Variants(ref_seq_id), bias_mod.EndVariant(ref.Variants(ref_seq_id), cur_end_position, allele));
		}
		else{
			// Forward
			CTConversion( sim_reads.org_seq_.at(strand), ref, ref_seq_id, cur_start_position, allele, cur_methylation_start, rdist, rgen, false);

			// Reverse
			CTConversion( sim_reads.org_seq_.at(!strand), ref, ref_seq_id, cur_end_position, allele, cur_methylation_start, rdist, rgen, true);
		}
	}
}

bool Simulator::SimulateFromGivenBlock(
		const SimBlock &block, // forward block
		const SimUnit &unit,
		const Reference &ref,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen){
	// Seed the random generator for this block
	rgen.seed( block.seed_ );

	SimPair sim_reads;

	uintFragCount read_number(0);

	VariantBiasVarModifiers bias_mod(block.first_variant_id_, ref.NumAlleles());

	uintDupCount fragment_counts;

	uintPercent gc_perc;

	uintSeqLen cur_end_position;
	Surrounding surrounding_start;

	intVariantId cur_methylation_start = block.first_methylation_id_;

	array<double, 2> probability_chosen;

	// Take the surrounding shifted by 1 as the first thing the loop does is shifting it back
	ref.ForwardSurrounding( surrounding_start, unit.ref_seq_id_, ( 0<block.start_pos_ ? block.start_pos_-1 : ref.SequenceLength(unit.ref_seq_id_)-1 ) );

	for( auto cur_start_position = block.start_pos_; cur_start_position < block.start_pos_ + kBlockSize && cur_start_position < ref.SequenceLength(unit.ref_seq_id_); ++cur_start_position ){
		surrounding_start.UpdateForward( ref.ReferenceSequence(unit.ref_seq_id_), cur_start_position );

		if(ref.MethylationLoaded()){
			if(cur_methylation_start < ref.UnmethylatedRegions(unit.ref_seq_id_).size() && ref.UnmethylatedRegions(unit.ref_seq_id_).at(cur_methylation_start).second <= cur_start_position){
				++cur_methylation_start;
			}
		}

		do{ //while(bias_mod.start_variant_pos_)
			auto frag_len_start = max(static_cast<uintSeqLen>(1),static_cast<uintSeqLen>(stats.FragmentDistribution().InsertLengths().from()));
			PrepareBiasModForCurrentStartPos( bias_mod, unit.ref_seq_id_, ref, cur_start_position, frag_len_start, surrounding_start );

			for( auto fragment_length=frag_len_start; fragment_length < stats.FragmentDistribution().InsertLengths().to(); ++fragment_length){
				for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
					if( !ref.VariantsLoaded() || !AlleleSkipped(bias_mod, allele, ref.Variants(unit.ref_seq_id_), cur_start_position) ){
						for( auto &prob : probability_chosen ){
							prob = rdist.ZeroToOne(rgen);
						}

						if(ProbabilityAboveThreshold(probability_chosen, unit.ref_seq_id_, fragment_length)){
							PrepareBiasModForCurrentFragmentLength( bias_mod, unit.ref_seq_id_, ref, cur_start_position, fragment_length, allele );

							cur_end_position = cur_start_position + fragment_length + bias_mod.end_pos_shift_.at(allele);
							if( cur_end_position < ref.SequenceLength(unit.ref_seq_id_)){
								for(bool strand : { false, true }){
									if(ProbabilityAboveThreshold(probability_chosen, strand, unit.ref_seq_id_, fragment_length)){
										// Determine how many read pairs are generated for this strand and allele at this position with this fragment_length
										gc_perc = GetGCPercent( bias_mod, unit.ref_seq_id_, ref, cur_end_position, fragment_length, allele );

										if(ref.VariantsLoaded()){
											fragment_counts = stats.FragmentDistribution().GetFragmentCounts(bias_normalization_, unit.ref_seq_id_, fragment_length, gc_perc, bias_mod.surrounding_start_.at(allele), bias_mod.surrounding_end_.at(allele), probability_chosen.at(strand));
										}
										else{
											fragment_counts = stats.FragmentDistribution().GetFragmentCounts(bias_normalization_, unit.ref_seq_id_, fragment_length, gc_perc, surrounding_start, bias_mod.surrounding_end_.at(allele), probability_chosen.at(strand));
										}

										if( fragment_counts ){
											GetOrgSeq(sim_reads, strand, allele, fragment_length, cur_start_position, cur_end_position, unit.ref_seq_id_, ref, stats, bias_mod);

											CTConversion(sim_reads, strand, allele, cur_start_position, cur_end_position, unit.ref_seq_id_, ref, bias_mod, cur_methylation_start, rdist, rgen);

											if(ref.VariantsLoaded()){
												if(!CreateReads( ref, stats, estimates, rdist, sim_reads, rgen, fragment_counts, strand, allele, unit.ref_seq_id_, fragment_length, read_number, &block, cur_start_position, cur_end_position, bias_mod.StartVariant(), bias_mod.EndVariant(ref.Variants(unit.ref_seq_id_), cur_end_position, allele) )){
													return false;
												}
											}
											else{
												if(!CreateReads( ref, stats, estimates, rdist, sim_reads, rgen, fragment_counts, strand, allele, unit.ref_seq_id_, fragment_length, read_number, &block, cur_start_position, cur_end_position )){
													return false;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			CheckForInsertedBasesToStartFrom(bias_mod, unit.ref_seq_id_, cur_start_position, ref);
		} while(bias_mod.start_variant_pos_);
	}

	return true;
}

bool Simulator::SimulateAdapterOnlyPairs(
		const Reference &ref,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen){
	if(num_adapter_only_pairs_){
		rdist.Reset();

		SimPair sim_reads;
		for(uintTempSeq template_segment=2; template_segment--; ){
			sim_reads.org_seq_.at(template_segment) = "";
		}

		rgen.seed( block_seed_gen_() );
		uintFragCount read_number(0);

		if( !CreateReads(ref, stats, estimates, rdist, sim_reads, rgen, num_adapter_only_pairs_, false, 0, 0, 0, read_number, NULL, 0, 0) ){
			return false;
		}
	}

	return true;
}

void Simulator::SimulationThread( Simulator &self, Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates ){
	mt19937_64 rgen;
	SimBlock *block;
	SimUnit *unit;
	GeneralRandomDistributions rdist(stats);

	while( !self.simulation_error_ && self.GetNextBlock(ref, stats, estimates, block, unit) ){
		// As long as there are new blocks to simulate from, simulate
		rdist.Reset();
		self.SimulateFromGivenBlock(*block, *unit, ref, stats, estimates, rdist, rgen);
		block->finished_ = true;
	}

	// Adapters are positioned after all blocks are created to avoid a race condition with the first block to guarantee reproducible seeds to each block and the adapters
	if(!self.simulation_error_ && !self.adapter_only_simulated_.test_and_set()){
		self.SimulateAdapterOnlyPairs(ref, stats, estimates, rdist, rgen);
	}
}

bool Simulator::WriteOutSystematicErrorProfile(const string &id, vector<pair<Dna5,uintPercent>> sys_errors, Dna5String dom_err, CharString err_perc){
	// Convert systematic errors into seq, qual for fastq format
	resize(dom_err, sys_errors.size());
	resize(err_perc, sys_errors.size());

	for(uintSeqLen pos=0; pos<sys_errors.size(); ++pos){
		at(dom_err, pos) = sys_errors.at(pos).first;
		uintPercent q_perc = sys_errors.at(pos).second;
		if(86 < q_perc){
			// Shorten percentage that it fits into 94 (0-93) values (odd values above 86 are joined with the even value before, so e.g. 87->86)
			q_perc -= (q_perc - 85)/2;
		}
		q_perc += 33; // Shift to readable ascii codes
		at(err_perc, pos) = q_perc;
	}

	// Write record
	try{
		writeRecord(dest_.at(0), id.c_str(), dom_err, err_perc);
	}
	catch(const Exception &e){
		printErr << "Could not write systematic error profile: " << e.what() << std::endl;
		return false;
	}

	return true;
}

Simulator::Simulator():
	written_records_(0),
	last_unit_(NULL),
	rdist_zero_to_one_(0,1)
	{
}

bool Simulator::CreateSystematicErrorProfile(
		const char *destination_file,
		const Reference &ref,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		uintSeed seed){

	if( ref.HasN() ){
		printErr << "Reference contains ambiguous bases(e.g. N). Please replace them, remove them from the scaffolds or split scaffolds into contigs." << std::endl;
		return false;
	}
	else{
		// create all missing directories
		CreateDir(destination_file);

		if( !open(dest_.at(0), destination_file) ){
			printErr << "Could not open '" << destination_file << "' for writing." << std::endl;
			return false;
		}
		else{
			block_seed_gen_.seed(seed);

			vector<pair<Dna5,uintPercent>> sys_errors;
			auto max_len = ref.MaxSequenceLength();
			sys_errors.reserve( max_len );
			Dna5String dom_err;
			reserve(dom_err, max_len);
			CharString err_perc;
			reserve(err_perc, max_len);

			for(uintRefSeqId ref_id=0; ref_id < ref.NumberSequences(); ++ref_id){
				// Get sytematic error for reverse strand (error is reversed, last error at first position, this is how it will be needed for simulation)
				ResetSystematicErrorCounters(ref);
				sys_errors.clear();
				ConstDna5StringReverseComplement reversed_ref( ref.ReferenceSequence(ref_id) );
				SetSystematicErrors(sys_errors, reversed_ref, 0, ref.SequenceLength(ref_id), stats, estimates);
				if( !WriteOutSystematicErrorProfile(string(toCString(ref.ReferenceId(ref_id))) + " reverse", sys_errors, dom_err, err_perc) ){
					return false;
				}

				// Get sytematic error for forward strand (is printed second, like it will be needed for simulation)
				ResetSystematicErrorCounters(ref);
				sys_errors.clear();
				SetSystematicErrors(sys_errors, ref.ReferenceSequence(ref_id), 0, ref.SequenceLength(ref_id), stats, estimates);
				if( !WriteOutSystematicErrorProfile(string(toCString(ref.ReferenceId(ref_id))) + " forward", sys_errors, dom_err, err_perc) ){
					return false;
				}
			}
		}
	}

	return true;
}

bool Simulator::Simulate(
		const char *destination_file_first,
		const char *destination_file_second,
		Reference &ref,
		DataStats &stats,
		const ProbabilityEstimates &estimates,
		uintNumThreads num_threads,
		uintSeed seed,
		uintFragCount num_read_pairs,
		double coverage,
		RefSeqBiasSimulation ref_bias_model,
		const string &ref_bias_file,
		const string &sys_error_file,
		const string &record_base_identifier,
		const string &var_file,
		const string &meth_file
		){
	// create all missing directories
	CreateDir(destination_file_first);
	CreateDir(destination_file_second);

	if( !open(dest_.at(0), destination_file_first) ){
		printErr << "Could not open '" << destination_file_first << "' for writing." << std::endl;
	}
	else if( !open(dest_.at(1), destination_file_second) ){
		printErr << "Could not open '" << destination_file_second << "' for writing." << std::endl;
		close(dest_.at(0));
	}
	else{
		// Set up random generator
		block_seed_gen_.seed(seed);

		// Replace N's in Reference to be consistent every time we simulate from a position
		ref.ReplaceN(seed);

		// Update the reference sequence bias so it fits to the potentially different reference
		if(stats.FragmentDistribution().UpdateRefSeqBias(ref_bias_model, ref_bias_file, ref, block_seed_gen_)){
			// Provide StringSets to store reads in before writing them out
			for( int i=2; i--; ){
				output_ids_.at(i) = new StringSet<CharString>;
				reserve(*output_ids_.at(i), kBatchSize, Exact());
				output_seqs_.at(i) = new StringSet<Dna5String>;
				reserve(*output_seqs_.at(i), kBatchSize, Exact());
				output_quals_.at(i) = new StringSet<CharString>;
				reserve(*output_quals_.at(i), kBatchSize, Exact());
			}

			// Store parameter in class
			if(record_base_identifier.size()){
				record_base_identifier_ = record_base_identifier;
			}
			else{
				record_base_identifier_ = "ReseqRead";
			}

			// Get average read length
			uintFragCount reads(0);
			uintRefLenCalc sum_read_length(0);
			for(auto template_segment=2; template_segment--;){
				for( auto len = stats.ReadLengths(template_segment).from(); len < stats.ReadLengths(template_segment).to(); ++len){
					reads += stats.ReadLengths(template_segment).at(len);
					sum_read_length += stats.ReadLengths(template_segment).at(len) * len;
				}
			}
			double average_read_length = static_cast<double>(sum_read_length)/reads;

			// Get number of read pairs to simulate
			auto total_ref_size = ref.TotalSize();
			double adapter_part = CoveragePropLostFromAdapters(stats);
			if(0 != num_read_pairs){
				total_pairs_ = num_read_pairs;
			}
			else{
				if(0.0 == coverage){
					// Keep the coverage from the original dataset
					coverage = stats.CorrectedCoverage();
				}

				total_pairs_ = CoverageToNumberPairs(coverage, total_ref_size, average_read_length, adapter_part);
			}

			// Keep the percentage of adapter only pairs
			num_adapter_only_pairs_ = round( static_cast<double>(total_pairs_)*stats.FragmentDistribution().InsertLengths()[0]/(stats.TotalNumberReads()/2) );

			printInfo << "Aiming for " << total_pairs_ << " read pairs, which corresponds to an approximated read depth of " << setprecision(4) << NumberPairsToCoverage(total_pairs_, total_ref_size, average_read_length, adapter_part) << setprecision(6) << "x" << std::endl;
			// Adapter only pairs(InsertLength == 0) are handled separately
			total_pairs_ -= num_adapter_only_pairs_;

			printInfo << "Preparing for simulation" << std::endl;
			simulation_error_ = false;

			// Prepare vcf file handle if needed
			if(!var_file.empty()){
				if(ref.PrepareVariantFile(var_file)){
					if(ref.ReadFirstVariants()){
						sys_dom_base_per_allele_.resize(ref.NumAlleles());
					}
					else{
						simulation_error_ = true;
					}
				}
				else{
					simulation_error_ = true;
				}
			}

			if(!simulation_error_ && !meth_file.empty()){
				if(ref.PrepareMethylationFile(meth_file)){
					if(!ref.ReadMethylation(2)){
						simulation_error_ = true;
					}
				}
				else{
					simulation_error_ = true;
				}
			}

			if(!simulation_error_){
				bias_normalization_ = stats.FragmentDistribution().CalculateBiasNormalization(coverage_groups_, non_zero_thresholds_, ref, num_threads, total_pairs_);

				if(0.0 == bias_normalization_){
					simulation_error_ = true;
				}
				else{
					// Chose systematic errors for adapters
					sys_gc_range_ = Divide(sum_read_length,reads)/2; // Set gc range for systematic errors to half the average read length

					for(auto template_segment=2; template_segment--;){
						adapter_sys_error_.at(template_segment).reserve( stats.Adapters().Counts(template_segment).size() );
						adapter_sys_error_.at(template_segment).resize( stats.Adapters().Counts(template_segment).size() );

						for(auto seq=stats.Adapters().Counts(template_segment).size(); seq--;){
							if( stats.Adapters().Counts(template_segment).at(seq) ){ // If the adapter does not appear, we don't need it
								ResetSystematicErrorCounters(ref);

								adapter_sys_error_.at(template_segment).at(seq).reserve( length(stats.Adapters().Sequence(template_segment, seq)) );

								SetSystematicErrors(adapter_sys_error_.at(template_segment).at(seq), stats.Adapters().Sequence(template_segment, seq), 0, length(stats.Adapters().Sequence(template_segment, seq)), stats, estimates);
							}
						}
					}

					// Prepare systematic error file for reading in case it is specified
					if(!sys_error_file.empty()){
						if( !open(sys_file_, sys_error_file.c_str()) ){
							printErr << "Could not open '" << sys_error_file << "' for reading." << std::endl;
							simulation_error_ = true;
						}
						else{
							sys_from_file_ = true;

							auto max_len = ref.MaxSequenceLength();
							reserve(sys_dom_err_, max_len);
							reserve(sys_err_perc_, max_len);
						}
					}
					else{
						sys_from_file_ = false;
					}
				}
			}

			if(!simulation_error_){
				adapter_only_simulated_.clear();

				// Create a buffer of blocks, so that the blocks into which long fragments reach already exist
				if( CreateBlock(ref, stats, estimates) ){
					for( auto n_blocks = stats.FragmentDistribution().InsertLengths().to()/kBlockSize; n_blocks--; ){
						CreateBlock(ref, stats, estimates);
					}

					printInfo << "Starting read generation" << std::endl;

					thread threads[num_threads];
					for(auto i = num_threads; i--; ){
						threads[i] = thread(SimulationThread, std::ref(*this), std::ref(ref), std::cref(stats), std::cref(estimates));
					}
					for(auto i = num_threads; i--; ){
						threads[i].join();
					}
				}
				else{
					simulation_error_ = true;
				}
			}

			if(simulation_error_){
				printErr << "An error occurred in the process: Terminating simulation" << std::endl;
			}

			this->Flush(); // Flush out all sequences that have not been written to disc yet

			// Clean up left over blocks and units
			if(first_unit_){
				SimBlock *del_block(first_unit_->first_block_);
				SimUnit *del_unit(first_unit_);
				while( del_block ){
					if(first_unit_->first_block_->next_block_){
						first_unit_->first_block_ = first_unit_->first_block_->next_block_;
					}
					else if(first_unit_->next_unit_){
						first_unit_ = first_unit_->next_unit_;
						delete del_unit;
						del_unit = first_unit_;
					}
					else{
						// Terminate loop
						first_unit_->first_block_ = NULL;
					}

					delete del_block->partner_block_;
					delete del_block;
					del_block = first_unit_->first_block_;
				}
				delete del_unit;
				first_unit_ = NULL;
				last_unit_ = NULL;
			}
		}

		close(dest_.at(0));
		close(dest_.at(1));
		if(sys_from_file_){
			close(sys_file_);
		}

		ref.ClearAllVariants();
		ref.ClearAllMethylation();
	}

	if(simulation_error_){
		// Remove the simulated files so it is clear that we encountered an error and do not accidentally continue with the old file
		DeleteFile( destination_file_first );
		DeleteFile( destination_file_second );
		return false;
	}
	else{
		printInfo << "Simulation finished succesfully" << std::endl;
		return true;
	}
}
