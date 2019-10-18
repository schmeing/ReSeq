#include "Simulator.h"
using reseq::Simulator;

#include <algorithm>
using std::max;
using std::min;
#include <array>
using std::array;
//include <atomic>
using std::atomic;
#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::steady_clock;
#include <cmath>
using std::round;
#include <exception>
using std::exception;
#include <iomanip>
using std::setprecision;
#include <iostream>
using std::cout;
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
using reseq::utilities::GetDominantLastX;
using reseq::utilities::IsGC;
using reseq::utilities::Percent;
using reseq::utilities::ReverseComplementorDna;
using reseq::utilities::SafePercent;
using reseq::utilities::SetDominantLastX;
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
					adapter_bases += stats.ReadLengthsByFragmentLength(seg).at(frag_len).at(read_len) * (read_len-frag_len);
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

				if( block->sys_errors_.size() <= ++block_pos ){
					block = block->next_block_;
					block_pos = 0;
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
	par.seq_qual_ = estimates.SequenceQuality(template_segment, tile_id).Draw(prob_seq_qual, prob_sum, {par.gc_seq_, mean_error_rate, fragment_length}, rdist.ZeroToOne(rgen));
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
	sys_dom_last5_ = 4;
	sys_seq_content_last5_.fill(0);

	if(ref.VariantsLoaded()){
		sys_last_var_pos_per_allele_.clear();
		sys_last_var_pos_per_allele_.resize(ref.NumAlleles(), numeric_limits<uintSeqLen>::max());
		sys_last_base_per_allele_.clear();
		sys_last_base_per_allele_.resize(ref.NumAlleles(), 4);
		sys_dom_last5_per_allele_.clear();
		sys_dom_last5_per_allele_.resize(ref.NumAlleles(), 4);
		sys_seq_content_last5_per_allele_.clear();
		sys_seq_content_last5_per_allele_.resize(ref.NumAlleles(), sys_seq_content_last5_);
		sys_seq_last5_per_allele_.clear();
		sys_seq_last5_per_allele_.resize(ref.NumAlleles());
	}

	sys_gc_ = 0;
	sys_gc_bases_ = 0;

	distance_to_start_of_error_region_ = 0;
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

void Simulator::SetSystematicErrorVariantsReverse(uintSeqLen &start_dist_error_region, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	if(ref.VariantsLoaded()){
		uintAlleleId chosen_allele;
		vector<pair<Dna5,uintPercent>> tmp_sys_errors;
		Dna5 last_base;
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
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele) <= var.position_+5 ){
				for( auto pos = sys_last_var_pos_per_allele_.at(chosen_allele)-1; pos > var.position_; --pos ){
					SetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos)), 5, sys_seq_last5_per_allele_.at(chosen_allele), length(sys_seq_last5_per_allele_.at(chosen_allele)) );
					if(length(sys_seq_last5_per_allele_.at(chosen_allele)) > 4){
						erase(sys_seq_last5_per_allele_.at(chosen_allele), 0);
					}
					sys_seq_last5_per_allele_.at(chosen_allele) += Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos));
				}
			}
			else{
				sys_seq_content_last5_per_allele_.at(chosen_allele).fill(0);
				if(rev_pos){
					GetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), 5, ConstDna5StringReverseComplement(ref.ReferenceSequence(ref_seq_id)), rev_pos );
					ref.ReferenceSequence(sys_seq_last5_per_allele_.at(chosen_allele), ref_seq_id, var.position_+min(rev_pos,static_cast<uintSeqLen>(5))+1, min(rev_pos,static_cast<uintSeqLen>(5)), true);
				}
				else{
					sys_dom_last5_per_allele_.at(chosen_allele) = 4;
					sys_seq_last5_per_allele_.at(chosen_allele) = "";
				}
			}

			// Get start_dist_error_region (Variants cannot start an error region)
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < block_pos; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = block.start_pos_+block.sys_errors_.size()-ref.Variants(ref_seq_id).at(var_id+1).position_-1; pos < block_pos; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
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
					DrawSystematicError(tmp_sys_errors, Complement::Dna(at(var.var_seq_, vpos)), last_base, sys_dom_last5_per_allele_.at(chosen_allele), SafePercent(gc, gc_bases), TransformDistanceToStartOfErrorRegion(start_dist_error_region), estimates);
					last_base = Complement::Dna(at(var.var_seq_, vpos));
					SetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), last_base, 5, sys_seq_last5_per_allele_.at(chosen_allele), length(sys_seq_last5_per_allele_.at(chosen_allele)));
					if(length(sys_seq_last5_per_allele_.at(chosen_allele)) > 4){
						erase(sys_seq_last5_per_allele_.at(chosen_allele), 0);
					}
					sys_seq_last5_per_allele_.at(chosen_allele) += last_base;
				}
			}

			// Insert variant
			block.err_variants_.emplace_back(block_pos, tmp_sys_errors, var.allele_);

			// Update information for next variants
			auto ref_allele = ref.NumAlleles(); // For this allele we already used the reference, so in case we have another allele identical to the reference up to the current variance we can copy the values
			if(sys_last_var_pos_per_allele_.at(chosen_allele) >= ref.SequenceLength(ref_seq_id) || sys_last_var_pos_per_allele_.at(chosen_allele)+length(var.var_seq_) > var.position_+5){
				ref_allele = chosen_allele;
			}
			for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
				if(var.InAllele(allele)){
					if( allele != chosen_allele){ // Dom base for chosen allele was already updated during the drawing process
						// Update dom base
						if(sys_last_var_pos_per_allele_.at(allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(allele)+length(var.var_seq_) <= var.position_+5){
							for( auto pos = sys_last_var_pos_per_allele_.at(allele)-1; pos > var.position_; --pos ){
								SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos)), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)) );
								if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
									erase(sys_seq_last5_per_allele_.at(allele), 0);
								}
								sys_seq_last5_per_allele_.at(allele) += Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), pos));
							}

							for(uintSeqLen vpos=length(var.var_seq_); vpos--; ){
								SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), Complement::Dna(at(var.var_seq_, vpos)), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)));
								if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
									erase(sys_seq_last5_per_allele_.at(allele), 0);
								}
								sys_seq_last5_per_allele_.at(allele) += Complement::Dna(at(var.var_seq_, vpos));
							}
						}
						else{
							if(ref_allele < ref.NumAlleles()){
								// We can copy from an identical allele
								sys_dom_last5_per_allele_.at(allele) = sys_dom_last5_per_allele_.at(ref_allele);
								for(auto base=sys_seq_content_last5_per_allele_.at(allele).size(); base--; ){
									sys_seq_content_last5_per_allele_.at(allele) = sys_seq_content_last5_per_allele_.at(ref_allele);
								}
								sys_seq_last5_per_allele_.at(allele) = sys_seq_last5_per_allele_.at(ref_allele);
							}
							else{
								ref_allele = allele;

								sys_seq_content_last5_per_allele_.at(allele).fill(0);
								if(rev_pos){
									GetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), 5, ConstDna5StringReverseComplement(ref.ReferenceSequence(ref_seq_id)), rev_pos );
									ref.ReferenceSequence(sys_seq_last5_per_allele_.at(allele), ref_seq_id, var.position_+min(rev_pos,static_cast<uintSeqLen>(5))+1, min(rev_pos,static_cast<uintSeqLen>(5)), true);
								}
								else{
									sys_dom_last5_per_allele_.at(allele) = 4;
									sys_seq_last5_per_allele_.at(allele) = "";
								}

								for(uintSeqLen vpos=length(var.var_seq_); vpos--; ){
									SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), Complement::Dna(at(var.var_seq_, vpos)), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)));
									if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
										erase(sys_seq_last5_per_allele_.at(allele), 0);
									}
									sys_seq_last5_per_allele_.at(allele) += Complement::Dna(at(var.var_seq_, vpos));
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
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = block.start_pos_+block.sys_errors_.size()-ref.Variants(ref_seq_id).at(var_id+1).position_-1; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
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

	// Set sytematic error for all reverse blocks of the unit
	ConstDna5StringReverseComplement reversed_ref( ref.ReferenceSequence(ref_id) );
	ResetSystematicErrorCounters(ref);
	uintSeqLen start_pos(0), end_pos(ref.SequenceLength(unit->ref_seq_id_) - block->start_pos_);
	while(block){
		block->sys_errors_.reserve(end_pos - start_pos);
		if(sys_from_file_){
			ReadSystematicErrors(block->sys_errors_, start_pos, end_pos);
			SetSystematicErrorVariantsReverse(distance_to_start_of_error_region_, *block, unit->ref_seq_id_, block->start_pos_, ref, stats, estimates);
		}
		else{
			auto tmp_distance_to_start_of_error_region = distance_to_start_of_error_region_;
			SetSystematicErrors(block->sys_errors_, reversed_ref, start_pos, end_pos, stats, estimates);
			SetSystematicErrorVariantsReverse(tmp_distance_to_start_of_error_region, *block, unit->ref_seq_id_, block->start_pos_, ref, stats, estimates);
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

void Simulator::SetSystematicErrorVariantsForward(uintSeqLen &start_dist_error_region, SimBlock &block, uintRefSeqId ref_seq_id, uintSeqLen end_pos, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	if(ref.VariantsLoaded()){
		uintAlleleId chosen_allele;
		vector<pair<Dna5,uintPercent>> tmp_sys_errors;
		Dna5 last_base;
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
			if( sys_last_var_pos_per_allele_.at(chosen_allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(chosen_allele)+5 >= var.position_ ){
				for( auto pos = sys_last_var_pos_per_allele_.at(chosen_allele)+1; pos < var.position_; ++pos ){
					SetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), at(ref.ReferenceSequence(ref_seq_id), pos), 5, sys_seq_last5_per_allele_.at(chosen_allele), length(sys_seq_last5_per_allele_.at(chosen_allele)) );
					if(length(sys_seq_last5_per_allele_.at(chosen_allele)) > 4){
						erase(sys_seq_last5_per_allele_.at(chosen_allele), 0);
					}
					sys_seq_last5_per_allele_.at(chosen_allele) += at(ref.ReferenceSequence(ref_seq_id), pos);
				}
			}
			else{
				sys_seq_content_last5_per_allele_.at(chosen_allele).fill(0);
				if(var.position_){
					GetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), 5, ref.ReferenceSequence(ref_seq_id), var.position_ );
					ref.ReferenceSequence(sys_seq_last5_per_allele_.at(chosen_allele), ref_seq_id, var.position_-min(var.position_,static_cast<uintSeqLen>(5)), min(var.position_,static_cast<uintSeqLen>(5)), false);
				}
				else{
					sys_dom_last5_per_allele_.at(chosen_allele) = 4;
					sys_seq_last5_per_allele_.at(chosen_allele) = "";
				}
			}

			// Get start_dist_error_region (Variants cannot start an error region)
			if( block.first_variant_id_ == var_id ){
				for(auto pos = 0; pos < var.position_-block.start_pos_; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = ref.Variants(ref_seq_id).at(var_id-1).position_-block.start_pos_; pos < var.position_-block.start_pos_; ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
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
					DrawSystematicError(tmp_sys_errors, at(var.var_seq_, vpos), last_base, sys_dom_last5_per_allele_.at(chosen_allele), SafePercent(gc, gc_bases), TransformDistanceToStartOfErrorRegion(start_dist_error_region), estimates);
					last_base = at(var.var_seq_, vpos);
					SetDominantLastX( sys_dom_last5_per_allele_.at(chosen_allele), sys_seq_content_last5_per_allele_.at(chosen_allele), last_base, 5, sys_seq_last5_per_allele_.at(chosen_allele), length(sys_seq_last5_per_allele_.at(chosen_allele)));
					if(length(sys_seq_last5_per_allele_.at(chosen_allele)) > 4){
						erase(sys_seq_last5_per_allele_.at(chosen_allele), 0);
					}
					sys_seq_last5_per_allele_.at(chosen_allele) += last_base;
				}
			}

			// Insert variant
			block.err_variants_.emplace_back(var.position_-block.start_pos_, tmp_sys_errors, var.allele_);

			// Update information for next variants
			auto ref_allele = ref.NumAlleles(); // For this allele we already used the reference, so in case we have another allele identical to the reference up to the current variance we can copy the values
			if(sys_last_var_pos_per_allele_.at(chosen_allele) >= ref.SequenceLength(ref_seq_id) || sys_last_var_pos_per_allele_.at(chosen_allele)+5 < var.position_+length(var.var_seq_)){
				ref_allele = chosen_allele;
			}
			for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
				if(var.InAllele(allele)){
					if( allele != chosen_allele){ // Dom base for chosen allele was already updated during the drawing process
						// Update dom base
						if(sys_last_var_pos_per_allele_.at(allele) < ref.SequenceLength(ref_seq_id) && sys_last_var_pos_per_allele_.at(allele)+5 >= var.position_+length(var.var_seq_)){
							for( auto pos = sys_last_var_pos_per_allele_.at(allele)+1; pos < var.position_; ++pos ){
								SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), at(ref.ReferenceSequence(ref_seq_id), pos), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)) );
								if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
									erase(sys_seq_last5_per_allele_.at(allele), 0);
								}
								sys_seq_last5_per_allele_.at(allele) += at(ref.ReferenceSequence(ref_seq_id), pos);
							}

							for(uintSeqLen vpos=0; vpos < length(var.var_seq_); ++vpos){
								SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), at(var.var_seq_, vpos), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)));
								if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
									erase(sys_seq_last5_per_allele_.at(allele), 0);
								}
								sys_seq_last5_per_allele_.at(allele) += at(var.var_seq_, vpos);
							}
						}
						else{
							if(ref_allele < ref.NumAlleles()){
								// We can copy from an identical allele
								sys_dom_last5_per_allele_.at(allele) = sys_dom_last5_per_allele_.at(ref_allele);
								for(auto base=sys_seq_content_last5_per_allele_.at(allele).size(); base--; ){
									sys_seq_content_last5_per_allele_.at(allele) = sys_seq_content_last5_per_allele_.at(ref_allele);
								}
								sys_seq_last5_per_allele_.at(allele) = sys_seq_last5_per_allele_.at(ref_allele);
							}
							else{
								ref_allele = allele;

								sys_seq_content_last5_per_allele_.at(allele).fill(0);
								if(var.position_){
									GetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), 5, ref.ReferenceSequence(ref_seq_id), var.position_ );
									ref.ReferenceSequence(sys_seq_last5_per_allele_.at(allele), ref_seq_id, var.position_-min(var.position_,static_cast<uintSeqLen>(5)), min(var.position_,static_cast<uintSeqLen>(5)), false);
								}
								else{
									sys_dom_last5_per_allele_.at(allele) = 4;
									sys_seq_last5_per_allele_.at(allele) = "";
								}

								for(uintSeqLen vpos=0; vpos < length(var.var_seq_); ++vpos){
									SetDominantLastX( sys_dom_last5_per_allele_.at(allele), sys_seq_content_last5_per_allele_.at(allele), at(var.var_seq_, vpos), 5, sys_seq_last5_per_allele_.at(allele), length(sys_seq_last5_per_allele_.at(allele)));
									if(length(sys_seq_last5_per_allele_.at(allele)) > 4){
										erase(sys_seq_last5_per_allele_.at(allele), 0);
									}
									sys_seq_last5_per_allele_.at(allele) += at(var.var_seq_, vpos);
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
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
				}
			}
			else{
				for(auto pos = ref.Variants(ref_seq_id).at(var_id-1).position_-block.start_pos_; pos < block.sys_errors_.size(); ++pos){
					stats.Coverage().UpdateDistances(start_dist_error_region, block.sys_errors_.at(pos).second);
				}
			}
		}
	}
}

bool Simulator::CreateBlock(Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, Benchmark &time){
	SimBlock *block;
	SimUnit *unit;

	steady_clock::time_point t1 = steady_clock::now();
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

			block = new SimBlock(last_unit_->last_block_->id_+1, last_unit_->last_block_->start_pos_ + kBlockSize, last_unit_->last_block_->partner_block_->partner_block_, block_seed_gen_());
			if(ref.VariantsLoaded()){
				// Find first variant that is in the new block (or one of the upcoming ones in case this block does not have variation), which is one after the last variant(reverse first) of the previous block
				block->first_variant_id_ = unit->last_block_->partner_block_->first_variant_id_+1;
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
			SetSystematicErrorVariantsForward(distance_to_start_of_error_region_, *block, unit->ref_seq_id_, end_pos, ref, stats, estimates);
		}
		else{
			auto tmp_distance_to_start_of_error_region = distance_to_start_of_error_region_;
			SetSystematicErrors(block->sys_errors_, ref.ReferenceSequence(unit->ref_seq_id_), block->start_pos_, end_pos, stats, estimates);
			SetSystematicErrorVariantsForward(tmp_distance_to_start_of_error_region, *block, unit->ref_seq_id_, end_pos, ref, stats, estimates);
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

		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		time.create_block_ += time_span.count();

		return true;
	}
	else{
		return false;
	}
}

bool Simulator::GetNextBlock( Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates, SimBlock *&block, SimUnit *&unit, Benchmark &time ){
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

	lock_guard<mutex> lock(block_creation_mutex_);

	// Create a new block to keep the buffer for long fragments
	if( !CreateBlock(ref, stats, estimates, time) ){
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
	while( VariantInsideCurrentFragment(bias_mod.unhandled_variant_id_.at(allele), cur_end_position, bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_, ref_seq_id, ref ) && 0==bias_mod.unhandled_bases_in_variant_.at(allele) ){
		auto &var( ref.Variants(ref_seq_id).at(bias_mod.unhandled_variant_id_.at(allele)) );

		if(!var.InAllele(allele)){
			++bias_mod.unhandled_variant_id_.at(allele);
		}
		else{
			// ---GC modification---
			// Add GC from variant
			for(uintSeqLen pos=0; pos < length(var.var_seq_) && var.position_ + pos < cur_end_position + bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_; ++pos){
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
				if(var.position_ + length(var.var_seq_) <= cur_end_position + bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_){
					// Completely in current fragment
					bias_mod.end_pos_shift_.at(allele) -= length(var.var_seq_) - 1;
					++bias_mod.unhandled_variant_id_.at(allele);
				}
				else{
					// Reaches out of fragment
					bias_mod.unhandled_bases_in_variant_.at(allele) = var.position_ + length(var.var_seq_) - (cur_end_position + bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_);
					bias_mod.end_pos_shift_.at(allele) -= length(var.var_seq_) - bias_mod.unhandled_bases_in_variant_.at(allele) - 1;
					// Variant is not handled completely so don't increase unhandled_variant_id
				}
			}
		}
	}
}

void Simulator::ChangeSurroundingBase(
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding,
		uintSurPos pos,
		Dna new_base){
	auto bit_in_block = 2*(Reference::surrounding_range_ - (pos%Reference::surrounding_range_) - 1);
	// ~ = bit-wise not
	surrounding.at(pos/Reference::surrounding_range_) = ( surrounding.at(pos/Reference::surrounding_range_) & ~(3 << bit_in_block) ) + ( static_cast<intSurrounding>(new_base) << bit_in_block );
}

void Simulator::DeleteSurroundingBaseShiftingOnRightSide(
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding,
		uintSurPos pos,
		Dna new_end_base){
	intSurrounding new_base = new_end_base;
	auto del_block = pos/Reference::surrounding_range_;

	// Add base at end and shift bases from end block to block with deletion
	for(auto block=Reference::num_surrounding_blocks_; --block > del_block; ){
		surrounding.at(block) = (surrounding.at(block) << 2) + new_base;
		new_base = surrounding.at(block) / Reference::SurroundingSize();
		surrounding.at(block) = surrounding.at(block) % Reference::SurroundingSize();
	}

	// Remove deleted base and add shifted base
	auto bit_in_block = 2*(Reference::surrounding_range_ - (pos%Reference::surrounding_range_) - 1);
	new_base += ( surrounding.at(del_block) % (1 << bit_in_block) ) << 2;
	surrounding.at(del_block) = ( surrounding.at(del_block) >> (bit_in_block+2) << (bit_in_block+2) ) + new_base;
}

void Simulator::DeleteSurroundingBaseShiftingOnLeftSide(
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding,
		uintSurPos pos,
		Dna new_end_base){
	intSurrounding new_base = new_end_base;
	auto del_block = pos/Reference::surrounding_range_;

	// Add base in front and shift bases from first block to block with deletion
	for(uintSurBlockId block=0; block < del_block; ++block){
		surrounding.at(block) += new_base*Reference::SurroundingSize();
		new_base = surrounding.at(block)%4;
		surrounding.at(block) = surrounding.at(block) >> 2;
	}

	// Add shifted base and remove deleted base
	auto bit_in_block = 2*(Reference::surrounding_range_ - (pos%Reference::surrounding_range_) - 1);
	surrounding.at(del_block) += new_base*Reference::SurroundingSize();
	surrounding.at(del_block) = ( surrounding.at(del_block) >> (bit_in_block+2) << (bit_in_block) ) + surrounding.at(del_block) % (1 << bit_in_block);
}

void Simulator::InsertSurroundingBasesShiftingOnRightSide(
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding,
		uintSurPos pos,
		DnaString new_bases){
	auto block = pos/Reference::surrounding_range_;
	uintSurPos bases_to_insert = min(length(new_bases), static_cast<size_t>(Reference::surrounding_range_*Reference::num_surrounding_blocks_-pos));

	// Shift right from current block to make room for new bases
	auto shift_blocks = bases_to_insert/Reference::surrounding_range_; // For every Reference::surrounding_range_ in bases_to_insert we can shift whole blocks
	auto shift_bases = bases_to_insert%Reference::surrounding_range_;
	for(auto cur_block=Reference::num_surrounding_blocks_-shift_blocks; cur_block-- > block+1; ){
		// First shift only the bases and ignore the blocks at the end that we completely remove
		surrounding.at(cur_block) >>= 2*shift_bases; // Remove bases not longer needed
		surrounding.at(cur_block) += surrounding.at(cur_block-1) % (1 << 2*shift_bases) * (1 << 2*(Reference::surrounding_range_-shift_bases)); // Add new bases from previous block
	}

	uintSurPos inv_pos_in_block = Reference::surrounding_range_ - (pos%Reference::surrounding_range_);
	intSurrounding tmp_sur = surrounding.at(block) % (1 << 2*inv_pos_in_block); // Store everything that needs to be shifted out of the start block
	surrounding.at(block) >>= 2*inv_pos_in_block; // And remove it from the block
	tmp_sur >>= 2*shift_bases; // Remove what was already shifted

	if(0 < shift_blocks && Reference::num_surrounding_blocks_ > block+shift_blocks){
		for(auto cur_block=Reference::num_surrounding_blocks_; cur_block-- > block+shift_blocks+1; ){
			// Now shift the blocks
			surrounding.at(cur_block) = surrounding.at(cur_block-shift_blocks);
		}

		// Shift everything that has to be shifted right from starting block and was not already shifted above (Happens if the bases need to go into an empty block freed by the block shift)
		surrounding.at(block+shift_blocks) = tmp_sur;
	}

	// Insert new bases in starting block
	uintSurPos ins_pos=0;
	uintSurPos bases_to_insert_into_this_block = min(bases_to_insert, inv_pos_in_block);
	bases_to_insert -= bases_to_insert_into_this_block;
	for(; bases_to_insert_into_this_block--; ){
		surrounding.at(block) <<= 2;
		surrounding.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos++));
	}

	if(0 == shift_blocks && inv_pos_in_block > shift_bases){
		// If the complete insertion happens in the starting block add the stuff back that was shifted out too much
		surrounding.at(block) <<= 2*(inv_pos_in_block-shift_bases);
		surrounding.at(block) += tmp_sur;
	}
	else{
		while(0 < bases_to_insert){
			bases_to_insert_into_this_block = min(bases_to_insert, Reference::surrounding_range_);
			bases_to_insert -= bases_to_insert_into_this_block;

			// Store the part of the block that is not overwritten by insertions
			tmp_sur = surrounding.at(++block) % (1 << 2*(Reference::surrounding_range_-bases_to_insert_into_this_block));

			// Insert new bases
			surrounding.at(block) = 0;
			for(auto i=bases_to_insert_into_this_block; i--; ){
				surrounding.at(block) <<= 2;
				surrounding.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos++));
			}

			// Restore the bases that should not be overwritten
			surrounding.at(block) <<= 2*(Reference::surrounding_range_-bases_to_insert_into_this_block);
			surrounding.at(block) += tmp_sur;
		}
	}
}

void Simulator::InsertSurroundingBasesShiftingOnLeftSide(
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding,
		uintSurPos pos,
		DnaString new_bases){
	auto block = pos/Reference::surrounding_range_;
	uintSurPos bases_to_insert = min(length(new_bases), static_cast<size_t>(pos+1));

	// Shift left from current block to make room for new bases
	auto shift_blocks = bases_to_insert/Reference::surrounding_range_; // For every Reference::surrounding_range_ in bases_to_insert we can shift whole blocks
	auto shift_bases = bases_to_insert%Reference::surrounding_range_;
	for(uintSurBlockId cur_block=shift_blocks; cur_block < block; ++cur_block ){
		// First shift only the bases and ignore the blocks at the beginning that we completely remove
		surrounding.at(cur_block) %= 1 << 2*(Reference::surrounding_range_-shift_bases);// Remove bases not longer needed
		surrounding.at(cur_block) <<= 2*shift_bases; // Shift remaining bases (do this after the removal to avoid overflow)
		surrounding.at(cur_block) += surrounding.at(cur_block+1) >> 2*(Reference::surrounding_range_-shift_bases); // Add new bases from next block
	}

	uintSurPos pos_in_block = pos%Reference::surrounding_range_+1; // 1-bases position in block
	intSurrounding tmp_sur = surrounding.at(block) % (1 << 2*(Reference::surrounding_range_-pos_in_block)); // Store everything that needs to stay in the start block
	if(pos_in_block > shift_bases){
		// Keep the stuff that is after the shift on the left before the insertion
		surrounding.at(block) >>= 2*(Reference::surrounding_range_-pos_in_block);
		surrounding.at(block) %= 1 << 2*(pos_in_block-shift_bases); // Remove what was already shifted
	}
	else{
		surrounding.at(block) = 0; // If everything was already shifted, we only need what's right of the insertion and that is in tmp_sur
	}

	if(0 < shift_blocks){
		if(block >= shift_blocks){
			for(uintSurBlockId cur_block=0; cur_block+shift_blocks < block; ++cur_block){
				// Now shift the blocks
				surrounding.at(cur_block) = surrounding.at(cur_block+shift_blocks);
			}

			// Shift everything that has to be shifted left from starting block and was not already shifted above (Happens if the bases need to go into an empty block freed by the block shift)
			surrounding.at(block-shift_blocks) = surrounding.at(block) << 2*(Reference::surrounding_range_-(pos_in_block-shift_bases));
		}
		surrounding.at(block) = 0; // The not-shifted stuff to keep was moved to the proper shifted block, so now we can clean the start block surrounding
	}

	// Insert new bases in starting block
	uintSurPos bases_to_insert_into_this_block = min(bases_to_insert, pos_in_block);
	uintSurPos ins_pos_to = length(new_bases); // In case we do not insert all bases from new_bases, we need to keep track of this separately to bases_to_insert
	for(uintSurPos ins_pos=ins_pos_to-bases_to_insert_into_this_block; ins_pos<ins_pos_to; ++ins_pos){
		surrounding.at(block) <<= 2;
		surrounding.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos));
	}
	bases_to_insert -= bases_to_insert_into_this_block;
	ins_pos_to -= bases_to_insert_into_this_block;

	// Add back what was right of insertion
	surrounding.at(block) <<= 2*(Reference::surrounding_range_-pos_in_block);
	surrounding.at(block) += tmp_sur;

	while(0 < bases_to_insert){
		bases_to_insert_into_this_block = min(bases_to_insert, Reference::surrounding_range_);

		// Leave only the part of the block that should not be overwritten by insertions
		surrounding.at(--block) >>= 2*bases_to_insert_into_this_block;

		// Insert new bases
		for(uintSurPos ins_pos=ins_pos_to-bases_to_insert_into_this_block; ins_pos<ins_pos_to; ++ins_pos){
			surrounding.at(block) <<= 2;
			surrounding.at(block) += static_cast<uintBaseCall>(at(new_bases, ins_pos));
		}
		bases_to_insert -= bases_to_insert_into_this_block;
		ins_pos_to -= bases_to_insert_into_this_block;
	}
}

void Simulator::HandleSurroundingVariantsBeforeCenter(
		array<intSurrounding, Reference::num_surrounding_blocks_> &mod_surrounding,
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

		if( 0 > sur_pos || Reference::num_surrounding_blocks_*Reference::surrounding_range_ <= sur_pos ){
			break;
		}

		if(ref.Variants(ref_seq_id).at(cur_var).InAllele(allele)){
			if(0 == length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
				// Deletion
				if(reverse){
					auto new_base_pos = static_cast<intSeqShift>(center_position) + pos_shift - static_cast<intSeqShift>(Reference::num_surrounding_blocks_ * Reference::surrounding_range_);
					if(0 > new_base_pos){
						DeleteSurroundingBaseShiftingOnRightSide(mod_surrounding, sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), ref.SequenceLength(ref_seq_id) + new_base_pos)));
					}
					else{
						DeleteSurroundingBaseShiftingOnRightSide(mod_surrounding, sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), new_base_pos)));
					}
					--pos_shift;
				}
				else{
					if(++pos_shift > center_position){
						DeleteSurroundingBaseShiftingOnLeftSide(mod_surrounding, sur_pos, at(ref.ReferenceSequence(ref_seq_id), ref.SequenceLength(ref_seq_id) + center_position - pos_shift));
					}
					else{
						DeleteSurroundingBaseShiftingOnLeftSide(mod_surrounding, sur_pos, at(ref.ReferenceSequence(ref_seq_id), center_position - pos_shift));
					}
				}
			}
			else{
				// Base modification (Insertions may also include a base modification and we skip checking first if we exchange it with the same)
				if(reverse){
					ChangeSurroundingBase(mod_surrounding, sur_pos, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)));
				}
				else{
					ChangeSurroundingBase(mod_surrounding, sur_pos, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0));
				}

				if(1 < length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
					// Insertion (First base is the already handled)
					if(reverse){
						InsertSurroundingBasesShiftingOnRightSide(mod_surrounding, sur_pos, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1)));
						pos_shift += length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
					}
					else{
						InsertSurroundingBasesShiftingOnLeftSide(mod_surrounding, sur_pos, suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1));
						pos_shift -= length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
					}
				}
			}
		}
	}
}

void Simulator::HandleSurroundingVariantsAfterCenter(
		array<intSurrounding, Reference::num_surrounding_blocks_> &mod_surrounding,
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

		if( 0 > sur_pos || Reference::num_surrounding_blocks_*Reference::surrounding_range_ <= sur_pos ){
			break;
		}

		if(ref.Variants(ref_seq_id).at(cur_var).InAllele(allele)){
			if(0 == length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
				// Deletion
				if(reverse){
					++pos_shift;
					DeleteSurroundingBaseShiftingOnLeftSide(mod_surrounding, sur_pos, Complement::Dna(at(ref.ReferenceSequence(ref_seq_id), (center_position + pos_shift)%ref.SequenceLength(ref_seq_id))));
				}
				else{
					DeleteSurroundingBaseShiftingOnRightSide(mod_surrounding, sur_pos, at(ref.ReferenceSequence(ref_seq_id), (center_position - pos_shift + Reference::num_surrounding_blocks_ * Reference::surrounding_range_)%ref.SequenceLength(ref_seq_id)));
					--pos_shift;
				}
			}
			else{
				// Base modification (Insertions may also include a base modification and we skip checking first if we exchange it with the same)
				if(reverse){
					ChangeSurroundingBase(mod_surrounding, sur_pos, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)));
				}
				else{
					ChangeSurroundingBase(mod_surrounding, sur_pos, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0));
				}

				if(1 < length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)){
					// Insertion (First base is the already handled)
					if(reverse){
						if(sur_pos){
							InsertSurroundingBasesShiftingOnLeftSide(mod_surrounding, sur_pos-1, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1)));
							pos_shift -= length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-1;
						}
					}
					else{
						if(sur_pos+1 < Reference::num_surrounding_blocks_*Reference::surrounding_range_){
							InsertSurroundingBasesShiftingOnRightSide(mod_surrounding, sur_pos+1, suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1));
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
		const array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_start) const{
	// Variants in the roll around are ignored for variant modifications of surrounding (but to fill in deletions still use the proper roll around base)
	// Copy reference surrounding into allele
	for(auto block=surrounding_start.size(); block--; ){
		bias_mod.mod_surrounding_start_.at(allele).at(block) = surrounding_start.at(block);
	}

	if(ref.Variants(ref_seq_id).size()){
		// Handle variants before center (Shifting to/from the left)
		intSeqShift pos_shift = -Reference::surrounding_start_pos_;
		auto cur_var = bias_mod.first_variant_id_;
		if(bias_mod.start_variant_pos_){
			// Handle the base substitution (as the insertion might include it)
			ChangeSurroundingBase(bias_mod.mod_surrounding_start_.at(allele), pos_shift, at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0) );
			// Partial insertion (as the rest is shifted to the right)
			InsertSurroundingBasesShiftingOnLeftSide(bias_mod.mod_surrounding_start_.at(allele), pos_shift, infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, bias_mod.start_variant_pos_+1));
			pos_shift -= bias_mod.start_variant_pos_;
		}

		HandleSurroundingVariantsBeforeCenter(bias_mod.mod_surrounding_start_.at(allele), cur_start_position, pos_shift, cur_var, allele, ref_seq_id, ref, false);

		// Handle variants after center (Shifting to/from the right)
		pos_shift = -Reference::surrounding_start_pos_;
		cur_var = bias_mod.first_variant_id_;
		if(bias_mod.start_variant_pos_){
			if(length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) > bias_mod.start_variant_pos_+1){
				// Partial insertion (as the rest is shifted to the left)
				if(pos_shift+1 < Reference::num_surrounding_blocks_*Reference::surrounding_range_){
					InsertSurroundingBasesShiftingOnRightSide(bias_mod.mod_surrounding_start_.at(allele), pos_shift+1, suffix(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, bias_mod.start_variant_pos_+1));
					pos_shift += length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) - bias_mod.start_variant_pos_ - 1;
				}
			}

			// Skip this insertion for the subsequent variant handling
			++cur_var;
		}

		HandleSurroundingVariantsAfterCenter(bias_mod.mod_surrounding_start_.at(allele), cur_start_position, pos_shift, cur_var, allele, ref_seq_id, ref, false);
	}
}

void Simulator::PrepareBiasModForCurrentStartPos(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		uintSeqLen cur_end_position,
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_start) const{
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
			uintSeqLen tmp_gc=0;
			for(uintSeqLen pos=bias_mod.start_variant_pos_; pos < length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ + pos - bias_mod.start_variant_pos_ < cur_end_position; ++pos){
				if( IsGC(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, pos) ){
					++tmp_gc;
				}
			}

			// Remove GC from reference as we are starting from the non-reference part
			if( ref.GC(ref_seq_id, cur_start_position) ){
				--tmp_gc;
			}

			// End pos shift cannot be longer than insert length or variant length
			auto tmp_shift = static_cast<intSeqShift>(1) - static_cast<intSeqShift>( min(cur_end_position-cur_start_position, static_cast<uintSeqLen>(length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_))-bias_mod.start_variant_pos_) );

			// Insert into all alleles (Alleles not having this insertion are anyways not included later)
			for(uintAlleleId all=0; all < ref.NumAlleles(); ++all){
				bias_mod.gc_mod_.at(all) = tmp_gc;
				bias_mod.end_pos_shift_.at(all) = tmp_shift;
				++bias_mod.unhandled_variant_id_.at(all); // Handled first_variant_id_
			}
		}

		// Check for every allele if variants cause changes to gc or surrounding
		for(uintAlleleId all=0; all < ref.NumAlleles(); ++all){
			if( !AlleleSkipped(bias_mod, all, ref.Variants(ref_seq_id), cur_start_position) ){
				VariantModStartSurrounding(bias_mod, all, ref_seq_id, ref, cur_start_position, surrounding_start); // Must be before HandleGCModAndEndPosShiftForNewVariants, so that bias_mod.unhandled_variant_id_ is untouched
				HandleGCModAndEndPosShiftForNewVariants(bias_mod, all, ref_seq_id, ref, cur_end_position);
			}
		}
	}
}

void Simulator::VariantModEndSurrounding(
		VariantBiasVarModifiers &bias_mod,
		uintAlleleId allele,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_end_position,
		const array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_end) const{
	// Variants in the roll around are ignored for variant modifications of surrounding (but to fill in deletions still use the proper roll around base)
	// Copy reference surrounding into allele if it doesn't need to much shifting afterwards
	auto last_position = cur_end_position - 1 + bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_;
	if( last_position < ref.SequenceLength(ref_seq_id) ){
		if( -5 < bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_ && bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_ < 5 ){
			for(auto block=surrounding_end.size(); block--; ){
				bias_mod.mod_surrounding_end_.at(allele).at(block) = surrounding_end.at(block);
			}

			if(0 < bias_mod.end_pos_shift_.at(allele) + bias_mod.fragment_length_extension_){
				for( auto up_pos = cur_end_position; up_pos <= last_position; ++up_pos ){
					ref.UpdateReverseSurrounding( bias_mod.mod_surrounding_end_.at(allele), ref.ReferenceSequence(ref_seq_id), up_pos );
				}
			}
			else{
				for( auto up_pos = cur_end_position-1; up_pos-- > last_position; ){
					ref.RollBackReverseSurrounding( bias_mod.mod_surrounding_end_.at(allele), ref.ReferenceSequence(ref_seq_id), up_pos );
				}
			}
		}
		else{
			ref.ReverseSurrounding( bias_mod.mod_surrounding_end_.at(allele), ref_seq_id, last_position );
		}

		if(ref.Variants(ref_seq_id).size()){
			// Handle variants after center (Shifting to/from the left)
			intSeqShift pos_shift = -Reference::surrounding_start_pos_;
			auto cur_var = bias_mod.unhandled_variant_id_.at(allele); // At the variants where pos == last_position if exists or after it otherwise
			if( bias_mod.unhandled_bases_in_variant_.at(allele) ){
				// Partial insertion (as the rest is shifted to the right)
				InsertSurroundingBasesShiftingOnLeftSide(bias_mod.mod_surrounding_end_.at(allele), pos_shift, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele))));
				pos_shift -= bias_mod.unhandled_bases_in_variant_.at(allele);
				++cur_var; // This insertion is handled, do not handle at again in HandleSurroundingVariantsAfterCenter
			}
			else if( bias_mod.start_variant_pos_ && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == last_position && length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) > bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1){
				if(pos_shift){
					// Partial insertion (as the rest is shifted to the right)
					InsertSurroundingBasesShiftingOnLeftSide(bias_mod.mod_surrounding_end_.at(allele), pos_shift-1, ReverseComplementorDna(suffix(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1)));
					pos_shift -= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_) - (bias_mod.start_variant_pos_ - bias_mod.end_pos_shift_.at(allele) + 1);
				}
			}

			HandleSurroundingVariantsAfterCenter(bias_mod.mod_surrounding_end_.at(allele), last_position, pos_shift, cur_var, allele, ref_seq_id, ref, true);

			// Handle variants before center (Shifting to/from the right)
			pos_shift = -Reference::surrounding_start_pos_;
			cur_var = bias_mod.unhandled_variant_id_.at(allele); // At the variants where pos == last_position if exists or after it otherwise
			if( bias_mod.unhandled_bases_in_variant_.at(allele) ){
				if(pos_shift+1 < Reference::num_surrounding_blocks_*Reference::surrounding_range_){
					// Handle the base substitution (as the insertion might include it)
					ChangeSurroundingBase(bias_mod.mod_surrounding_end_.at(allele), pos_shift+1, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)) );

					// Partial insertion (as the rest is shifted to the left)
					InsertSurroundingBasesShiftingOnRightSide(bias_mod.mod_surrounding_end_.at(allele), pos_shift+1, ReverseComplementorDna(infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele))));
					pos_shift += length(ref.Variants(ref_seq_id).at(cur_var).var_seq_)-bias_mod.unhandled_bases_in_variant_.at(allele)-1;
				}
			}
			else if( bias_mod.start_variant_pos_ && ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).position_ == last_position ){
				cur_var = bias_mod.first_variant_id_;
				// Handle the base substitution (as the insertion might include it)
				ChangeSurroundingBase(bias_mod.mod_surrounding_end_.at(allele), pos_shift, Complement::Dna(at(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 0)) );

				// Partial insertion (as the rest is shifted to the left)
				InsertSurroundingBasesShiftingOnRightSide(bias_mod.mod_surrounding_end_.at(allele), pos_shift, ReverseComplementorDna(infix(ref.Variants(ref_seq_id).at(cur_var).var_seq_, 1, bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele)+1)));
				pos_shift += bias_mod.start_variant_pos_-bias_mod.end_pos_shift_.at(allele);
			}

			HandleSurroundingVariantsBeforeCenter(bias_mod.mod_surrounding_end_.at(allele), last_position, pos_shift, cur_var, allele, ref_seq_id, ref, true);
		}
	}
}

void Simulator::UpdateBiasModForCurrentFragmentLength(
		VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		const Reference &ref,
		uintSeqLen cur_start_position,
		uintSeqLen cur_end_position,
		array<intSurrounding, Reference::num_surrounding_blocks_> &surrounding_end) const{
	if(ref.VariantsLoaded()){
		if(bias_mod.start_variant_pos_ && cur_end_position-cur_start_position <= length(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_)-bias_mod.start_variant_pos_){
			// Still in start variant insertion
			auto pos = bias_mod.start_variant_pos_ + cur_end_position-cur_start_position - 1;
			if( IsGC(ref.Variants(ref_seq_id).at(bias_mod.first_variant_id_).var_seq_, pos) ){
				for(uintAlleleId all=0; all < ref.NumAlleles(); ++all){
					++bias_mod.gc_mod_.at(all);
				}
			}

			for(uintAlleleId all=0; all < ref.NumAlleles(); ++all){
				if( !AlleleSkipped(bias_mod, all, ref.Variants(ref_seq_id), cur_start_position) ){
					--bias_mod.end_pos_shift_.at(all);
					VariantModEndSurrounding(bias_mod, all, ref_seq_id, ref, cur_end_position, surrounding_end);
				}
			}
		}
		else{
			for(uintAlleleId all=0; all < ref.NumAlleles(); ++all){
				if( !AlleleSkipped(bias_mod, all, ref.Variants(ref_seq_id), cur_start_position) ){
					VariantModEndSurrounding(bias_mod, all, ref_seq_id, ref, cur_end_position, surrounding_end); // Must be called before HandleGCModAndEndPosShiftForNewVariants, so that bias_mod.unhandled_variant_id_ and bias_mod.unhandled_bases_in_variant_ have not been updated yet

					if(bias_mod.unhandled_bases_in_variant_.at(all)){
						// Unhandled bases only occur for insertions and the reference base has already been covered
						auto &var( ref.Variants(ref_seq_id).at(bias_mod.unhandled_variant_id_.at(all)) );

						// ---GC modification---
						// Add GC from variant
						uintSeqLen pos = length(var.var_seq_)-bias_mod.unhandled_bases_in_variant_.at(all);
						if( IsGC(var.var_seq_, pos) ){
							++bias_mod.gc_mod_.at(all);
						}

						// ---End pos shift---
						--bias_mod.end_pos_shift_.at(all);

						if(0 == --bias_mod.unhandled_bases_in_variant_.at(all)){
							++bias_mod.unhandled_variant_id_.at(all);
						}
					}
					else{
						HandleGCModAndEndPosShiftForNewVariants(bias_mod, all, ref_seq_id, ref, cur_end_position);
					}
				}
			}
		}
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

void Simulator::CheckForFragmentLengthExtension( VariantBiasVarModifiers &bias_mod, uintSeqLen fragment_length, uintSeqLen fragment_length_to, const Reference &ref, const DataStats &stats ) const{
	if(ref.VariantsLoaded() && fragment_length+1 == fragment_length_to && fragment_length+1 < stats.FragmentDistribution().InsertLengths().to()){
		// Fragment length is limited by Reference Sequence end, so we can extend it for variants with a negative end_position_shift
		auto max_end_shift = - *min_element(bias_mod.end_pos_shift_.begin(), bias_mod.end_pos_shift_.end());
		SetToMax(max_end_shift, 0); // Shifts that reduce the maximum fragment_length don't need extra attention

		if(bias_mod.fragment_length_extension_ < max_end_shift && fragment_length+bias_mod.fragment_length_extension_+1 < stats.FragmentDistribution().InsertLengths().to()){
			++bias_mod.fragment_length_extension_;
		}
		else{
			bias_mod.fragment_length_extension_ = 0;
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
		const VariantBiasVarModifiers &bias_mod) const{
	if(ref.VariantsLoaded()){
		// Forward
		ref.ReferenceSequence( sim_reads.org_seq_.at(strand), ref_seq_id, cur_start_position, min(fragment_length+bias_mod.fragment_length_extension_, static_cast<uintSeqLen>(stats.ReadLengths(strand).to()+stats.Errors().MaxLenDeletion())), false, ref.Variants(ref_seq_id), bias_mod.StartVariant(), allele );

		// Reverse
		ref.ReferenceSequence( sim_reads.org_seq_.at(!strand), ref_seq_id, cur_end_position+bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_, min(fragment_length+bias_mod.fragment_length_extension_, static_cast<uintSeqLen>(stats.ReadLengths(!strand).to()+stats.Errors().MaxLenDeletion())), true, ref.Variants(ref_seq_id), bias_mod.EndVariant(ref.Variants(ref_seq_id), allele, cur_end_position), allele );
	}
	else{
		// Forward
		ref.ReferenceSequence( sim_reads.org_seq_.at(strand), ref_seq_id, cur_start_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(strand).to()+stats.Errors().MaxLenDeletion())), false );

		// Reverse
		ref.ReferenceSequence( sim_reads.org_seq_.at(!strand), ref_seq_id, cur_end_position, min(fragment_length, static_cast<uintSeqLen>(stats.ReadLengths(!strand).to()+stats.Errors().MaxLenDeletion())), true );
	}
}

bool Simulator::SimulateFromGivenBlock(
		const SimBlock &block, // forward block
		const SimUnit &unit,
		const Reference &ref,
		const DataStats &stats,
		const ProbabilityEstimates &estimates,
		GeneralRandomDistributions &rdist,
		mt19937_64 &rgen,
		Benchmark &time){
	steady_clock::time_point t1 = steady_clock::now();

	// Seed the random generator for this block
	rgen.seed( block.seed_ );

	SimPair sim_reads;

	uintFragCount read_number(0);

	VariantBiasVarModifiers bias_mod(block.first_variant_id_, (ref.VariantsLoaded() ? ref.NumAlleles() : 0));

	uintDupCount fragment_counts;

	uintSeqLen gc, tmp_gc(0);

	uintSeqLen cur_end_position, fragment_length_to;
	array<intSurrounding, Reference::num_surrounding_blocks_> surrounding_start, surrounding_end;

	// Take the surrounding shifted by 1 as the first thing the loop does is shifting it back
	ref.ForwardSurrounding( surrounding_start, unit.ref_seq_id_, ( 0<block.start_pos_ ? block.start_pos_-1 : ref.SequenceLength(unit.ref_seq_id_)-1 ) );

	for( auto cur_start_position = block.start_pos_; cur_start_position < block.start_pos_ + kBlockSize && cur_start_position < ref.SequenceLength(unit.ref_seq_id_); ++cur_start_position ){
		fragment_length_to = min(static_cast<uintSeqLen>(stats.FragmentDistribution().InsertLengths().to()), ref.SequenceLength(unit.ref_seq_id_)-cur_start_position+1); //+1 because it is one after the last acceptable fragment length

		ref.UpdateForwardSurrounding( surrounding_start, ref.ReferenceSequence(unit.ref_seq_id_), cur_start_position );

		do{ //while(bias_mod.start_variant_pos_)
			cur_end_position = max(static_cast<uintSeqLen>(1),static_cast<uintSeqLen>(stats.FragmentDistribution().InsertLengths().from()))-1 + cur_start_position;
			if( cur_end_position < ref.SequenceLength(unit.ref_seq_id_) || ref.VariantsLoaded()){ // Don't check gc if it's unnecessary (and will crash)
				if(cur_end_position >= ref.SequenceLength(unit.ref_seq_id_) && ref.VariantsLoaded()){
					bias_mod.fragment_length_extension_ = cur_end_position-ref.SequenceLength(unit.ref_seq_id_)+1;
					cur_end_position = ref.SequenceLength(unit.ref_seq_id_)-1;
					fragment_length_to -= bias_mod.fragment_length_extension_;
				}

				gc=ref.GCContentAbsolut( unit.ref_seq_id_, cur_start_position, cur_end_position );

				ref.ReverseSurrounding( surrounding_end, unit.ref_seq_id_, ( 0<cur_end_position ? cur_end_position-1 : ref.SequenceLength(unit.ref_seq_id_)-1 ) ); // Take the surrounding shifted by 1 as the first thing the loop does is shifting it back

				PrepareBiasModForCurrentStartPos( bias_mod, unit.ref_seq_id_, ref, cur_start_position, cur_end_position, surrounding_start );

				for( auto fragment_length=max(static_cast<uintSeqLen>(1),static_cast<uintSeqLen>(stats.FragmentDistribution().InsertLengths().from()))-bias_mod.fragment_length_extension_; fragment_length < fragment_length_to; ++fragment_length){
					if( ref.GC(unit.ref_seq_id_, cur_end_position) ){
						++gc;
					}

					ref.UpdateReverseSurrounding( surrounding_end, ref.ReferenceSequence(unit.ref_seq_id_), cur_end_position++ );

					do{ //while( bias_mod.fragment_length_extension_ )
						UpdateBiasModForCurrentFragmentLength( bias_mod, unit.ref_seq_id_, ref, cur_start_position, cur_end_position, surrounding_end );

						if(stats.FragmentDistribution().InsertLengths().at(fragment_length+bias_mod.fragment_length_extension_)){
							for(uintAlleleId allele=0; allele < ref.NumAlleles(); ++allele){
								if( !ref.VariantsLoaded() || (!AlleleSkipped(bias_mod, allele, ref.Variants(unit.ref_seq_id_), cur_start_position) && cur_end_position+bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_ <= ref.SequenceLength(unit.ref_seq_id_)) ){
									if(ref.VariantsLoaded()){
										// Add/subtract the gc from the additional.left out reference bases due to variant indels
										if(bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_ > 0){
											tmp_gc = ref.GCContentAbsolut( unit.ref_seq_id_, cur_end_position, cur_end_position+bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_ );
										}
										else{
											tmp_gc = - static_cast<intSeqShift>(ref.GCContentAbsolut( unit.ref_seq_id_, cur_end_position+bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_, cur_end_position ));
										}
									}

									for(bool strand : { false, true }){
										// Determine how many read pairs are generated for this strand and allele at this position with this fragment_length
										steady_clock::time_point t1a = steady_clock::now();
										if(ref.VariantsLoaded()){
											fragment_counts = stats.FragmentDistribution().GetFragmentCounts(ref, bias_normalization_, unit.ref_seq_id_, fragment_length+bias_mod.fragment_length_extension_, Percent(gc+bias_mod.gc_mod_.at(allele)+tmp_gc,fragment_length+bias_mod.fragment_length_extension_), bias_mod.mod_surrounding_start_.at(allele), bias_mod.mod_surrounding_end_.at(allele), rdist.ZeroToOne(rgen), non_zero_threshold_);
										}
										else{
											fragment_counts = stats.FragmentDistribution().GetFragmentCounts(ref, bias_normalization_, unit.ref_seq_id_, fragment_length, Percent(gc,fragment_length), surrounding_start, surrounding_end, rdist.ZeroToOne(rgen), non_zero_threshold_);
										}
										steady_clock::time_point t2a = steady_clock::now();
										duration<double> time_span = duration_cast<duration<double>>(t2a - t1a);
										time.get_counts_ += time_span.count();

										if( fragment_counts ){
											t1a = steady_clock::now();

											GetOrgSeq(sim_reads, strand, allele, fragment_length, cur_start_position, cur_end_position, unit.ref_seq_id_, ref, stats, bias_mod);

											if(ref.VariantsLoaded()){
												if(!CreateReads( ref, stats, estimates, rdist, sim_reads, rgen, fragment_counts, strand, allele, unit.ref_seq_id_, fragment_length+bias_mod.fragment_length_extension_, read_number, &block, cur_start_position, cur_end_position+bias_mod.end_pos_shift_.at(allele)+bias_mod.fragment_length_extension_, bias_mod.StartVariant(), bias_mod.EndVariant(ref.Variants(unit.ref_seq_id_), allele, cur_end_position) )){
													return false;
												}
											}
											else{
												if(!CreateReads( ref, stats, estimates, rdist, sim_reads, rgen, fragment_counts, strand, allele, unit.ref_seq_id_, fragment_length, read_number, &block, cur_start_position, cur_end_position )){
													return false;
												}
											}

											t2a = steady_clock::now();
											time_span = duration_cast<duration<double>>(t2a - t1a);
											time.create_reads_ += time_span.count();
										}
									}
								}
							}
						}

						CheckForFragmentLengthExtension( bias_mod, fragment_length, fragment_length_to, ref, stats );
					} while(bias_mod.fragment_length_extension_);
				}
			}

			CheckForInsertedBasesToStartFrom(bias_mod, unit.ref_seq_id_, cur_start_position, ref);
		} while(bias_mod.start_variant_pos_);
	}

	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	time.handle_block_ += time_span.count();

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
	Benchmark time;
	GeneralRandomDistributions rdist(stats);

	steady_clock::time_point t1 = steady_clock::now();

	while( !self.simulation_error_ && self.GetNextBlock(ref, stats, estimates, block, unit, time) ){
		// As long as there are new blocks to simulate from, simulate
		rdist.Reset();
		self.SimulateFromGivenBlock(*block, *unit, ref, stats, estimates, rdist, rgen, time);
		block->finished_ = true;
	}

	// Adapters are positioned after all blocks are created to avoid a race condition with the first block to guarantee reproducible seeds to each block and the adapters
	if(!self.simulation_error_ && !self.adapter_only_simulated_.test_and_set()){
		self.SimulateAdapterOnlyPairs(ref, stats, estimates, rdist, rgen);
	}

	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	time.total_ += time_span.count();

	self.time_mutex_.lock();
	self.time_.total_ += time.total_;
	self.time_.create_block_ += time.create_block_;
	self.time_.handle_block_ += time.handle_block_;
	self.time_.get_counts_ += time.get_counts_;
	self.time_.create_reads_ += time.create_reads_;
	self.time_mutex_.unlock();
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

void Simulator::Simulate(
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
		const string &var_file
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
			if(0 != num_read_pairs || 0.0 != coverage){
				if(0 != num_read_pairs){
					total_pairs_ = num_read_pairs;
				}
				else{ //0.0 != coverage
					total_pairs_ = CoverageToNumberPairs(coverage, total_ref_size, average_read_length, adapter_part);
				}
				// Keep the percentage of adapter only pairs
				num_adapter_only_pairs_ = round( static_cast<double>(total_pairs_)*stats.FragmentDistribution().InsertLengths()[0]/(stats.TotalNumberReads()/2) );
			}
			else{
				// Adapter only pairs(InsertLength == 0) are handled separately
				num_adapter_only_pairs_ = stats.FragmentDistribution().InsertLengths()[0];
				 // We only need to generate a number of pairs equal to half of the reads
				total_pairs_ = stats.TotalNumberReads()/2;
			}

			printInfo << "Aiming for " << total_pairs_ << " read pairs, which corresponds to an approximated read depth of " << setprecision(4) << NumberPairsToCoverage(total_pairs_, total_ref_size, average_read_length, adapter_part) << setprecision(6) << "x" << std::endl;
			// Adapter only pairs(InsertLength == 0) are handled separately
			total_pairs_ -= num_adapter_only_pairs_;

			printInfo << "Preparing for simulation" << std::endl;
			simulation_error_ = false;

			// Prepare vcf file handle if needed
			if(!var_file.empty()){
				if(ref.PrepareVariantFile(var_file)){
					if(!ref.ReadFirstVariants()){
						simulation_error_ = true;
					}
				}
				else{
					simulation_error_ = true;
				}
			}

			if(!simulation_error_){
				double max_bias;
				bias_normalization_ = stats.FragmentDistribution().CalculateBiasNormalization(max_bias, ref, num_threads, total_pairs_) / ref.NumAlleles();
				non_zero_threshold_ = stats.FragmentDistribution().CalculateNonZeroThreshold(bias_normalization_, max_bias);

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

			if(!simulation_error_){
				adapter_only_simulated_.clear();

				time_.total_ = 0.0;
				time_.create_block_ = 0.0;
				time_.handle_block_ = 0.0;
				time_.get_counts_ = 0.0;
				time_.create_reads_ = 0.0;

				// Create a buffer of blocks, so that the blocks into which long fragments reach already exist
				if( CreateBlock(ref, stats, estimates, time_) ){
					for( auto n_blocks = stats.FragmentDistribution().InsertLengths().to()/kBlockSize; n_blocks--; ){
						CreateBlock(ref, stats, estimates, time_);
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

			printDebug << "Time total: " << time_.total_ << std::endl;
			printDebug << "Time create block_: " << time_.create_block_ << std::endl;
			printDebug << "Time handle block_: " << time_.handle_block_ << std::endl;
			printDebug << "Time get counts_: " << time_.get_counts_ << std::endl;
			printDebug << "Time create reads_: " << time_.create_reads_ << std::endl;

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
	}

	if(simulation_error_){
		// Remove the simulated files so it is clear that we encountered an error and do not accidentally continue with the old file
		DeleteFile( destination_file_first );
		DeleteFile( destination_file_second );
	}
}
