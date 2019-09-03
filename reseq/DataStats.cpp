#include "DataStats.h"
using reseq::SeqQualityStats;
using reseq::DataStats;
using reseq::Vect;

#include <algorithm>
using std::max;
using std::min;
//include <array>
using std::array;
#include <exception>
using std::exception;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <limits>
using std::numeric_limits;
#include <mutex>
using std::lock_guard;
using std::mutex;
using std::unique_lock;
#include <stdexcept>
using std::runtime_error;
//include <string>
using std::string;
//include <unordered_set>
using std::unordered_set;
//include <utility>
using std::pair;
#include <thread>
using std::thread;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/bam_io.h>
using seqan::atEnd;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamHeader;
using seqan::CharString;
using seqan::CigarElement;
using seqan::contigNames;
using seqan::Dna;
using seqan::Dna5;
using seqan::Dna5String;
using seqan::Exception;
using seqan::FormattedFileContext;
using seqan::FunctorComplement;
using seqan::IupacString;
using seqan::ModComplementDna5;
using seqan::ModifiedString;
using seqan::ModReverse;
using seqan::readRecord;
using seqan::readHeader;
using seqan::Size;
using seqan::String;

#include "CMakeConfig.h"
//include "utilities.h"
using reseq::utilities::CreateDir;
using reseq::utilities::Divide;
using reseq::utilities::MeanWithRoundingToFirst;
using reseq::utilities::Percent;
using reseq::utilities::SafePercent;
using reseq::utilities::SetDominantLastX;
using reseq::utilities::SetToMin;
using reseq::utilities::SetToMax;

uint64_t DataStats::RecordHasher::operator() (CoverageStats::FullRecord * const &record) const{
	// FNV-1a hash: https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
	const auto cdata = reinterpret_cast<unsigned char *>( toCString(record->record_.qName) );
	uint64_t hash = 14695981039346656037ull;
	for( auto i=length(record->record_.qName); i--; ){
		const auto next = std::size_t {cdata[i]};
		hash = (hash ^ next) * 1099511628211ull;
	}

	return hash;
}


bool  DataStats::RecordEqual::operator() (CoverageStats::FullRecord * const &lhs, CoverageStats::FullRecord * const &rhs) const{
	return lhs->record_.qName == rhs->record_.qName;
}

inline bool DataStats::PotentiallyValid( const BamAlignmentRecord &record_first ) const{
	return !hasFlagNextUnmapped(record_first) && !hasFlagUnmapped(record_first) && record_first.rNextId == record_first.rID && maximum_insert_length_ >= record_first.pNext - record_first.beginPos && reference_->SequenceLength(record_first.rID) > maximum_insert_length_+2*reference_->MinDistToRefSeqEnds();
}

bool DataStats::IsSecondRead( CoverageStats::FullRecord *record, CoverageStats::FullRecord *&record_first, CoverageStats::CoverageBlock *&block ){
	auto insert_info = first_read_records_.insert( record );
	if( insert_info.second ){
		// Insertion worked, so it is the first read
		if( PotentiallyValid(record->record_) ){
			coverage_.AddFragment( record->record_.rID, record->record_.beginPos, block );
		}

		return false;
	}
	else{
		// Insertion failed, so it is the second read
		record_first = *insert_info.first;
		first_read_records_.erase(insert_info.first);

		return true;
	}
}

bool DataStats::CheckForAdapters(const seqan::BamAlignmentRecord &record_first, const seqan::BamAlignmentRecord &record_second){
	uint32_t adapter_position_first, adapter_position_second;
	bool adapter_detected = adapters_.Detect(adapter_position_first, adapter_position_second, record_first, record_second, *reference_);

	if(adapter_detected){
		auto insert_length = MeanWithRoundingToFirst(adapter_position_first,adapter_position_second);
		fragment_distribution_.AddInsertLengths(insert_length);
		++tmp_read_lengths_by_fragment_length_.at(0).at(insert_length).at(length(record_first.seq));
		++tmp_read_lengths_by_fragment_length_.at(1).at(insert_length).at(length(record_second.seq));

		return true;
	}
	else{
		return false;
	}
}

void DataStats::EvalBaseLevelStats( CoverageStats::FullRecord *full_record, uint16_t template_segment, uint16_t tile_id, uint8_t &paired_seq_qual ){
	const BamAlignmentRecord &record(full_record->record_);

	Size<CharString>::Type pos_reversed(length(record.seq));
	FunctorComplement<Dna5> complementor;

	uint16_t strand = hasFlagRC(record) && hasFlagFirst(record) || !hasFlagRC(record) && hasFlagLast(record);

	SeqQualityStats<uint64_t> seq_qual_stats;
	seq_qual_stats[maximum_quality_]; // Resize to maximum quality
	for( auto qual : record.qual ){
		++seq_qual_stats.at(qual - phred_quality_offset_);
	}
	seq_qual_stats.Calculate();
	seq_qual_stats.CalculateProbabilityMean();
	full_record->sequence_quality_ = seq_qual_stats.mean_;

	Dna5 base;
	uint16_t qual, last_qual(1);
	array<uint64_t, 5> read_bases = {0,0,0,0,0};
	uint32_t homoquality_length = 0;
	uint32_t homopolymer_length = 0;
	Dna5 homopolymer_nucleotide = 0;

	for( Size<CharString>::Type pos = 0; pos < length(record.seq); ++pos){
		// In case of reversed sequences base and qual are read in reverse from bamfile so pos represents the position in the fq file
		if( hasFlagRC(record) ){
			base = complementor(record.seq[--pos_reversed]);
			qual = record.qual[pos_reversed]-phred_quality_offset_;
		}
		else{
			base = record.seq[pos];
			qual = record.qual[pos]-phred_quality_offset_;
		}

		qualities_.AddRawBase(template_segment, base, tile_id, strand, qual, seq_qual_stats.mean_, last_qual, pos);

		++read_bases.at(base);
		++tmp_sequence_content_.at(template_segment).at(base).at(pos);

		if( last_qual == qual ){
			++homoquality_length;
		}
		else{
			if( pos ){
				qualities_.AddRawHomoqualimer(last_qual, homoquality_length);
			}
			homoquality_length = 1;
		}

		if( homopolymer_nucleotide == base ){
			++homopolymer_length;
		}
		else{
			++tmp_homopolymer_distribution_.at(homopolymer_nucleotide).at(homopolymer_length);
			homopolymer_nucleotide = base;
			homopolymer_length = 1;
		}

		// Quality based on preceding quality
		last_qual = qual;
	}

	qualities_.AddRawHomoqualimer(last_qual, homoquality_length);
	++tmp_homopolymer_distribution_.at(homopolymer_nucleotide).at(homopolymer_length); // Add the homopolymer at read end

	// Read level summaries of base level stats
	qualities_.AddRawRead(paired_seq_qual, seq_qual_stats, template_segment, tile_id, read_bases, length(record.seq));

	++tmp_gc_read_content_.at(template_segment).at(Percent( read_bases.at(1)+read_bases.at(2), length(record.seq) ));
	++tmp_n_content_.at(template_segment).at(Percent( read_bases.at(4), length(record.seq) ));
}

bool DataStats::EvalReferenceStatistics(
		CoverageStats::FullRecord *record,
		uint16_t template_segment,
		CoverageStats::CoverageBlock *coverage_block,
		ThreadData& thread_data,
		uint32_t paired_read_end ){
	// Get gc on reference
	array<uint64_t, 5> seq_content_reference = {0,0,0,0,0};
	auto cur_var = coverage_block->first_variant_id_;
	CoverageStats::PrepareVariantPositionCheck( cur_var, record->record_.rID, record->from_ref_pos_, *reference_, false );
	for( auto ref_pos = record->from_ref_pos_; ref_pos < record->to_ref_pos_; ++ref_pos ){
		if( CoverageStats::IsVariantPosition( cur_var, record->record_.rID, ref_pos, *reference_, false )){
			++seq_content_reference.at(4);
		}
		else{
			++seq_content_reference.at( reference_->ReferenceSequence(record->record_.rID)[ref_pos] );
		}
	}
	auto gc_percent = SafePercent( seq_content_reference.at(1)+seq_content_reference.at(2), record->to_ref_pos_-record->from_ref_pos_-seq_content_reference.at(4) );
	record->reference_gc_ = gc_percent;

	uint32_t coverage_pos;
	if( hasFlagRC(record->record_) ){
		coverage_pos = coverage_.GetStartPos(record->to_ref_pos_ - 1, coverage_block);
	}
	else{
		coverage_pos = coverage_.GetStartPos(record->from_ref_pos_, coverage_block);
	}

	cur_var = coverage_block->first_variant_id_;
	if( reference_->VariantPositionsLoaded() && hasFlagRC(record->record_) ){
		CoverageStats::PrepareVariantPositionCheck( cur_var, record->record_.rID, record->to_ref_pos_-1, *reference_, true );
	}

	ModifiedString< ModifiedString<const IupacString, ModComplementDna5>, ModReverse> reversed_seq(record->record_.seq);
	ModifiedString< const CharString, ModReverse> reversed_qual(record->record_.qual);
	FunctorComplement<Dna5> complementor;
	ModifiedString<const String<CigarElement<>>, ModReverse> rev_cigar(record->record_.cigar);

	uint32_t read_pos(0), ref_pos(0), full_ref_pos;
	Dna ref_base;
	Dna5 base;
	uint8_t last_base(5);
	uint8_t qual;

	array<uint64_t, 5> seq_content_mapped = {0,0,0,0,0};
	uint32_t num_errors(0), indel_pos(0);
	uint16_t indel_type(0);

	for( const auto &cigar_element : (hasFlagRC(record->record_)?rev_cigar:record->record_.cigar) ){
		switch( cigar_element.operation ){
		case 'M':
		case '=':
		case 'X':
			for(auto i=cigar_element.count; i--; ){
				if(ref_pos < record->to_ref_pos_-record->from_ref_pos_){
					// We are currently not comparing the adapter to the reference
					if( hasFlagRC(record->record_) ){
						base = reversed_seq[read_pos];
						qual = reversed_qual[read_pos]-phred_quality_offset_;
						full_ref_pos = record->to_ref_pos_-1-ref_pos;
						ref_base = complementor(reference_->ReferenceSequence(record->record_.rID)[full_ref_pos]);
					}
					else{
						base = record->record_.seq[read_pos];
						qual = record->record_.qual[read_pos]-phred_quality_offset_;
						full_ref_pos = record->from_ref_pos_+ref_pos;
						ref_base = reference_->ReferenceSequence(record->record_.rID)[full_ref_pos];
					}

					if(static_cast<uint16_t>(ref_base) < 4 && !CoverageStats::IsVariantPosition( cur_var, record->record_.rID, full_ref_pos, *reference_, hasFlagRC(record->record_) )){ // Don't do this stuff when we have an N
						++tmp_sequence_content_reference_.at(template_segment).at(hasFlagRC(record->record_)).at(ref_base).at(read_pos);

						errors_.AddBasePlotting(template_segment, ref_base, base, qual, last_base);
						errors_.AddInDel( indel_type, last_base, ErrorStats::kNoInDel, indel_pos, read_pos, gc_percent);

						if(ref_base != base){
							++num_errors;
						}

						++seq_content_mapped.at(base);
					}
					else{
						++seq_content_mapped.at(4);
					}

					if( hasFlagRC(record->record_) ){
						coverage_.AddReverse(coverage_pos, coverage_block, base);
					}
					else{
						coverage_.AddForward(coverage_pos, coverage_block, base);
					}

					last_base = base;

					++ref_pos;
					++read_pos;
					if( hasFlagRC(record->record_) ){
						coverage_.DecrementPos(coverage_pos, coverage_block);
					}
					else{
						coverage_.IncrementPos(coverage_pos, coverage_block);
					}

					indel_type = 0;
					indel_pos = 0;
				}
			}
			break;
		case 'N':
		case 'D':
			if(ref_pos < record->to_ref_pos_-record->from_ref_pos_){
				// We are currently not comparing the adapter to the reference
				for(indel_pos=0; indel_pos<cigar_element.count; ++indel_pos ){
					if( hasFlagRC(record->record_) ){
						full_ref_pos = record->to_ref_pos_-1-ref_pos;
						ref_base = complementor(reference_->ReferenceSequence(record->record_.rID)[full_ref_pos]);
					}
					else{
						full_ref_pos = record->from_ref_pos_+ref_pos;
						ref_base = reference_->ReferenceSequence(record->record_.rID)[full_ref_pos];
					}

					if(static_cast<uint16_t>(ref_base) < 4 && !CoverageStats::IsVariantPosition( cur_var, record->record_.rID, full_ref_pos, *reference_, hasFlagRC(record->record_) )){ // Don't do this stuff when we have an N
						errors_.AddInDel( indel_type, last_base, ErrorStats::kDeletion, indel_pos, read_pos, gc_percent);
					}

					++ref_pos;
					indel_type = 1;
				}

				num_errors += cigar_element.count;

				if( hasFlagRC(record->record_) ){
					coverage_.SubtractPos(coverage_pos, coverage_block, cigar_element.count);
				}
				else{
					coverage_.AddPos(coverage_pos, coverage_block, cigar_element.count);
				}
			}
			break;
		case 'I':
			if(ref_pos < record->to_ref_pos_-record->from_ref_pos_){
				if( hasFlagRC(record->record_) ){
					full_ref_pos = record->to_ref_pos_-1-ref_pos;
					ref_base = complementor(reference_->ReferenceSequence(record->record_.rID)[full_ref_pos]);
				}
				else{
					full_ref_pos = record->from_ref_pos_+ref_pos;
					ref_base = reference_->ReferenceSequence(record->record_.rID)[full_ref_pos];
				}

				if(static_cast<uint16_t>(ref_base) >= 4 || CoverageStats::IsVariantPosition( cur_var, record->record_.rID, full_ref_pos, *reference_, hasFlagRC(record->record_) )){
					read_pos += cigar_element.count;
					indel_type = 0;
				}
				else{
					// Don't do this stuff when we have an N
					for(indel_pos=0; indel_pos<cigar_element.count; ++indel_pos ){
						if( hasFlagRC(record->record_) ){
							base = reversed_seq[read_pos];
						}
						else{
							base = record->record_.seq[read_pos];
						}

						errors_.AddInDel( indel_type, last_base, static_cast<ErrorStats::InDelDef>(static_cast<uint16_t>(base)+2), indel_pos, read_pos, gc_percent);
						++read_pos;
						indel_type = 0;
					}
					num_errors += cigar_element.count;
				}
				seq_content_mapped.at(4) += cigar_element.count; // Ignore inserted bases for gc_percent
			}
			break;
		}
	}

	++tmp_gc_read_content_reference_.at(template_segment).at(gc_percent);
	++tmp_gc_read_content_mapped_.at(template_segment).at(SafePercent( seq_content_mapped.at(1) + seq_content_mapped.at(2), length(record->record_.seq)-seq_content_mapped.at(4) ));

	errors_.AddRead(template_segment, num_errors);

	return true;
}

bool DataStats::EvalRecord( pair<CoverageStats::FullRecord *, CoverageStats::FullRecord *> &record, ThreadData &thread_data ){
	// Start handling pair by determining tile(and counting it) and template segment
	uint16_t tile_id;
	if( !tiles_.GetTileId(tile_id, record.first->record_.qName) ){
		return false;
	}

	uint16_t template_segment(0);
	if( hasFlagLast(record.second->record_) ){
		template_segment = 1;
	}

	// Base level stats
	unsigned char paired_seq_qual(0);
	EvalBaseLevelStats( record.first, (template_segment+1)%2, tile_id, paired_seq_qual );
	EvalBaseLevelStats( record.second, template_segment, tile_id, paired_seq_qual );

	// Mapping and Reference stats
	if( hasFlagUnmapped(record.first->record_) || hasFlagNextUnmapped(record.first->record_) ){
		if( !hasFlagUnmapped(record.first->record_) ){
			++tmp_single_read_mapping_quality_.at(record.first->record_.mapQ);
		}
		else if( !hasFlagUnmapped(record.second->record_) ){
			++tmp_single_read_mapping_quality_.at(record.second->record_.mapQ);
		}

		if(CheckForAdapters(record.first->record_, record.second->record_)){
			reads_in_unmapped_pairs_with_adapters_ += 2;
		}
		else{
			reads_in_unmapped_pairs_without_adapters_ += 2;
		}
	}
	else{
		if( !hasFlagAllProper(record.first->record_) ){
			++tmp_improper_pair_mapping_quality_.at(record.first->record_.mapQ);
			++tmp_improper_pair_mapping_quality_.at(record.second->record_.mapQ);
		}
		else{
			++tmp_proper_pair_mapping_quality_.at(record.first->record_.mapQ);
			++tmp_proper_pair_mapping_quality_.at(record.second->record_.mapQ);
		}

		if( !QualitySufficient( record.first->record_ ) || !QualitySufficient( record.second->record_ ) ){
			CheckForAdapters(record.first->record_, record.second->record_); // We don't trust the mapping, so look for adapters like in unmapped case

			reads_with_low_quality_ += 2;
		}
		else{
			if( record.first->record_.rID == record.second->record_.rID ){
				if( reference_->SequenceLength(record.first->record_.rID) > maximum_insert_length_+2*reference_->MinDistToRefSeqEnds() ){
					auto read_length_on_reference_first = GetReadLengthOnReference(record.first->record_);
					auto end_pos_first = record.first->record_.beginPos + read_length_on_reference_first;

					auto read_length_on_reference_second = GetReadLengthOnReference(record.second->record_);
					auto end_pos_second = record.second->record_.beginPos + read_length_on_reference_second;

					if(min_dist_to_ref_seq_ends_ <= min(record.first->record_.beginPos, record.second->record_.beginPos) && reference_->SequenceLength(record.first->record_.rID)-max(end_pos_first, end_pos_second) >= min_dist_to_ref_seq_ends_ ){
						bool adapter_detected(false);
						if(hasFlagRC(record.first->record_) && !hasFlagRC(record.second->record_) && end_pos_first > record.second->record_.beginPos && end_pos_first - record.second->record_.beginPos > read_length_on_reference_first/2 ){
							uint32_t adapter_position_first, adapter_position_second;
							// Reads in proper direction for adapters (reverse-forward) and overlap is large enough (minimum half of the read length)
							adapters_.Detect(adapter_position_first, adapter_position_second, record.first->record_, record.second->record_, *reference_, true); // Call it just to insert the adapter stats, position is better known from mapping
							adapter_detected = true;
						}

						if(InProperDirection( record.first->record_, end_pos_second ) || adapter_detected){
							reads_used_ += 2;

							CoverageStats::CoverageBlock *coverage_block;
							// Pair stats
							uint32_t insert_length;
							if(adapter_detected){
								insert_length = end_pos_first - record.second->record_.beginPos;

								fragment_distribution_.AddGCContent(reference_->GCContent(record.second->record_.rID, record.second->record_.beginPos, end_pos_first));
								fragment_distribution_.FillInOutskirtContent(*reference_, record.second->record_, end_pos_first); // Use second record here as it is forward, so if it is sequenced second, fragment is reversed
								coverage_block = coverage_.FindBlock( record.second->record_.rID, record.second->record_.beginPos );

								record.first->from_ref_pos_ = record.second->record_.beginPos;
								record.second->from_ref_pos_ = record.second->record_.beginPos;
								record.first->to_ref_pos_ = end_pos_first;
								record.second->to_ref_pos_ = end_pos_first;

								fragment_distribution_.AddFragmentSite(record.first->record_.rID, insert_length, record.first->from_ref_pos_, template_segment, *reference_); // If second record(must be forward to be accepted due to adapter) is sequenced second, fragment is reversed
							}
							else{
								insert_length = end_pos_second - record.first->record_.beginPos;

								fragment_distribution_.AddGCContent(reference_->GCContent(record.first->record_.rID, record.first->record_.beginPos, end_pos_second));
								fragment_distribution_.FillInOutskirtContent(*reference_, record.first->record_, end_pos_second);
								coverage_block = coverage_.FindBlock( record.first->record_.rID, record.first->record_.beginPos );

								record.first->from_ref_pos_ = record.first->record_.beginPos;
								record.second->from_ref_pos_ = record.second->record_.beginPos;
								record.first->to_ref_pos_ = end_pos_first;
								record.second->to_ref_pos_ = end_pos_second;

								fragment_distribution_.AddFragmentSite(record.first->record_.rID, insert_length, record.first->from_ref_pos_, (template_segment+1)%2, *reference_); // If first record(must be forward to be accepted without adapters) is sequenced second, fragment is reversed
							}

							fragment_distribution_.AddAbundance(record.first->record_.rID);
							fragment_distribution_.AddInsertLengths(insert_length); // We do not just use the template length from field 9 of the bam file, because bwa stores there the right-most position of the reverse read - the left-most position of the forward read, which is not the fragment length in the case of the reverse read being positioned before the forward read. Bowtie2 gives the correct length here, but I don't want to sacrifice compatibility so easily

							++tmp_read_lengths_by_fragment_length_.at(0).at(insert_length).at(length(record.first->record_.seq));
							++tmp_read_lengths_by_fragment_length_.at(1).at(insert_length).at(length(record.second->record_.seq));

							record.first->tile_id_ = tile_id;
							record.second->tile_id_ = tile_id;
							record.first->fragment_length_ = insert_length;
							record.second->fragment_length_ = insert_length;

							if( !EvalReferenceStatistics( record.first, (template_segment+1)%2, coverage_block, thread_data, 0 ) ){
								return false;
							}
							if( !EvalReferenceStatistics( record.second, template_segment, coverage_block, thread_data, end_pos_first ) ){
								return false;
							}
						}
					}
					else{
						reads_to_close_to_ref_seq_ends_ += 2;
					}
				}
				else{
					reads_on_too_short_fragments_ += 2;
				}
			}
		}
	}

	return true;
}

bool DataStats::SignsOfPairsWithNamesNotIdentical(){
	// All reads have been processed, so first_read_records_ must be empty
	if( first_read_records_.size() ){
		printErr << "Read names of reads from the same pair are not identical. Cannot continue.\nExample position: " << (*first_read_records_.begin())->record_.rID << ":" << (*first_read_records_.begin())->record_.beginPos << "\nExample read name: " << (*first_read_records_.begin())->record_.qName << '\n';
		return true;
	}

	return false;
}

void DataStats::PrepareReadIn(uint8_t size_mapping_quality, uint32_t size_indel){
	uint32_t size_pos = max(read_lengths_.at(0).to(), read_lengths_.at(1).to());

	uint64_t num_reads(0);
	uint64_t num_bases(0);
	for( int template_segment = 2; template_segment--; ){
		for(auto len=read_lengths_.at(template_segment).from(); len < read_lengths_.at(template_segment).to(); ++len ){
			num_reads += read_lengths_.at(template_segment).at(len);
			num_bases += read_lengths_.at(template_segment).at(len)*len;
		}
	}

	coverage_.Prepare( Divide(num_bases,reference_->TotalSize()), Divide(num_bases,num_reads) );
	duplicates_.PrepareTmpDuplicationVector(maximum_insert_length_);
	errors_.Prepare(tiles_.NumTiles(), maximum_quality_+1, size_pos, size_indel);
	fragment_distribution_.Prepare(*reference_, maximum_insert_length_, reads_per_ref_seq_bin_);
	qualities_.Prepare(tiles_.NumTiles(), maximum_quality_+1, size_pos, maximum_insert_length_);

	// Prepare vector in this class
	tmp_proper_pair_mapping_quality_.resize(size_mapping_quality);
	tmp_improper_pair_mapping_quality_.resize(size_mapping_quality);
	tmp_single_read_mapping_quality_.resize(size_mapping_quality);

	for( auto template_segment=2; template_segment--; ){
		SetDimensions( tmp_read_lengths_by_fragment_length_.at(template_segment), maximum_insert_length_+1, size_pos );

		tmp_gc_read_content_.at(template_segment).resize(101);
		tmp_gc_read_content_reference_.at(template_segment).resize(101);
		tmp_gc_read_content_mapped_.at(template_segment).resize(101);

		tmp_n_content_.at(template_segment).resize(101);
		for( auto base=5; base--; ){
			tmp_sequence_content_.at(template_segment).at(base).resize(size_pos);
		}
		for( auto strand=2; strand--; ){
			for( auto base=4; base--; ){
				tmp_sequence_content_reference_.at(template_segment).at(strand).at(base).resize(maximum_read_length_on_reference_+1);
			}
		}
	}

	for( auto base=5; base--; ){
		tmp_homopolymer_distribution_.at(base).resize(size_pos);
	}
}

bool DataStats::FinishReadIn(){
	// Collect ambigous adapters into a single one
	adapters_.Finalize();

	if( !coverage_.Finalize(*reference_, qualities_, errors_, phred_quality_offset_, maximum_insert_length_+2*reference_->MinDistToRefSeqEnds()) ){
		return false;
	}
	if( !errors_.Finalize() ){
		return false;
	}
	fragment_distribution_.Finalize();
	qualities_.Finalize( SumVect(read_lengths_.at(0)) );

	// Copy vectors to final ones
	proper_pair_mapping_quality_.Acquire(tmp_proper_pair_mapping_quality_);
	improper_pair_mapping_quality_.Acquire(tmp_improper_pair_mapping_quality_);
	single_read_mapping_quality_.Acquire(tmp_single_read_mapping_quality_);

	for( auto template_segment=2; template_segment--; ){
		read_lengths_by_fragment_length_.at(template_segment).Acquire(tmp_read_lengths_by_fragment_length_.at(template_segment));

		gc_read_content_.at(template_segment).Acquire(tmp_gc_read_content_.at(template_segment));
		gc_read_content_reference_.at(template_segment).Acquire(tmp_gc_read_content_reference_.at(template_segment));
		gc_read_content_mapped_.at(template_segment).Acquire(tmp_gc_read_content_mapped_.at(template_segment));

		n_content_.at(template_segment).Acquire(tmp_n_content_.at(template_segment));
		for( auto base=5; base--; ){
			sequence_content_.at(template_segment).at(base).Acquire(tmp_sequence_content_.at(template_segment).at(base));
		}
		for( auto strand=2; strand--; ){
			for( auto base=4; base--; ){
				sequence_content_reference_.at(template_segment).at(strand).at(base).Acquire(tmp_sequence_content_reference_.at(template_segment).at(strand).at(base));
			}
		}
	}

	tmp_homopolymer_distribution_.at(0).at(0) = 0; // Remove homopolymers introduced by initialization values
	for( auto base=5; base--; ){
		homopolymer_distribution_.at(base).Acquire(tmp_homopolymer_distribution_.at(base));
	}

	return true;
}

void DataStats::Shrink(){
	adapters_.Shrink();
	coverage_.Shrink();
	errors_.Shrink();
	fragment_distribution_.Shrink();
	qualities_.Shrink();
	tiles_.Shrink();
	
	for( uint16_t template_segment=2; template_segment--; ){
		ShrinkVect(read_lengths_by_fragment_length_.at(template_segment));

		gc_read_content_.at(template_segment).Shrink();
		gc_read_content_reference_.at(template_segment).Shrink();
		gc_read_content_mapped_.at(template_segment).Shrink();

		n_content_.at(template_segment).Shrink();
		for( int called_base=5; called_base--; ){
			sequence_content_.at(template_segment).at(called_base).Shrink();
			homopolymer_distribution_.at(called_base).Shrink();
		}
		for( int ref_base=4; ref_base--; ){
			for( int strand=2; strand--; ){
				sequence_content_reference_.at(template_segment).at(strand).at(ref_base).Shrink();
			}
		}
	}

	proper_pair_mapping_quality_.Shrink();
	improper_pair_mapping_quality_.Shrink();
	single_read_mapping_quality_.Shrink();
}

// Deactivation of bias calculation only for speeding up of tests
bool DataStats::Calculate(uint16_t num_threads){
	// Calculate the duplicates from observation of fragments at same position and reference characteristics
	if(!fragment_distribution_.FinalizeBiasCalculation(*reference_, num_threads, duplicates_)){
		return false;
	}

	return true;
}

void DataStats::PrepareGeneral(){
	// Set total_number_reads_ to the number of accepted reads as it is the number of read records after read in
	total_number_reads_ = SumVect(read_lengths_.at(0)) + SumVect(read_lengths_.at(1));

	adapters_.SumCounts();
}

bool DataStats::OrderOfBamFileCorrect(const seqan::BamAlignmentRecord &record, pair< uint32_t, uint32_t > last_record_pos){
	if( !hasFlagUnmapped(record) ){
		if( record.rID < last_record_pos.first || record.rID == last_record_pos.first && record.beginPos < last_record_pos.second ){
			printErr << "Bamfile is not sorted by position. Detected at position " << record.beginPos << " with read name: " << record.qName << '\n';
			return false;
		}
		last_record_pos = {record.rID, record.beginPos};
	}

	return true;
}

bool DataStats::PreRun(BamFileIn &bam, const char *bam_file, BamHeader &header, uint8_t &size_mapping_quality, uint32_t &size_indel, uint64_t max_ref_seq_bin_size){
	printInfo << "Starting PreRun\n";

	bool error = false;
	BamAlignmentRecord record;
	pair< uint32_t, uint32_t > last_record_pos{0,0};

	uint32_t read_length_on_reference;
	uint8_t min_mapq(255), max_mapq(0);
	size_indel = 0;

	auto num_ref_bins = fragment_distribution_.CreateRefBins(*reference_, max_ref_seq_bin_size);
	reads_per_ref_seq_bin_.resize(num_ref_bins);

	do{
		try{
			++total_number_reads_; // Counter total_number_reads_ for read in, later in Calculate correct it with the number of accepted reads from read_length_

			readRecord(record, bam);

			if( IsValidRecord(record) ){
				if( OrderOfBamFileCorrect(record, last_record_pos) ){
					if(hasFlagFirst(record)){
						tiles_.EnterTile(record.qName); // Names of records in a pair are identical, so tile must only be identified once

						++read_lengths_.at(0)[ length(record.qual) ];
					}
					else{
						++read_lengths_.at(1)[ length(record.qual) ];
					}

					// Determine read length range on reference
					if(!hasFlagUnmapped(record)){
						read_length_on_reference = GetReadLengthOnReference(record, size_indel);
						SetToMax(maximum_read_length_on_reference_, read_length_on_reference);
						SetToMin(minimum_read_length_on_reference_, read_length_on_reference);

						auto ref_bin = fragment_distribution_.GetRefSeqBin(record.rID, record.beginPos, *reference_);
						++reads_per_ref_seq_bin_.at(ref_bin);

						// Check if it is close to the border of a bin in which case add it to the other bin as well
						if(record.beginPos+maximum_insert_length_ < reference_->SequenceLength(record.rID)){
							auto ref_bin2 = fragment_distribution_.GetRefSeqBin(record.rID, record.beginPos+maximum_insert_length_, *reference_);
							if(ref_bin2 != ref_bin){
								++reads_per_ref_seq_bin_.at(ref_bin2);
							}
						}
						if(record.beginPos>=maximum_insert_length_){
							auto ref_bin2 = fragment_distribution_.GetRefSeqBin(record.rID, record.beginPos-maximum_insert_length_, *reference_);
							if(ref_bin2 != ref_bin){
								++reads_per_ref_seq_bin_.at(ref_bin2);
							}
						}
					}

					// Determine quality range
					for( auto i=length(record.qual); i--;){
						SetToMax(maximum_quality_, record.qual[i]);
						SetToMin(minimum_quality_, record.qual[i]);
					}

					SetToMax(max_mapq, record.mapQ);
					SetToMin(min_mapq, record.mapQ);

				}
				else{
					error = true;
				}
			}
			else{
				error = true;
			}
		}
		catch(const Exception &e){
			error = true;
			printErr << "Could not read record " << total_number_reads_ << " in " << bam_file << ": " << e.what() << '\n';
		}
	}while( !atEnd(bam) );

	if(error){
		return false;
	}

	// Jump back to beginning of bam file
	seqan::close(bam);
	if( !open(bam, bam_file) ){
		printErr << "Could not open " << bam_file << " for reading.\n";
		return false;
	}

	try{
		readHeader(header, bam);
	}
	catch(const Exception &e){
		printErr << "Could not read header in " << bam_file << ": " << e.what() << '\n';
		return false;
	}

	for( uint16_t template_segment=2; template_segment--; ){
		read_lengths_.at(template_segment).Shrink();
	}

	if( minimum_quality_ < 64 ){
		if( minimum_quality_ < 33 ){
			printErr << "Minimum quality is " << minimum_quality_ << ", which is lower than the 33 from Sanger encoding. Are those values correct?\n";
			return false;
		}
		phred_quality_offset_ = 33; // Sanger format
	}
	else{
		phred_quality_offset_ = 64; // Old Illumina format
	}

	minimum_quality_ -= phred_quality_offset_;
	maximum_quality_ -= phred_quality_offset_;

	printInfo << "Finished PreRun\nTotal number of reads: " << total_number_reads_ << '\n'
			<< "Read length: " << min(read_lengths_.at(0).from(), read_lengths_.at(1).from()) << " - " << (max(read_lengths_.at(0).to(), read_lengths_.at(1).to())-1) << '\n'
			<< "Read length on reference: " << minimum_read_length_on_reference_ << " - " << maximum_read_length_on_reference_ << '\n'
			<< "Quality: " << static_cast<uint16_t>(minimum_quality_) << " - " << static_cast<uint16_t>(maximum_quality_) << " ( " << static_cast<uint16_t>(phred_quality_offset_) << "-based encoding )\n"
			<< "Mapping Quality: " << static_cast<uint16_t>(min_mapq) << " - " << static_cast<uint16_t>(max_mapq) << '\n'
			<< "Maximum InDel Length: " << size_indel << '\n';

	size_mapping_quality = max_mapq+1;
	++size_indel; // Add one to go from index of last value to size

	return true;
}

bool DataStats::ReadRecords( BamFileIn &bam, bool &not_done, ThreadData &thread_data ){
	CoverageStats::FullRecord *record;
	CoverageStats::CoverageBlock *cov_block(NULL);
	lock_guard<mutex> lock(read_mutex_);
	try{
		while( thread_data.rec_store_.size() < batch_size_ && !atEnd(bam) && reading_success_){
			++read_records_;

			record = new CoverageStats::FullRecord;
			readRecord(record->record_, bam);

			if(!hasFlagUnmapped(record->record_) && reference_->SequenceLength(record->record_.rID) > maximum_insert_length_+2*reference_->MinDistToRefSeqEnds()){
				if( !coverage_.EnsureSpace(record->record_.rID, record->record_.beginPos, record->record_.beginPos+maximum_read_length_on_reference_, record, *reference_, maximum_insert_length_+2*reference_->MinDistToRefSeqEnds()) ){
					return false;
				}
			}

			// Work with the pairs and not individual reads, so first combine reads to pairs, by storing first and then processing both reads at the same time
			CoverageStats::FullRecord *record_first;
			if( IsSecondRead(record, record_first, cov_block) ){
				thread_data.rec_store_.emplace_back(record_first, record);
			}

			if( total_number_reads_ > batch_size_ && !(read_records_%(total_number_reads_/20)) ){
				lock_guard<mutex> lock(print_mutex_);
				printInfo << "Read " << static_cast<uint16_t>(Percent(read_records_, total_number_reads_)) << "% of the reads.\n";
			}
		}
		not_done = !atEnd(bam);
	}
	catch(const Exception &e){
		lock_guard<mutex> lock(print_mutex_);
		printErr << "Could not read record " << read_records_ << ": " << e.what() << '\n';
		return false;
	}

	return true;
}

void DataStats::ReadThread( DataStats &self, BamFileIn &bam, uint32_t max_seq_bin_len ){
	CoverageStats::CoverageBlock *cov_block;
	uint64_t processed_fragments(0);
	bool not_done(true);
	ThreadData thread_data( self.maximum_insert_length_, max_seq_bin_len );
	thread_data.rec_store_.reserve(self.batch_size_);
	uint32_t still_needed_reference_sequence(0);
	uint32_t still_needed_position(0);
	while(not_done && self.reading_success_){
		// ReadRecords:
		// Reads in records and bundles the to batches of paired records
		// Registers reads as unprocessed in coverage blocks
		// Adds pointers of all potentially valid reads to each block they potentially overlap (read length on reference yet unknown)
		if( self.ReadRecords(bam, not_done, thread_data) ){
			cov_block = NULL;
			// Process batch
			for( auto &rec : thread_data.rec_store_ ){
				// EvalRecord:
				// Apply multiple filters to check whether fragments are valid
				// Set to_ref_pos for valid reads, that to_ref_pos==0 is used in coverage_.CleanUp to remove the non-valid reads
				// Enters valid fragments into duplication check
				// Calculate all the statistics that are independent of other reads
				// Fills the coverage information in the coverage blocks
				if(self.EvalRecord(rec, thread_data)){
					if( self.PotentiallyValid(rec.first->record_) ){
						// coverage_.RemoveFragment
						// Marks the reads as processed, so the fully processed coverage blocks can be handled after all reads are processed
						self.coverage_.RemoveFragment( rec.first->record_.rID, rec.first->record_.beginPos, cov_block, processed_fragments );
					}

					// Reduce count of in use pointers to records and if 0 remove them
					self.coverage_.DeleteRecord( rec.first, self.qualities_ );
					self.coverage_.DeleteRecord( rec.second, self.qualities_ );
				}
				else{
					self.reading_success_ = false;
				}
			}

			if(cov_block){ // To make sure we don't have a whole batch of non-valid reads as coverage_.CleanUp already handled blocks at the beginning and we need at least one block still existing to keep track of the position
				// coverage_.CleanUp:
				// Once all reads in it have been handled processes the coverage blocks: The coverage itself and all the read statistics that depend on the other reads (systematic error)
				// Removes the processed blocks
				still_needed_reference_sequence = self.coverage_.CleanUp(still_needed_position, *self.reference_, self.qualities_, self.errors_, self.phred_quality_offset_, cov_block, processed_fragments);
			}

			self.coverage_.PreLoadVariants(*self.reference_);

			self.fragment_distribution_.HandleReferenceSequencesUntil(still_needed_reference_sequence, still_needed_position, thread_data.fragment_distribution_, *self.reference_, self.duplicates_, self.print_mutex_);
		}
		else{
			self.reading_success_ = false;
		}

		thread_data.rec_store_.clear();
	}

	if(0 == --self.running_threads_){
		{
			lock_guard<mutex> lock(self.finish_threads_mutex_);
			self.finish_threads_ = true;
		}
		self.finish_threads_cv_.notify_all();
	}
	else{
		unique_lock<mutex> lock(self.finish_threads_mutex_);
		self.finish_threads_cv_.wait(lock, [&self]{return self.finish_threads_;});
    }

	self.fragment_distribution_.FinishThreads(thread_data.fragment_distribution_, *self.reference_, self.duplicates_, self.print_mutex_);
	self.fragment_distribution_.AddThreadData(thread_data.fragment_distribution_);
}

DataStats::DataStats(Reference *ref, uint32_t maximum_insert_length, uint16_t minimum_mapping_quality):
	reference_(ref),
	batch_size_(10000), // 10,000 read pairs are read in each batch
	maximum_insert_length_(maximum_insert_length),
	minimum_mapping_quality_(minimum_mapping_quality), // For arabidopsis thaliana mapping with bowtie2 a strange clustering of fragments at same positions was observed for lower qualities than 10
	block_size_(1000),
	reading_success_(true),
	read_records_(0),
	minimum_quality_(255),
	maximum_quality_(0),
	minimum_read_length_on_reference_(numeric_limits<uint32_t>::max()),
	maximum_read_length_on_reference_(0),
	total_number_reads_(0),
	reads_in_unmapped_pairs_without_adapters_(0),
	reads_in_unmapped_pairs_with_adapters_(0),
	reads_with_low_quality_(0),
	reads_on_too_short_fragments_(0),
	reads_to_close_to_ref_seq_ends_(0),
	reads_used_(0){
	for( uint16_t template_segment=2; template_segment--; ){
		read_lengths_.at(template_segment).SetOffset(1); // As a read with length of 0 cannot be considered a read, this length can be excluded right from the start
	}
	if(ref){
		min_dist_to_ref_seq_ends_ = ref->MinDistToRefSeqEnds();
	}
}

bool DataStats::IsValidRecord(const BamAlignmentRecord &record){
	// Non-valid with warnings
	if( !hasFlagMultiple(record) ){
		printErr << "Read '" << record.qName << "' is not paired.\n";
		return false;
	}
	if( hasFlagFirst(record) == hasFlagLast(record) ){
		printErr << "Read '" << record.qName << "' is either first and second in its template or none of it.\n";
		return false;
	}
	if(0==length(record.seq)){
		printErr << "Read '" << record.qName << "' does not have the sequence stored.\n";
		return false;
	}
	if(length(record.seq)!=length(record.qual)){
		printErr << "Read '" << record.qName << "' does not have the same length for the sequence and its quality.\n";
		return false;
	}
	if( hasFlagSecondary(record) || hasFlagSupplementary(record) ){ // Only primary mappings are valid
		printErr << "Read '" << record.qName << "' is not a primary mapping.\n";
		return false;
	}

	return true;
}

uint32_t DataStats::GetReadLengthOnReference(const BamAlignmentRecord &record){
	uint32_t real_read_length = length(record.seq);

	// Update real_read_length based on insertion and deletion
	for( const auto &cigar_element : record.cigar ){
		switch( cigar_element.operation ){
		case 'N':
		case 'D':
			real_read_length += cigar_element.count;
			break;
		case 'I':
			real_read_length -= cigar_element.count;
			break;
		}
	}
	return real_read_length;
}

uint32_t DataStats::GetReadLengthOnReference(const BamAlignmentRecord &record, uint32_t &max_indel){
	uint32_t real_read_length = length(record.seq);

	// Update real_read_length based on insertion and deletion
	for( const auto &cigar_element : record.cigar ){
		switch( cigar_element.operation ){
		case 'N':
		case 'D':
			SetToMax(max_indel, cigar_element.count);
			real_read_length += cigar_element.count;
			break;
		case 'I':
			SetToMax(max_indel, cigar_element.count);
			real_read_length -= cigar_element.count;
			break;
		}
	}
	return real_read_length;
}

// Deactivation of pcr calculation only for speeding up of tests
bool DataStats::ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const string &variant_file, uint64_t max_ref_seq_bin_size, uint16_t num_threads, bool calculate_bias ){
	if(calculate_bias){
		fragment_distribution_.ActivateBiasCalculation();
	}
	else{
		fragment_distribution_.DeactivateBiasCalculation(); // Only for speeding up tests
	}

	bool success = true;
	BamFileIn bam;

	if(reference_){
		if( !open(bam, bam_file) ){
			printErr << "Could not open " << bam_file << " for reading.\n";
			success = false;
		}
		else if(atEnd(bam)){
			printErr << bam_file << " does not contain any sequences.\n";
			success = false;
		}
		else{
			BamHeader header;
			try{
				readHeader(header, bam);
			}
			catch(const Exception &e){
				printErr << "Could not read header in " << bam_file << ": " << e.what() << '\n';
				success = false;
			}

			if( success ){
				const FormattedFileContext<BamFileIn, void>::Type &bam_context = context(bam);

				if( length(contigNames(bam_context)) == reference_->NumberSequences() ){
					bool reference_ordering_identical = true;
					for(auto i = length(contigNames(bam_context)); i--; ){
						if( !(reference_->ReferenceIdFirstPart(i) == contigNames(bam_context)[i]) ){
							reference_ordering_identical = false;
							printErr << "The ordering of the reference sequences in the bam file is not identical to the reference file\n" << contigNames(bam_context)[i] << '\n' << reference_->ReferenceId(i) << '\n';
							success = false;
							break;
						}
					}

					if(reference_ordering_identical){
						uint8_t size_mapping_quality;
						uint32_t size_indel;
						if( PreRun(bam, bam_file, header, size_mapping_quality, size_indel, max_ref_seq_bin_size) ){
							if( adapters_.LoadAdapters( adapter_file, adapter_matrix, phred_quality_offset_, std::max(read_lengths_.at(0).to(), read_lengths_.at(1).to())) ){
								if(!variant_file.empty()){
									if(reference_->PrepareVariantFile(variant_file)){
										if(!reference_->ReadFirstVariantPositions()){
											success = false;
										}
									}
									else{
										success = false;
									}
								}

								if(success){
									PrepareReadIn(size_mapping_quality, size_indel);

									printInfo << "Starting main read-in\n";

									auto max_seq_bin_len = fragment_distribution_.MaxRefSeqBinLength(*reference_);
									thread threads[num_threads];
									running_threads_ = num_threads;
									finish_threads_ = false;
									for(auto i = num_threads; i--; ){
										threads[i] = thread(ReadThread, std::ref(*this), std::ref(bam), max_seq_bin_len);
									}
									for(auto i = num_threads; i--; ){
										threads[i].join();
									}

									if(reading_success_){
										printInfo << "Finished processing all reads.\n";
									}
									else{
										success = false;
									}
								}
							}
							else{
								success = false;
							}
						}
						else{
							success = false;
						}
					}
				}
				else{
					printErr << "The number of the reference sequences in the bam file is not identical to the number of sequences in the reference file\n";
					success = false;
				}
			}
		}
	}
	else{
		printErr << "No reference provided\n";
		success = false;
	}
	seqan::close(bam);

	reference_->ClearAllVariantPositions();

	if( 0 == reads_used_ && calculate_bias ){ // In case we are testing and not calculating the bias afterwards, we can continue running
		printErr << "No reads in the file passed all criteria to be used for the statistics\n";
		success = false;
	}

	if(!success || SignsOfPairsWithNamesNotIdentical() || !FinishReadIn() ){
		total_number_reads_ = 0; // Marking that the reading in was not successful
		return false;
	}

	uint64_t reads_in_wrong_place = total_number_reads_ - reads_in_unmapped_pairs_without_adapters_ - reads_in_unmapped_pairs_with_adapters_ - reads_with_low_quality_ - reads_on_too_short_fragments_ - reads_to_close_to_ref_seq_ends_ - reads_used_;
	printInfo << "Of the " << total_number_reads_ << " reads in the file\n";
	printInfo << reads_used_ << " (" << static_cast<uint16_t>(Percent(reads_used_, total_number_reads_)) << "\%) could be used for all statistics\n";
	printInfo << reads_with_low_quality_ << " (" << static_cast<uint16_t>(Percent(reads_with_low_quality_, total_number_reads_)) << "\%) had too low mapping qualities\n";
	printInfo << reads_on_too_short_fragments_ << " (" << static_cast<uint16_t>(Percent(reads_on_too_short_fragments_, total_number_reads_)) << "\%) are on reference sequences that are too short\n";
	printInfo << reads_to_close_to_ref_seq_ends_ << " (" << static_cast<uint16_t>(Percent(reads_to_close_to_ref_seq_ends_, total_number_reads_)) << "\%) mapped to close to the ends of a reference sequence\n";
	printInfo << reads_in_wrong_place << " (" << static_cast<uint16_t>(Percent(reads_in_wrong_place, total_number_reads_)) << "\%) were mapping too far apart from their partner or in wrong direction\n";
	printInfo << reads_in_unmapped_pairs_with_adapters_ << " (" << static_cast<uint16_t>(Percent(reads_in_unmapped_pairs_with_adapters_, total_number_reads_)) << "\%) were in an unmapped pair with adapters detected\n";
	printInfo << reads_in_unmapped_pairs_without_adapters_ << " (" << static_cast<uint16_t>(Percent(reads_in_unmapped_pairs_without_adapters_, total_number_reads_)) << "\%) were in an unmapped pair without adapters detected\n";

	Shrink();
	if(!Calculate(num_threads)){
		total_number_reads_ = 0; // Marking that the reading in was not successful
		return false;
	}

	return true;
}

bool DataStats::Load( const char *archive_file ){
	try{
		// create and open an archive for input
		ifstream ifs(archive_file);
		boost::archive::text_iarchive ia(ifs);

		// read class state from archive
		ia >> *this;
	}
	catch(const exception& e){
		printErr << "Could not load data statistics: " << e.what() << "\n";
		return false;
	}

	return true;
}

bool DataStats::Save( const char *archive_file ) const{
	try{
		// create all missing directories
		CreateDir(archive_file);

		// create and open a character archive for output
		ofstream ofs(archive_file);
		boost::archive::text_oarchive oa(ofs);

		// save data to archive
		oa << *this;
	}
	catch(const exception& e){
		printErr<< "Could not save data statistics: " << e.what() << "\n";
		return false;
	}

	return true;
}

void DataStats::PrepareProcessing(){
	PrepareGeneral();

	adapters_.PrepareSimulation();
	qualities_.PrepareEstimation();
	errors_.PrepareSimulation();
}

void DataStats::PreparePlotting(){
	PrepareGeneral();
	coverage_.PreparePlotting();
	fragment_distribution_.PreparePlotting();
	duplicates_.PreparePlotting();
	errors_.PreparePlotting();
	qualities_.PreparePlotting();
}

void DataStats::PrepareTesting(){
	PrepareGeneral();

	coverage_.PreparePlotting();
	duplicates_.PrepareTesting();
	errors_.PreparePlotting();
	qualities_.PrepareTesting();
}
