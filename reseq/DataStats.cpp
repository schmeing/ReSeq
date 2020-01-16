#include "DataStats.h"
using reseq::SeqQualityStats;
using reseq::DataStats;
using reseq::Vect;

#include <algorithm>
using std::max;
using std::min;
//include <array>
using std::array;
#include <cstring>
using std::strlen;
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
using seqan::contigNames;
using seqan::Dna;
using seqan::Dna5;
using seqan::Exception;
using seqan::readRecord;
using seqan::readHeader;

#include "CMakeConfig.h"
//include "utilities.hpp"
using reseq::utilities::ConstIupacStringReverseComplement;
using reseq::utilities::ReversedConstCharString;
using reseq::utilities::ReversedConstCigarString;
using reseq::utilities::Complement;
using reseq::utilities::at;
using reseq::utilities::CreateDir;
using reseq::utilities::Divide;
using reseq::utilities::FileExists;
using reseq::utilities::IsN;
using reseq::utilities::MeanWithRoundingToFirst;
using reseq::utilities::Percent;
using reseq::utilities::SafePercent;
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

inline bool DataStats::PotentiallyValidGeneral( const seqan::BamAlignmentRecord &record ) const{
	return !hasFlagNextUnmapped(record) && !hasFlagUnmapped(record) && record.rNextId == record.rID && !reference_->ReferenceSequenceExcluded(record.rID);
}
inline bool DataStats::PotentiallyValidFirst( const seqan::BamAlignmentRecord &record_first ) const{
	return PotentiallyValidGeneral(record_first) && record_first.beginPos + maximum_insert_length_ >= record_first.pNext;
}
inline bool DataStats::PotentiallyValidSecond( const seqan::BamAlignmentRecord &record_second ) const{
	return PotentiallyValidGeneral(record_second) && record_second.pNext + maximum_insert_length_ >= record_second.beginPos;
}
inline bool DataStats::PotentiallyValid( const seqan::BamAlignmentRecord &record ) const{
	if(record.beginPos < record.pNext){
		return PotentiallyValidFirst(record);
	}
	else{
		return PotentiallyValidSecond(record);
	}
}

bool DataStats::IsSecondRead( CoverageStats::FullRecord *record, CoverageStats::FullRecord *&record_first, CoverageStats::CoverageBlock *&block ){
	auto insert_info = first_read_records_.insert( record );
	if( insert_info.second ){
		// Insertion worked, so it is the first read
		if( PotentiallyValidFirst(record->record_) ){
			// Use record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially soft-clipped bases at the beginning
			coverage_.AddFragment( record->record_.rID, (record->record_.beginPos>maximum_read_length_on_reference_?record->record_.beginPos-maximum_read_length_on_reference_:0), block );
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
	uintReadLen adapter_position_first, adapter_position_second;
	bool adapter_detected = adapters_.Detect(adapter_position_first, adapter_position_second, record_first, record_second, *reference_);

	if(adapter_detected){
		auto insert_length = MeanWithRoundingToFirst(adapter_position_first,adapter_position_second);
		fragment_distribution_.AddInsertLengths(insert_length);
		++tmp_non_mapped_read_lengths_by_fragment_length_.at(0).at(insert_length).at(length(record_first.seq));
		++tmp_non_mapped_read_lengths_by_fragment_length_.at(1).at(insert_length).at(length(record_second.seq));

		return true;
	}
	else{
		return false;
	}
}

void DataStats::EvalBaseLevelStats( CoverageStats::FullRecord *full_record, uintTempSeq template_segment, uintTempSeq strand, uintTileId tile_id, uintQual &paired_seq_qual ){
	const BamAlignmentRecord &record(full_record->record_);

	uintSeqLen pos_reversed(length(record.seq));

	SeqQualityStats<uintNucCount> seq_qual_stats;
	seq_qual_stats[maximum_quality_]; // Resize to maximum quality
	for( auto qual : record.qual ){
		++seq_qual_stats.at(qual - phred_quality_offset_);
	}
	seq_qual_stats.Calculate();
	seq_qual_stats.CalculateProbabilityMean();
	full_record->sequence_quality_ = seq_qual_stats.mean_;

	Dna5 base;
	uintQual qual, last_qual(1);
	array<uintSeqLen, 5> read_bases = {0,0,0,0,0};
	uintReadLen homoquality_length = 0;
	uintReadLen homopolymer_length = 0;
	Dna5 homopolymer_nucleotide = 0;

	for( uintReadLen pos = 0; pos < length(record.seq); ++pos){
		// In case of reversed sequences base and qual are read in reverse from bamfile so pos represents the position in the fq file
		if( hasFlagRC(record) ){
			base = Complement::Dna5(at(record.seq, --pos_reversed));
			qual = at(record.qual, pos_reversed)-phred_quality_offset_;
		}
		else{
			base = at(record.seq, pos);
			qual = at(record.qual, pos)-phred_quality_offset_;
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

	++tmp_gc_read_content_.at(template_segment).at(Percent( read_bases.at(1)+read_bases.at(2), static_cast<uintSeqLen>(length(record.seq)) ));
	++tmp_n_content_.at(template_segment).at(Percent( read_bases.at(4), static_cast<uintSeqLen>(length(record.seq)) ));
}

bool DataStats::EvalReferenceStatistics(
		CoverageStats::FullRecord *record,
		uintTempSeq template_segment,
		CoverageStats::CoverageBlock *coverage_block ){
	// Get gc on reference
	auto cov_block = coverage_block; // coverage_block might be changed in the next step so keep the original for later
	uintSeqLen coverage_pos = CoverageStats::GetStartPos(record->from_ref_pos_, cov_block);
	array<uintNucCount, 5> seq_content_reference = {0,0,0,0,0};
	for( auto ref_pos = record->from_ref_pos_; ref_pos < record->to_ref_pos_; ++ref_pos ){
		if( cov_block->coverage_.at(ref_pos-cov_block->start_pos_).valid_ ){
			++seq_content_reference.at( at(reference_->ReferenceSequence(record->record_.rID), ref_pos) );
		}
		else{
			++seq_content_reference.at(4);
		}
		CoverageStats::IncrementPos(coverage_pos, cov_block);
	}
	auto gc_percent = SafePercent( seq_content_reference.at(1)+seq_content_reference.at(2), record->to_ref_pos_-record->from_ref_pos_-seq_content_reference.at(4) );
	record->reference_gc_ = gc_percent;

	uintSeqLen ref_pos;
	if( hasFlagRC(record->record_) ){
		coverage_pos = CoverageStats::GetStartPos(record->to_ref_pos_ - 1, coverage_block);
		ref_pos = record->to_ref_pos_-1;
	}
	else{
		coverage_pos = CoverageStats::GetStartPos(record->from_ref_pos_, coverage_block);
		ref_pos = record->from_ref_pos_;
	}

	ConstIupacStringReverseComplement reversed_seq(record->record_.seq);
	ReversedConstCharString reversed_qual(record->record_.qual);
	ReversedConstCigarString rev_cigar(record->record_.cigar);

	uintReadLen read_pos(0);
	Dna ref_base;
	Dna5 base;
	uintBaseCall last_base(5);
	uintQual qual;

	array<uintNucCount, 5> seq_content_mapped = {0,0,0,0,0};
	uintReadLen num_errors(0), indel_pos(0);
	uintInDelType indel_type(0);

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
						qual = at(reversed_qual,read_pos)-phred_quality_offset_;
						ref_base = Complement::Dna5(at(reference_->ReferenceSequence(record->record_.rID), ref_pos));
					}
					else{
						base = at(record->record_.seq, read_pos);
						qual = at(record->record_.qual, read_pos)-phred_quality_offset_;
						ref_base = at(reference_->ReferenceSequence(record->record_.rID), ref_pos);
					}

					if( !coverage_block->coverage_.at(coverage_pos).valid_ ){
						++seq_content_mapped.at(4);
					}
					else{
						++tmp_sequence_content_reference_.at(template_segment).at(hasFlagRC(record->record_)).at(ref_base).at(read_pos);

						errors_.AddBasePlotting(template_segment, ref_base, base, qual, last_base);
						errors_.AddInDel( indel_type, last_base, ErrorStats::kNoInDel, indel_pos, read_pos, gc_percent);

						++seq_content_mapped.at(base);

						if(ref_base != base){
							++num_errors;
						}
					}

					last_base = base;

					++read_pos;

					if( hasFlagRC(record->record_) ){
						coverage_.AddReverse(coverage_pos, coverage_block, base);

						CoverageStats::DecrementPos(coverage_pos, coverage_block);
						--ref_pos;
					}
					else{
						coverage_.AddForward(coverage_pos, coverage_block, base);

						CoverageStats::IncrementPos(coverage_pos, coverage_block);
						++ref_pos;
					}

					indel_type = 0;
					indel_pos = 0;
				}
			}
			break;
		case 'N':
		case 'D':
			for(indel_pos=0; indel_pos<cigar_element.count; ++indel_pos ){
				if(record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_){
					// We are currently not comparing the adapter to the reference
					if( hasFlagRC(record->record_) ){
						ref_base = Complement::Dna5(at(reference_->ReferenceSequence(record->record_.rID), ref_pos));
					}
					else{
						ref_base = at(reference_->ReferenceSequence(record->record_.rID), ref_pos);
					}

					if( coverage_block->coverage_.at(coverage_pos).valid_ ){
						errors_.AddInDel( indel_type, last_base, ErrorStats::kDeletion, indel_pos, read_pos, gc_percent);

						++num_errors;
					}

					if( hasFlagRC(record->record_) ){
						CoverageStats::DecrementPos(coverage_pos, coverage_block);
						--ref_pos;
					}
					else{
						CoverageStats::IncrementPos(coverage_pos, coverage_block);
						++ref_pos;
					}
					indel_type = 1;
				}
			}
			break;
		case 'I':
			if(record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_){
				if( hasFlagRC(record->record_) ){
					ref_base = Complement::Dna5(at(reference_->ReferenceSequence(record->record_.rID), ref_pos));
				}
				else{
					ref_base = at(reference_->ReferenceSequence(record->record_.rID), ref_pos);
				}

				if( !coverage_block->coverage_.at(coverage_pos).valid_ ){
					read_pos += cigar_element.count;
					indel_type = 0;
				}
				else{
					for(indel_pos=0; indel_pos<cigar_element.count; ++indel_pos ){
						if( hasFlagRC(record->record_) ){
							base = at(reversed_seq, read_pos);
						}
						else{
							base = at(record->record_.seq, read_pos);
						}

						errors_.AddInDel( indel_type, last_base, static_cast<ErrorStats::InDelDef>(static_cast<uintBaseCall>(base)+2), indel_pos, read_pos, gc_percent);
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

bool DataStats::EvalRecord( pair<CoverageStats::FullRecord *, CoverageStats::FullRecord *> record, ThreadData &thread ){
	// Start handling pair by determining tile(and counting it) and template segment
	uintTileId tile_id;
	if( !tiles_.GetTileId(tile_id, record.first->record_.qName) ){
		return false;
	}

	uintTempSeq template_segment(0);
	if( hasFlagLast(record.second->record_) ){
		template_segment = 1;
	}

	uintTempSeq strand;
	if(hasFlagUnmapped(record.first->record_)){
		if(hasFlagUnmapped(record.second->record_)){
			strand = 2; // Invalid strand as we cannot know it without mapping
		}
		else{
			strand = (hasFlagRC(record.second->record_) && hasFlagFirst(record.second->record_)) || (!hasFlagRC(record.second->record_) && hasFlagLast(record.second->record_));
		}
	}
	else{
		strand = (hasFlagRC(record.first->record_) && hasFlagFirst(record.first->record_)) || (!hasFlagRC(record.first->record_) && hasFlagLast(record.first->record_));
	}

	// Base level stats
	uintQual paired_seq_qual(0);
	EvalBaseLevelStats( record.first, (template_segment+1)%2, strand, tile_id, paired_seq_qual );
	EvalBaseLevelStats( record.second, template_segment, strand, tile_id, paired_seq_qual );

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
			if( reference_->ReferenceSequenceExcluded(record.first->record_.rID) || reference_->ReferenceSequenceExcluded(record.second->record_.rID) ){
				reads_on_too_short_fragments_ += 2;
			}
			else{
				uintSeqLen start_pos_first, end_pos_first, start_pos_second, end_pos_second;
				GetReadPosOnReference(start_pos_first, end_pos_first, record.first->record_);
				GetReadPosOnReference(start_pos_second, end_pos_second, record.second->record_);
				if( reference_->FragmentExcluded( thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_, record.first->record_.rID, start_pos_first, end_pos_first ) || reference_->FragmentExcluded( thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_, record.second->record_.rID, start_pos_second, end_pos_second ) ){
					reads_in_excluded_regions_ += 2;
				}
				else{
					if( CheckForAdapters(record.first->record_, record.second->record_) ){ // We don't trust the mapping, so look for adapters like in unmapped case
						reads_with_low_quality_with_adapters_ += 2;
					}
					else{
						reads_with_low_quality_without_adapters_ += 2;
					}
				}
			}
		}
		else{
			if( record.first->record_.rID == record.second->record_.rID ){
				if( reference_->ReferenceSequenceExcluded(record.first->record_.rID) ){
					reads_on_too_short_fragments_ += 2;
				}
				else{
					uintSeqLen start_pos_first, start_pos_second, end_pos_first, end_pos_second;
					GetReadPosOnReference(start_pos_first, end_pos_first, record.first->record_);
					GetReadPosOnReference(start_pos_second, end_pos_second, record.second->record_);

					if(start_pos_second < start_pos_first){
						// Switch records as the soft-trimming changed order in bam sort
						auto tmp = record.first;
						record.first = record.second;
						record.second = tmp;

						template_segment = (template_segment+1)%2;

						auto tmp2 = start_pos_first;
						start_pos_first = start_pos_second;
						start_pos_second = tmp2;
						tmp2 = end_pos_first;
						end_pos_first = end_pos_second;
						end_pos_second = tmp2;
					}

					bool adapter_detected(false);
					if(hasFlagRC(record.first->record_) && !hasFlagRC(record.second->record_) && end_pos_first > start_pos_second && end_pos_first - start_pos_second > (end_pos_first - start_pos_first)/2 ){
						adapter_detected = true;
					}

					if(InProperDirection( record.first->record_, end_pos_second, start_pos_first, start_pos_second ) || adapter_detected){
						if( (!adapter_detected && reference_->FragmentExcluded( thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_, record.first->record_.rID, start_pos_first, end_pos_second )) || (adapter_detected && reference_->FragmentExcluded( thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_, record.first->record_.rID, start_pos_second, end_pos_first )) ){
							reads_in_excluded_regions_ += 2;
						}
						else{
							reads_used_ += 2;

							if(adapter_detected ){
								uintReadLen adapter_position_first, adapter_position_second;
								// Reads in proper direction for adapters (reverse-forward) and overlap is large enough (minimum half of the read length)
								adapters_.Detect(adapter_position_first, adapter_position_second, record.first->record_, record.second->record_, *reference_, true); // Call it just to insert the adapter stats, position is better known from mapping
							}

							CoverageStats::CoverageBlock *coverage_block;
							// Pair stats
							uintSeqLen insert_length;
							if(adapter_detected){
								insert_length = end_pos_first - start_pos_second;

								fragment_distribution_.AddGCContent(reference_->GCContent(record.second->record_.rID, start_pos_second, end_pos_first));
								fragment_distribution_.FillInOutskirtContent(*reference_, record.second->record_, start_pos_second, end_pos_first); // Use second record here as it is forward, so if it is sequenced second, fragment is reversed
								coverage_block = coverage_.FindBlock( record.second->record_.rID, start_pos_second );

								record.first->from_ref_pos_ = start_pos_second;
								record.second->from_ref_pos_ = start_pos_second;
								record.first->to_ref_pos_ = end_pos_first;
								record.second->to_ref_pos_ = end_pos_first;

								fragment_distribution_.AddFragmentSite(record.first->record_.rID, insert_length, record.first->from_ref_pos_, template_segment, *reference_, thread.fragment_distribution_); // If second record(must be forward to be accepted due to adapter) is sequenced second, fragment is reversed
							}
							else{
								insert_length = end_pos_second - start_pos_first;

								fragment_distribution_.AddGCContent(reference_->GCContent(record.first->record_.rID, start_pos_first, end_pos_second));
								fragment_distribution_.FillInOutskirtContent(*reference_, record.first->record_, start_pos_first, end_pos_second);
								coverage_block = coverage_.FindBlock( record.first->record_.rID, start_pos_first );

								record.first->from_ref_pos_ = start_pos_first;
								record.second->from_ref_pos_ = start_pos_second;
								record.first->to_ref_pos_ = end_pos_first;
								record.second->to_ref_pos_ = end_pos_second;

								fragment_distribution_.AddFragmentSite(record.first->record_.rID, insert_length, record.first->from_ref_pos_, (template_segment+1)%2, *reference_, thread.fragment_distribution_); // If first record(must be forward to be accepted without adapters) is sequenced second, fragment is reversed
							}

							fragment_distribution_.AddAbundance(record.first->record_.rID);
							fragment_distribution_.AddInsertLengths(insert_length); // We do not just use the template length from field 9 of the bam file, because bwa stores there the right-most position of the reverse read - the left-most position of the forward read, which is not the fragment length in the case of the reverse read being positioned before the forward read. Bowtie2 gives the correct length here, but I don't want to sacrifice compatibility so easily

							++tmp_read_lengths_by_fragment_length_.at(0).at(insert_length).at(length(record.first->record_.seq));
							++tmp_read_lengths_by_fragment_length_.at(1).at(insert_length).at(length(record.second->record_.seq));

							record.first->tile_id_ = tile_id;
							record.second->tile_id_ = tile_id;
							record.first->fragment_length_ = insert_length;
							record.second->fragment_length_ = insert_length;

							if( !EvalReferenceStatistics( record.first, (template_segment+1)%2, coverage_block ) ){
								return false;
							}
							if( !EvalReferenceStatistics( record.second, template_segment, coverage_block ) ){
								return false;
							}
						}
					}
				}
			}
		}
	}

	return true;
}

bool DataStats::SignsOfPairsWithNamesNotIdentical(){
	// All reads have been processed, so first_read_records_ must be empty
	if( first_read_records_.size() ){
		printErr << "Read names of reads from the same pair are not identical. Cannot continue.\nExample position: " << (*first_read_records_.begin())->record_.rID << ":" << (*first_read_records_.begin())->record_.beginPos << "\nExample read name: " << (*first_read_records_.begin())->record_.qName << std::endl;
		return true;
	}

	return false;
}

void DataStats::PrepareReadIn(uintQual size_mapping_quality, uintReadLen size_indel, uintSeqLen max_ref_seq_bin_size){
	uintReadLen size_pos = max(read_lengths_.at(0).to(), read_lengths_.at(1).to());

	uintFragCount num_reads(0);
	uintNucCount num_bases(0);
	for( int template_segment = 2; template_segment--; ){
		for(auto len=read_lengths_.at(template_segment).from(); len < read_lengths_.at(template_segment).to(); ++len ){
			num_reads += read_lengths_.at(template_segment).at(len);
			num_bases += read_lengths_.at(template_segment).at(len)*len;
		}
	}

	adapters_.PrepareAdapters(size_pos, phred_quality_offset_);
	coverage_.Prepare( Divide(num_bases,reference_->TotalSize()), Divide(num_bases,num_reads), maximum_read_length_on_reference_ );
	duplicates_.PrepareTmpDuplicationVector(maximum_insert_length_);
	errors_.Prepare(tiles_.NumTiles(), maximum_quality_+1, size_pos, size_indel);
	fragment_distribution_.Prepare(*reference_, maximum_insert_length_, max_ref_seq_bin_size, reads_per_frag_len_bin_, lowq_reads_per_frag_len_bin_);
	reads_per_frag_len_bin_.clear();
	reads_per_frag_len_bin_.shrink_to_fit();
	lowq_reads_per_frag_len_bin_.clear();
	lowq_reads_per_frag_len_bin_.shrink_to_fit();
	qualities_.Prepare(tiles_.NumTiles(), maximum_quality_+1, size_pos, maximum_insert_length_);

	// Prepare vector in this class
	tmp_proper_pair_mapping_quality_.resize(size_mapping_quality);
	tmp_improper_pair_mapping_quality_.resize(size_mapping_quality);
	tmp_single_read_mapping_quality_.resize(size_mapping_quality);

	for( auto template_segment=2; template_segment--; ){
		SetDimensions( tmp_read_lengths_by_fragment_length_.at(template_segment), maximum_insert_length_+1, size_pos );
		SetDimensions( tmp_non_mapped_read_lengths_by_fragment_length_.at(template_segment), maximum_insert_length_+1, size_pos );

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
	adapters_.Finalize(); // Collect ambigous adapters into a single one
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
		non_mapped_read_lengths_by_fragment_length_.at(template_segment).Acquire(tmp_non_mapped_read_lengths_by_fragment_length_.at(template_segment));
		ShrinkVect(non_mapped_read_lengths_by_fragment_length_.at(template_segment));
		for(auto frag_len = non_mapped_read_lengths_by_fragment_length_.at(template_segment).from(); frag_len < non_mapped_read_lengths_by_fragment_length_.at(template_segment).to(); ++frag_len){
			for(auto read_len = non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).from(); read_len < non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).to(); ++read_len){
				tmp_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).at(read_len) += non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).at(read_len);
			}
		}
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
	qualities_.Shrink();
	tiles_.Shrink();
	
	for( uintTempSeq template_segment=2; template_segment--; ){
		ShrinkVect(read_lengths_by_fragment_length_.at(template_segment));
	}
}

// Deactivation of bias calculation only for speeding up of tests
bool DataStats::Calculate(uintNumThreads num_threads){
	// Calculate the duplicates from observation of fragments at same position and reference characteristics
	if(!fragment_distribution_.FinalizeBiasCalculation(*reference_, num_threads, duplicates_)){
		return false;
	}

	// Calculate coverage corrected for low quality sites
	uintFragCount num_reads(0);
	uintNucCount num_bases(0);
	for( int template_segment = 2; template_segment--; ){
		for(auto len=read_lengths_.at(template_segment).from(); len < read_lengths_.at(template_segment).to(); ++len ){
			num_reads += read_lengths_.at(template_segment).at(len);
			num_bases += read_lengths_.at(template_segment).at(len)*len;
		}
	}
	corrected_coverage_ = fragment_distribution_.CorrectedCoverage(*reference_, static_cast<double>(num_bases)/num_reads); // Must be called after FinalizeBiasCalculation, where the coverage is corrected
	printInfo << "Estimated a corrected mean coverage of " << corrected_coverage_ << "x" << std::endl;

	return true;
}

void DataStats::PrepareGeneral(){
	// Set total_number_reads_ to the number of accepted reads as it is the number of read records after read in
	total_number_reads_ = SumVect(read_lengths_.at(0)) + SumVect(read_lengths_.at(1));

	adapters_.SumCounts();
}

bool DataStats::OrderOfBamFileCorrect(const seqan::BamAlignmentRecord &record, pair< uintRefSeqId, uintSeqLen > last_record_pos){
	if( !hasFlagUnmapped(record) ){
		if( record.rID < last_record_pos.first || (record.rID == last_record_pos.first && record.beginPos < last_record_pos.second) ){
			printErr << "Bamfile is not sorted by position. Detected at position " << record.beginPos << " with read name: " << record.qName << std::endl;
			return false;
		}
		last_record_pos = {record.rID, record.beginPos};
	}

	return true;
}

bool DataStats::PreRun(BamFileIn &bam, const char *bam_file, BamHeader &header, uintQual &size_mapping_quality, uintReadLen &size_indel){
	printInfo << "Starting PreRun" << std::endl;

	bool error = false;
	BamAlignmentRecord record;
	pair< uintRefSeqId, uintSeqLen > last_record_pos{0,0};

	uintReadLen read_length_on_reference;
	uintQual min_mapq(255), max_mapq(0);
	size_indel = 0;

	uintRefLenCalc num_bins(0);
	for(auto ref_seq = reference_->NumberSequences(); ref_seq--; ){
		num_bins += reference_->SequenceLength(ref_seq)/maximum_insert_length_ + 1;
	}
	reads_per_frag_len_bin_.resize(num_bins);
	lowq_reads_per_frag_len_bin_.resize(num_bins);
	uintRefSeqId cur_ref_seq(0);
	uintRefLenCalc ref_seq_start_bin(0);

	adapters_.PrepareAdapterPrediction();
	do{
		try{
			++total_number_reads_; // Counter total_number_reads_ for read in, later in Calculate correct it with the number of accepted reads from read_length_

			readRecord(record, bam);

			if( IsValidRecord(record) ){
				if( OrderOfBamFileCorrect(record, last_record_pos) ){
					if( !hasFlagSecondary(record) && !hasFlagSupplementary(record) ){ // Ignore supplementary reads
						if(hasFlagFirst(record)){
							tiles_.EnterTile(record.qName); // Names of records in a pair are identical, so tile must only be identified once

							++read_lengths_.at(0)[ length(record.qual) ];
						}
						else{
							++read_lengths_.at(1)[ length(record.qual) ];
						}

						// Determine read length range on reference
						if(PotentiallyValid(record)){
							read_length_on_reference = GetReadLengthOnReference(record, size_indel);
							SetToMax(maximum_read_length_on_reference_, read_length_on_reference);
							SetToMin(minimum_read_length_on_reference_, read_length_on_reference);

							while(record.rID > cur_ref_seq){
								ref_seq_start_bin += reference_->SequenceLength(cur_ref_seq)/maximum_insert_length_ + 1;
								++cur_ref_seq;
							}
							++reads_per_frag_len_bin_.at( ref_seq_start_bin + record.beginPos/maximum_insert_length_ );
						}
						if( !hasFlagUnmapped(record) && record.mapQ < minimum_mapping_quality_ ){
							while(record.rID > cur_ref_seq){
								ref_seq_start_bin += reference_->SequenceLength(cur_ref_seq)/maximum_insert_length_ + 1;
								++cur_ref_seq;
							}
							++lowq_reads_per_frag_len_bin_.at( ref_seq_start_bin + record.beginPos/maximum_insert_length_ );
						}

						// Determine quality range
						for( auto i=length(record.qual); i--;){
							SetToMax(maximum_quality_, at(record.qual, i));
							SetToMin(minimum_quality_, at(record.qual, i));
						}

						SetToMax(max_mapq, record.mapQ);
						SetToMin(min_mapq, record.mapQ);

						adapters_.ExtractAdapterPart(record);
					}
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
			printErr << "Could not read record " << total_number_reads_ << " in " << bam_file << ": " << e.what() << std::endl;
		}
	}while( !atEnd(bam) );

	if(error){
		return false;
	}

	// Jump back to beginning of bam file
	seqan::close(bam);
	if( !open(bam, bam_file) ){
		printErr << "Could not open " << bam_file << " for reading." << std::endl;
		return false;
	}

	try{
		readHeader(header, bam);
	}
	catch(const Exception &e){
		printErr << "Could not read header in " << bam_file << ": " << e.what() << std::endl;
		return false;
	}

	for( uintTempSeq template_segment=2; template_segment--; ){
		read_lengths_.at(template_segment).Shrink();
	}

	if( minimum_quality_ < 64 ){
		if( minimum_quality_ < 33 ){
			printErr << "Minimum quality is " << minimum_quality_ << ", which is lower than the 33 from Sanger encoding. Are those values correct?" << std::endl;
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
			<< "Quality: " << static_cast<uintQualPrint>(minimum_quality_) << " - " << static_cast<uintQualPrint>(maximum_quality_) << " ( " << static_cast<uintQualPrint>(phred_quality_offset_) << "-based encoding )\n"
			<< "Mapping Quality: " << static_cast<uintQualPrint>(min_mapq) << " - " << static_cast<uintQualPrint>(max_mapq) << '\n'
			<< "Maximum InDel Length: " << size_indel << std::endl;

	size_mapping_quality = max_mapq+1;
	++size_indel; // Add one to go from index of last value to size

	return adapters_.PredictAdapters();
}

bool DataStats::ReadRecords( BamFileIn &bam, bool &not_done, ThreadData &thread_data ){
	CoverageStats::FullRecord *record;
	CoverageStats::CoverageBlock *cov_block(NULL);
	lock_guard<mutex> lock(read_mutex_);
	try{
		while( thread_data.rec_store_.size() < kBatchSize && !atEnd(bam) && reading_success_){
			++read_records_;

			record = new CoverageStats::FullRecord;
			readRecord(record->record_, bam);

			if( hasFlagSecondary(record->record_) || hasFlagSupplementary(record->record_) ){
				delete record;
			}
			else{
				if( PotentiallyValid(record->record_) ){
					// Use record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially soft-clipped bases at the beginning
					if( !coverage_.EnsureSpace(record->record_.rID, (record->record_.beginPos>maximum_read_length_on_reference_?record->record_.beginPos-maximum_read_length_on_reference_:0), record->record_.beginPos+maximum_read_length_on_reference_, record, *reference_) ){
						return false;
					}
				}

				// Add low q sites
				if( !QualitySufficient( record->record_ ) && !hasFlagUnmapped(record->record_) && !reference_->ReferenceSequenceExcluded(record->record_.rID) && // low quality and not excluded
					(hasFlagNextUnmapped(record->record_) || record->record_.rID != record->record_.rNextId || (!hasFlagRC(record->record_) && (record->record_.beginPos < record->record_.pNext || record->record_.beginPos > record->record_.pNext+length(record->record_.seq))) || (hasFlagRC(record->record_) && (record->record_.beginPos > record->record_.pNext || record->record_.beginPos+length(record->record_.seq) < record->record_.pNext)) ) ){ // not in an adapter pair
					uintSeqLen start_pos, end_pos;
					GetReadPosOnReference(start_pos, end_pos, record->record_);

					if( !reference_->FragmentExcluded( thread_data.last_exclusion_region_id_, thread_data.last_exclusion_ref_seq_, record->record_.rID, start_pos, end_pos ) ){
						if(hasFlagRC(record->record_)){
							fragment_distribution_.AddLowQSiteEnd(record->record_.rID, end_pos, *reference_, thread_data.fragment_distribution_);
						}
						else{
							fragment_distribution_.AddLowQSiteStart(record->record_.rID, start_pos, *reference_, thread_data.fragment_distribution_);
						}
					}
				}

				// Work with the pairs and not individual reads, so first combine reads to pairs, by storing first and then processing both reads at the same time
				CoverageStats::FullRecord *record_first;
				if( IsSecondRead(record, record_first, cov_block) ){
					thread_data.rec_store_.emplace_back(record_first, record);
				}
			}

			if( total_number_reads_ > kBatchSize && !(read_records_%(total_number_reads_/20)) ){
				lock_guard<mutex> lock(print_mutex_);
				printInfo << "Read " << static_cast<uintPercentPrint>(Percent(read_records_, total_number_reads_)) << "% of the reads." << std::endl;
			}
		}
		not_done = !atEnd(bam);
	}
	catch(const Exception &e){
		lock_guard<mutex> lock(print_mutex_);
		printErr << "Could not read record " << read_records_ << ": " << e.what() << std::endl;
		return false;
	}

	return true;
}

void DataStats::ReadThread( DataStats &self, BamFileIn &bam, uintSeqLen max_seq_bin_len ){
	CoverageStats::CoverageBlock *cov_block;
	uintFragCount processed_fragments(0);
	bool not_done(true);
	auto first_ref_seq = self.reference_->FirstNotExcludedSequence();
	if( first_ref_seq >= self.reference_->NumberSequences() ){
		printErr << "No contig longer than " << self.maximum_insert_length_ << std::endl;
		self.reading_success_ = false;
		first_ref_seq = 0; // Just so that the next command does not crash, before we shut down
	}
	ThreadData thread_data( self.maximum_insert_length_, max_seq_bin_len, self.reference_->StartExclusion(first_ref_seq) );
	thread_data.rec_store_.reserve(self.kBatchSize);
	uintRefSeqId still_needed_reference_sequence(0);
	uintSeqLen still_needed_position(0);
	while(not_done && self.reading_success_){
		// ReadRecords:
		// Reads in records and bundles them to batches of paired records
		// Registers reads as unprocessed in first coverage block, to make sure blocks are not removed before all reads in them are handled
		// Adds pointers of potentially valid reads to the last block they potentially overlap (read length on reference yet unknown), so the coverage part can be handled after all coverage information are gathered
		if( self.ReadRecords(bam, not_done, thread_data) ){
			cov_block = NULL;
			// Process batch
			for( auto &rec : thread_data.rec_store_ ){
				// EvalRecord:
				// Apply multiple filters to check whether fragments are valid
				// Set to_ref_pos for valid reads, that to_ref_pos==0 is used in coverage_.CleanUp to remove the non-valid reads
				// Calculate all the statistics that are independent of other reads
				// Fills the coverage information in the coverage blocks
				if(self.EvalRecord(rec, thread_data)){
					if( self.PotentiallyValidFirst(rec.first->record_) ){
						// coverage_.RemoveFragment
						// Marks the reads as processed in the first block of first read, so the fully processed coverage blocks can be handled after all reads are processed
						// Use record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially soft-clipped bases at the beginning
						self.coverage_.RemoveFragment( rec.first->record_.rID, (rec.first->record_.beginPos>self.maximum_read_length_on_reference_?rec.first->record_.beginPos-self.maximum_read_length_on_reference_:0), cov_block, processed_fragments );
					}
					else{
						// Remove reads that have not been added to a coverage block
						delete rec.first;
						delete rec.second;
					}
				}
				else{
					self.reading_success_ = false;
				}
			}

			if(cov_block){ // To make sure we don't have a whole batch of non-valid reads as coverage_.CleanUp already handled blocks at the beginning and we need at least one block still existing to keep track of the position
				// coverage_.CleanUp:
				// Once all reads in it have been handled processes the coverage blocks: The coverage itself and all the read statistics that depend on the other reads (systematic error)
				// Removes the processed blocks
				still_needed_reference_sequence = self.coverage_.CleanUp(still_needed_position, *self.reference_, self.qualities_, self.errors_, self.phred_quality_offset_, cov_block, processed_fragments, thread_data.coverage_);
			}

			self.coverage_.PreLoadVariants(*self.reference_);

			self.fragment_distribution_.HandleReferenceSequencesUntil(still_needed_reference_sequence, still_needed_position, thread_data.fragment_distribution_, *self.reference_, self.duplicates_, self.print_mutex_);
		}
		else{
			self.reading_success_ = false;
		}

		thread_data.rec_store_.clear();
	}

	if(self.reading_success_){
		if(0 == --self.running_threads_){
			{
				lock_guard<mutex> lock(self.finish_threads_mutex_);
				self.finish_threads_ = true;
			}
			self.finish_threads_cv_.notify_all();

			// This finalization might take a bit of time so do it already here
			if( !self.coverage_.Finalize(*self.reference_, self.qualities_, self.errors_, self.phred_quality_offset_, self.print_mutex_, thread_data.coverage_) ){
				self.reading_success_ = false;
			}
		}
		else{
			unique_lock<mutex> lock(self.finish_threads_mutex_);
			self.finish_threads_cv_.wait(lock, [&self]{return self.finish_threads_;});
		}

		if(self.reading_success_){
			self.fragment_distribution_.FinishThreads(thread_data.fragment_distribution_, *self.reference_, self.duplicates_, self.print_mutex_);
			self.fragment_distribution_.AddThreadData(thread_data.fragment_distribution_);
		}
	}
}

DataStats::DataStats(Reference *ref, uintSeqLen maximum_insert_length, uintQual minimum_mapping_quality):
	reference_(ref),
	maximum_insert_length_(maximum_insert_length),
	minimum_mapping_quality_(minimum_mapping_quality), // For arabidopsis thaliana mapping with bowtie2 a strange clustering of fragments at same positions was observed for lower qualities than 10
	reading_success_(true),
	read_records_(0),
	minimum_quality_(255),
	maximum_quality_(0),
	minimum_read_length_on_reference_(numeric_limits<uintReadLen>::max()),
	maximum_read_length_on_reference_(0),
	total_number_reads_(0),
	reads_in_unmapped_pairs_without_adapters_(0),
	reads_in_unmapped_pairs_with_adapters_(0),
	reads_with_low_quality_with_adapters_(0),
	reads_with_low_quality_without_adapters_(0),
	reads_on_too_short_fragments_(0),
	reads_in_excluded_regions_(0),
	reads_used_(0){
	for( uintTempSeq template_segment=2; template_segment--; ){
		read_lengths_.at(template_segment).SetOffset(1); // As a read with length of 0 cannot be considered a read, this length can be excluded right from the start
	}
}

bool DataStats::IsValidRecord(const BamAlignmentRecord &record){
	if( !hasFlagMultiple(record) ){
		printErr << "Read '" << record.qName << "' is not paired." << std::endl;
		return false;
	}
	if( hasFlagFirst(record) == hasFlagLast(record) ){
		printErr << "Read '" << record.qName << "' is either first and second in its template or none of it." << std::endl;
		return false;
	}
	if(0==length(record.seq)){
		printErr << "Read '" << record.qName << "' does not have the sequence stored." << std::endl;
		return false;
	}
	if(length(record.seq)!=length(record.qual)){
		printErr << "Read '" << record.qName << "' does not have the same length for the sequence and its quality." << std::endl;
		return false;
	}

	return true;
}

reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord &record){
	uintReadLen real_read_length = length(record.seq);

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

reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord &record, uintReadLen &max_indel){
	uintReadLen real_read_length = length(record.seq);

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

inline void DataStats::GetReadPosOnReference(uintSeqLen &start_pos, uintSeqLen &end_pos, const BamAlignmentRecord &record) const{
	end_pos = GetReadLengthOnReference(record); // start_pos will be added after it is corrected for soft-clipping
	start_pos = record.beginPos;

	if('S' == at(record.cigar, 0).operation){
		if(at(record.cigar, 0).count < start_pos){
			start_pos -= at(record.cigar, 0).count;
		}
		else{
			start_pos = 0;
		}
	}
	else if('H' == at(record.cigar, 0).operation){
		end_pos -= at(record.cigar, 0).count;

		if( 2 <= length(record.cigar) && 'S' == at(record.cigar, 1).operation){
			if(at(record.cigar, 1).count < start_pos){
				start_pos -= at(record.cigar, 1).count;
			}
			else{
				start_pos = 0;
			}
		}
	}

	end_pos += start_pos;

	if( 'H' == at(record.cigar, length(record.cigar)-1).operation && 2 <= length(record.cigar) ){
		end_pos -= at(record.cigar, length(record.cigar)-1).count;
	}

	if( reference_->SequenceLength( record.rID ) < end_pos ){
		end_pos = reference_->SequenceLength( record.rID ); // Can occur due to soft-clipping at the end
	}

	return;
}

// Deactivation of pcr calculation only for speeding up of tests
bool DataStats::ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const string &variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads, bool calculate_bias ){
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
			printErr << "Could not open " << bam_file << " for reading." << std::endl;
			success = false;
		}
		else if(atEnd(bam)){
			printErr << bam_file << " does not contain any sequences." << std::endl;
			success = false;
		}
		else{
			BamHeader header;
			try{
				readHeader(header, bam);
			}
			catch(const Exception &e){
				printErr << "Could not read header in " << bam_file << ": " << e.what() << std::endl;
				success = false;
			}

			if( success ){
				const auto &bam_context = context(bam);

				if( length(contigNames(bam_context)) == reference_->NumberSequences() ){
					bool reference_ordering_identical = true;
					for(auto i = length(contigNames(bam_context)); i--; ){
						if( !(reference_->ReferenceIdFirstPart(i) == at(contigNames(bam_context), i)) ){
							reference_ordering_identical = false;
							printErr << "The ordering of the reference sequences in the bam file is not identical to the reference file\n" << at(contigNames(bam_context), i) << '\n' << reference_->ReferenceId(i) << std::endl;
							success = false;
							break;
						}
					}

					if(reference_ordering_identical){
						uintQual size_mapping_quality;
						uintReadLen size_indel;

						reference_->PrepareExclusionRegions();
						reference_->ObtainExclusionRegions(reference_->NumberSequences(), maximum_insert_length_); // At the end we need all of them together anyways, so it makes no sense to read them in one after another

						if( (0 == strlen(adapter_file) && 0 == strlen(adapter_matrix)) || adapters_.LoadAdapters(adapter_file, adapter_matrix) ){
							if( PreRun(bam, bam_file, header, size_mapping_quality, size_indel) ){
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
									PrepareReadIn(size_mapping_quality, size_indel, max_ref_seq_bin_size);

									printInfo << "Starting main read-in" << std::endl;

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
										printInfo << "Finished processing all reads." << std::endl;
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
					printErr << "The number of the reference sequences in the bam file is not identical to the number of sequences in the reference file" << std::endl;
					success = false;
				}
			}
		}
	}
	else{
		printErr << "No reference provided" << std::endl;
		success = false;
	}
	seqan::close(bam);

	reference_->ClearAllVariantPositions();

	if(!success){
		return false;
	}

	uintFragCount reads_in_wrong_place = total_number_reads_ - reads_in_unmapped_pairs_without_adapters_ - reads_in_unmapped_pairs_with_adapters_ - reads_with_low_quality_with_adapters_ - reads_with_low_quality_without_adapters_ - reads_on_too_short_fragments_ - reads_in_excluded_regions_ - reads_used_;
	printInfo << "Of the " << total_number_reads_ << " reads in the file" << std::endl;
	printInfo << reads_used_ << " (" << static_cast<uintPercentPrint>(Percent(reads_used_, total_number_reads_)) << "\%) could be used for all statistics" << std::endl;
	printInfo << reads_with_low_quality_with_adapters_ << " (" << static_cast<uintPercentPrint>(Percent(reads_with_low_quality_with_adapters_, total_number_reads_)) << "\%) had too low mapping qualities with adapters detected" << std::endl;
	printInfo << reads_with_low_quality_without_adapters_ << " (" << static_cast<uintPercentPrint>(Percent(reads_with_low_quality_without_adapters_, total_number_reads_)) << "\%) had too low mapping qualities without adapters detected" << std::endl;
	printInfo << reads_in_excluded_regions_ << " (" << static_cast<uintPercentPrint>(Percent(reads_in_excluded_regions_, total_number_reads_)) << "\%) mapped to excluded regions of the reference" << std::endl;
	printInfo << reads_in_wrong_place << " (" << static_cast<uintPercentPrint>(Percent(reads_in_wrong_place, total_number_reads_)) << "\%) were mapping too far apart from their partner or in wrong direction" << std::endl;
	printInfo << reads_on_too_short_fragments_ << " (" << static_cast<uintPercentPrint>(Percent(reads_on_too_short_fragments_, total_number_reads_)) << "\%) are on reference sequences that are too short" << std::endl;
	printInfo << reads_in_unmapped_pairs_with_adapters_ << " (" << static_cast<uintPercentPrint>(Percent(reads_in_unmapped_pairs_with_adapters_, total_number_reads_)) << "\%) were in an unmapped pair with adapters detected" << std::endl;
	printInfo << reads_in_unmapped_pairs_without_adapters_ << " (" << static_cast<uintPercentPrint>(Percent(reads_in_unmapped_pairs_without_adapters_, total_number_reads_)) << "\%) were in an unmapped pair without adapters detected" << std::endl;

	if( 0 == reads_used_ && calculate_bias ){ // In case we are testing and not calculating the bias afterwards, we can continue running
		printErr << "No reads in the file passed all criteria to be used for the statistics" << std::endl;
		success = false;
	}

	if(!success || SignsOfPairsWithNamesNotIdentical() || !FinishReadIn() ){
		total_number_reads_ = 0; // Marking that the reading in was not successful
		return false;
	}

	Shrink();
	if(!Calculate(num_threads)){
		total_number_reads_ = 0; // Marking that the reading in was not successful
		return false;
	}

	reference_->ClearAllExclusionRegions(); // Don't remove any of it earlier because we still need all of them in Calculate

	creation_time_ = std::chrono::time_point_cast<std::chrono::seconds>( std::chrono::system_clock::now() ).time_since_epoch().count();

	return true;
}

bool DataStats::Load( const char *archive_file ){
	if( !FileExists(archive_file) ){
		printErr << "File '" << archive_file << "' does not exists or no read permission given." << std::endl;
		return false;
	}

	try{
		// create and open an archive for input
		ifstream ifs(archive_file);
		boost::archive::text_iarchive ia(ifs);

		// read class state from archive
		ia >> *this;
	}
	catch(const exception& e){
		printErr << "Could not load data statistics: " << e.what() << std::endl;
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
		printErr<< "Could not save data statistics: " << e.what() << std::endl;
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
