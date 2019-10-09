#include "AdapterStats.h"
using reseq::AdapterStats;

#include <algorithm>
using std::max;
using std::max_element;
using std::min;
#include <array>
using std::array;
#include <cmath>
using std::ceil;
#include <fstream>
using std::ifstream;
#include <list>
using std::list;
//include <string>
using std::string;
//include <vector>
using std::vector;
#include <utility>
using std::pair;

#include "reportingUtils.hpp"

//include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::CharString;
using seqan::Dna5;
using seqan::DnaString;
using seqan::Exception;
using seqan::reverseComplement;
using seqan::StringSet;

#include <seqan/seq_io.h>
using seqan::readRecords;
using seqan::SeqFileIn;

#include "skewer/src/matrix.h"

//include "utilities.hpp"
using reseq::utilities::SetToMax;

reseq::uintReadLen AdapterStats::CountErrors(uintReadLen &last_error_pos, const BamAlignmentRecord &record, const Reference &reference){
	uintReadLen num_errors(0), read_pos(0);
	auto ref_pos = record.beginPos;
	auto &ref_seq = reference.ReferenceSequence(record.rID);

	for( const auto &cigar_element : record.cigar ){
		switch( cigar_element.operation ){
		case 'M':
		case '=':
		case 'X':
			for(auto i=cigar_element.count; i--; ){
				if( ref_seq[ref_pos++] != record.seq[read_pos++] ){
					if( !hasFlagRC(record) || !num_errors ){ // In case of a normal read always overwrite errors to get the last; in case of a reversed read stop after the first
						last_error_pos = read_pos - 1;
					}
					++num_errors;
				}
			}
			break;
		case 'N':
		case 'D':
			// Deletions can be attributed to the base before or after the deletion, choose the one resulting in the lowest last_error_pos in the end (so shortest adapter required length)
			if( hasFlagRC(record) ){
				if( !num_errors ){
					last_error_pos = read_pos;
				}
			}
			else{
				last_error_pos = read_pos-1;
			}

			num_errors += cigar_element.count;
			ref_pos += cigar_element.count;
			break;
		case 'I':
			// Insertions can be attributed to the first or last base of the insertion, choose the one resulting in the lowest last_error_pos in the end (so shortest adapter required length)
			if( hasFlagRC(record) ){
				if( !num_errors ){
					last_error_pos = read_pos;
				}
			}
			else{
				last_error_pos = read_pos+cigar_element.count-1;
			}

			num_errors += cigar_element.count;
			read_pos += cigar_element.count;
			break;
		default:
			printErr << "Unknown cigar operation: " << cigar_element.operation << std::endl;
		}
	}

	// Reverse position for normal reads to get them counted from the end of the read
	if( !hasFlagRC(record) ){
		last_error_pos = length(record.seq) - 1 - last_error_pos;
	}

	return num_errors;
}

bool AdapterStats::VerifyInsert(const CharString &read1, const CharString &read2, intSeqShift pos1, intSeqShift pos2){
	if( 1 > pos1 || 1 > pos2){
		return true;
	}

	if( pos1 == pos2 ){
		uintReadLen matches(0), i2(0);
		for(uintReadLen i1=pos1; i1--; ){
			switch(read1[i1]){
			case 'A':
				if( 'T' == read2[i2] ){
					++matches;
				}
				break;
			case 'C':
				if( 'G' == read2[i2] ){
					++matches;
				}
				break;
			case 'G':
				if( 'C' == read2[i2] ){
					++matches;
				}
				break;
			case 'T':
				if( 'A' == read2[i2] ){
					++matches;
				}
				break;
			}
			++i2;
		}

		if(matches > pos1*95/100){
			return true;
		}
	}

	return false;
}

bool AdapterStats::AdaptersAmbigous(const seqan::DnaString &adaper1, const seqan::DnaString &adaper2, uintReadLen max_length){
	auto compare_until = max( static_cast<uintReadLen>(max(length(adaper1), length(adaper2))), max_length );
	uintReadLen num_diff = 0;
	for( auto pos=compare_until; pos--; ){
		if(adaper1[pos] != adaper2[pos]){
			++num_diff;
		}
	}
	return num_diff < 2+compare_until/10;
}

AdapterStats::AdapterStats(){
	// Clear adapter_overrun_bases_
	for( auto i = tmp_overrun_bases_.size(); i--; ){
		tmp_overrun_bases_.at(i) = 0;
	}
}

bool AdapterStats::LoadAdapters(const char *adapter_file, const char *adapter_matrix, uintQual phred_quality_offset, uintReadLen size_read_length){
	StringSet<CharString> ids;
	StringSet<DnaString> seqs;

	printInfo << "Loading adapters from: " << adapter_file << std::endl;
	printInfo << "Loading adapter combination matrix from: " << adapter_matrix << std::endl;

	try{
		SeqFileIn seq_file_in(adapter_file);
		readRecords(ids, seqs, seq_file_in);
	}
	catch(const Exception &e){
		printErr << "Could not load adapter sequences: " << e.what() << std::endl;
		return false;
	}

	vector<string> matrix;
	matrix.resize(length(seqs));

	ifstream matrix_file;
	matrix_file.open(adapter_matrix);

	uintAdapterId nline=0;
	while( length(seqs) > nline && matrix_file.good() ){
		getline(matrix_file, matrix.at(nline++));
	}

	matrix_file.close();

	// Check if matrix fits to adapter file
	bool error = false;
	if(length(seqs) != nline ){
		printErr << "Matrix has less rows than adapters loaded from file." << std::endl;
		error = true;
	}
	matrix_file.get(); // Test if file is at the end
	if( !matrix_file.eof() ){
		printErr << "Matrix has more rows than adapters loaded from file." << std::endl;
		error = true;
	}
	for( ; nline--; ){
		if( length(seqs) > matrix.at(nline).size() ){
			printErr << "Matrix row " << nline << " has less columns than adapters loaded from file." << std::endl;
			error = true;
		}
		else if( length(seqs) < matrix.at(nline).size() ){
			printErr << "Matrix row " << nline << " has more columns than adapters loaded from file." << std::endl;
			error = true;
		}
		for(auto nchar = matrix.at(nline).size(); nchar--; ){
			if( '0' != matrix.at(nline).at(nchar) && '1' != matrix.at(nline).at(nchar)){
				printErr << "Not allowed character '" << matrix.at(nline).at(nchar) << "' in matrix at line " << nline << ", position " << nchar << std::endl;
				error = true;
			}
		}
	}
	if(error){
		return false;
	}

	// Get adapters that are allowed for first or second
	array<vector<pair<DnaString, uintAdapterId>>, 2> adapter_list;
	for(uintTempSeq seg=2; seg--; ){
		adapter_list.at(seg).reserve(length(seqs));
	}

	array<uintReadLen, 2> size_adapter = {0,0};
	for( nline = length(seqs); nline--; ){
		for(auto nchar = length(seqs); nchar--; ){
			if( '1' == matrix.at(nline).at(nchar) ){
				// If one adapter pair has this adaptor as first, add it to first and continue with the next adapter
				adapter_list.at(0).push_back({seqs[nline],nline});
				SetToMax(size_adapter.at(0), length(seqs[nline]));
				break;
			}
		}
	}

	for(auto nchar = length(seqs); nchar--; ){
		for( nline = length(seqs); nline--; ){
			if( '1' == matrix.at(nline).at(nchar) ){
				adapter_list.at(1).push_back({seqs[nchar],nchar});
				SetToMax(size_adapter.at(1), length(seqs[nchar]));
				break;
			}
		}
	}

	for(uintTempSeq seg=2; seg--; ){
		// Also add reverse complement of adapters (as the direction is different depending on the sequencing machine Hiseq2000 vs. 4000)
		adapter_list.at(seg).reserve(adapter_list.at(seg).size()*2);
		for( uintAdapterId i=adapter_list.at(seg).size(); i--;){
			adapter_list.at(seg).push_back(adapter_list.at(seg).at(i));
			reverseComplement(adapter_list.at(seg).back().first);
			adapter_list.at(seg).back().second += length(ids);
		}

		// Sort adapters by sequence (necessary for SolveAmbiguities function later) and fill the class variables
		sort(adapter_list.at(seg));
		seqs_.at(seg).reserve(adapter_list.at(seg).size());
		names_.at(seg).reserve(adapter_list.at(seg).size());

		for( auto &adapter : adapter_list.at(seg) ){
			seqs_.at(seg).push_back(adapter.first);
			if(adapter.second < length(ids)){
				names_.at(seg).push_back((string(toCString(ids[adapter.second]))+"_f").c_str());
			}
			else{
				names_.at(seg).push_back((string(toCString(ids[adapter.second-length(ids)]))+"_r").c_str());
			}
		}
	}

	// Fill adapter_combinations_ with the valid combinations of both adapters
	combinations_.clear();
	combinations_.resize(adapter_list.at(0).size(), vector<bool>(adapter_list.at(1).size(), true) ); // Set all combinations by default to valid
	for( auto adapter1 = adapter_list.at(0).size(); adapter1--; ){
		for( auto adapter2 = adapter_list.at(1).size(); adapter2--; ){
			if('0' == matrix.at( adapter_list.at(0).at(adapter1).second%length(ids) ).at( adapter_list.at(1).at(adapter2).second%length(ids) )){
				// Combination not valid
				combinations_.at(adapter1).at(adapter2) = false;
			}
		}
	}

	// Resize adapter vectors
	tmp_counts_.resize(adapter_list.at(0).size());
	for( auto &dim1 : tmp_counts_ ){
		dim1.resize(adapter_list.at(1).size());
		for( auto &dim2 : dim1 ){
			dim2.resize(size_read_length);
			for( auto &dim3 : dim2 ){
				dim3.resize(size_read_length);
			}
		}
	}

	for(uintTempSeq seg=2; seg--; ){
		tmp_start_cut_.at(seg).resize(adapter_list.at(seg).size());
		for( auto &dim1 : tmp_start_cut_.at(seg) ){
			dim1.resize(size_adapter.at(seg));
		}
	}

	tmp_polya_tail_length_.resize(size_read_length);

	// Prepare skewer matrix for adapter identification
	skewer::cMatrix::InitParameters(skewer::TRIM_PE, 0.1, 0.03, phred_quality_offset, false); // 0.1 and 0.03 are the default values from the skewer software

	for( auto &adapter : seqs_.at(0) ){
		skewer::cMatrix::AddAdapter(skewer::cMatrix::firstAdapters, toCString(static_cast<CharString>(adapter)), length(adapter), skewer::TRIM_TAIL);
	}
	for( auto &adapter : seqs_.at(1) ){
		skewer::cMatrix::AddAdapter(skewer::cMatrix::secondAdapters, toCString(static_cast<CharString>(adapter)), length(adapter), skewer::TRIM_TAIL);
	}
	skewer::cMatrix::CalculateIndices(combinations_, adapter_list.at(0).size(), adapter_list.at(1).size());

	return true;
}

bool AdapterStats::Detect(uintReadLen &adapter_position_first, uintReadLen &adapter_position_second, const BamAlignmentRecord &record_first, const BamAlignmentRecord &record_second, const Reference &reference, bool properly_mapped){
	bool search_adapters(true);
	// Determine minimum adapter length
	int i_min_overlap;
	if( hasFlagUnmapped(record_first) || hasFlagUnmapped(record_second) ){
		uintReadLen last_error_pos(kMinimumAdapterLengthUnmapped-1); // Case of both reads unmapped
		if( (hasFlagUnmapped(record_first) || CountErrors(last_error_pos, record_first, reference)) && (hasFlagUnmapped(record_second) || CountErrors(last_error_pos, record_second, reference)) ){
			i_min_overlap = last_error_pos+1; // +1 because it is a position starting at 0 and we need a length
		}
		else{
			// In case one of the reads is a perfect match don't look for adapters
			search_adapters = false;
		}
	}
	else{
		uintReadLen last_error_pos1, last_error_pos2;
		if( CountErrors(last_error_pos1, record_first, reference) && CountErrors(last_error_pos2, record_second, reference) ){
			i_min_overlap = max(last_error_pos1, last_error_pos2)+1; // A perfect matching piece cannot be an adapter, so adapter must be at least reach last error in both reads
		}
		else{
			search_adapters = false;
		}
	}

	if(search_adapters){
		// Fill read1, read2 with first/second measured in sequencer(by fastq file) not by position in bam file and convert to measurement/fastq direction (not reference/bam direction)
		const CharString *qual1, *qual2;
		CharString read1, read2, reversed_qual;

		if( hasFlagLast(record_second) ){
			read1 = record_first.seq;
			if( hasFlagRC(record_first) ){
				reverseComplement(read1); // Reverses in place, so don't let it act on the record sequence itself
				reversed_qual = record_first.qual;
				reverse(reversed_qual);
				qual1 = &reversed_qual;
			}
			else{
				qual1 = &record_first.qual;
			}

			read2 = record_second.seq;
			if( hasFlagRC(record_second) ){
				reverseComplement(read2); // Reverses in place, so don't let it act on the record sequence itself
				reversed_qual = record_second.qual;
				reverse(reversed_qual);
				qual2 = &reversed_qual;
			}
			else{
				qual2 = &record_second.qual;
			}
		}
		else{
			read1 = record_second.seq;
			if( hasFlagRC(record_second) ){
				reverseComplement(read1); // Reverses in place, so don't let it act on the record sequence itself
				reversed_qual = record_second.qual;
				reverse(reversed_qual);
				qual1 = &reversed_qual;
			}
			else{
				qual1 = &record_second.qual;
			}

			read2 = record_first.seq;
			if( hasFlagRC(record_first) ){
				reverseComplement(read2); // Reverses in place, so don't let it act on the record sequence itself
				reversed_qual = record_first.qual;
				reverse(reversed_qual);
				qual2 = &reversed_qual;
			}
			else{
				qual2 = &record_first.qual;
			}
		}

		auto index1 = skewer::cMatrix::findAdapter(toCString(read1), length(read1), reinterpret_cast<unsigned char *>(toCString(*qual1)), length(*qual1), i_min_overlap);
		if( index1.bc-- ){ // Adapter has been detected for first read (bc is now index of adapter)
			auto index2 = skewer::cMatrix::findAdapter2(toCString(read2), length(read2), reinterpret_cast<unsigned char *>(toCString(*qual2)), length(*qual2), i_min_overlap);
			if( index2.bc-- && combinations_.at(index1.bc).at(index2.bc)){ // Adapter has been detected also for second read and it is a valid pair with the first
				if( 2 > max(index1.pos,index2.pos) - min(index1.pos,index2.pos) ){ // Allow a single one-base indel in the insert between the adapters
					if( properly_mapped || VerifyInsert(read1, read2, index1.pos, index2.pos) ){
						// Stats that need to be aware of cut length
						bool poly_a_tail = true;
						uintReadLen poly_a_tail_length = 0;
						for(uintReadLen pos=index1.pos + length(seqs_.at(0).at(index1.bc)); pos < length(read1); ++pos){
							if(poly_a_tail){
								if('A' == read1[pos]){
									++poly_a_tail_length;
								}
								else{
									poly_a_tail=false;
									++tmp_polya_tail_length_.at(poly_a_tail_length);
								}
							}

							if(!poly_a_tail){
								++tmp_overrun_bases_.at(Dna5(read1[pos]));
							}
						}

						poly_a_tail = true;
						poly_a_tail_length = 0;
						for(uintReadLen pos=index2.pos + length(seqs_.at(1).at(index2.bc)); pos < length(read2); ++pos){
							if(poly_a_tail){
								if('A' == read2[pos]){
									++poly_a_tail_length;
								}
								else{
									poly_a_tail=false;
									++tmp_polya_tail_length_.at(poly_a_tail_length);
								}
							}

							if(!poly_a_tail){
								++tmp_overrun_bases_.at(Dna5(read2[pos]));
							}
						}

						// If the read starts with an adapter note down if and how much it has been cut of at the beginning (position below 0) and set position to 0 for the further stats
						if(1>index1.pos){
							++tmp_start_cut_.at(0).at(index1.bc).at(-index1.pos);
							index1.pos = 0;
						}
						if(1>index2.pos){
							++tmp_start_cut_.at(1).at(index2.bc).at(-index2.pos);
							index2.pos = 0;
						}

						// Set return values
						if( hasFlagLast(record_first) ){
							adapter_position_first = index2.pos;
							adapter_position_second = index1.pos;
						}
						else{
							adapter_position_first = index1.pos;
							adapter_position_second = index2.pos;
						}

						// Fill in stats
						++tmp_counts_.at(index1.bc).at(index2.bc).at(length(read1)-index1.pos).at(length(read2)-index2.pos);

						return true;
					}
				}
			}
		}
	}

	// If no adapter where found fill in read length
	adapter_position_first = length(record_first.seq);
	adapter_position_second = length(record_second.seq);

	return false;
}

void AdapterStats::Finalize(){
	// Copy to final vectors
	counts_.resize(tmp_counts_.size());
	for( auto i = tmp_counts_.size(); i--; ){
		counts_.at(i).resize(tmp_counts_.at(i).size());
		for( auto j = tmp_counts_.at(i).size(); j--; ){
			counts_.at(i).at(j).Acquire(tmp_counts_.at(i).at(j));
			counts_.at(i).at(j).Shrink();
		}
	}

	for(uintTempSeq seg=2; seg--; ){
		start_cut_.at(seg).resize(tmp_start_cut_.at(seg).size());
		for( auto i = tmp_start_cut_.at(seg).size(); i--; ){
			start_cut_.at(seg).at(i).Acquire(tmp_start_cut_.at(seg).at(i));
		}
	}

	polya_tail_length_.Acquire(tmp_polya_tail_length_);

	for( auto i = tmp_overrun_bases_.size(); i--; ){
		overrun_bases_.at(i) = tmp_overrun_bases_.at(i);
	}
}

void AdapterStats::SumCounts(){
	for(uintTempSeq seg=2; seg--; ){
		count_sum_.at(seg).clear();
		count_sum_.at(seg).resize(start_cut_.at(seg).size(), 0);
	}

	// Sum adapters starting from length where adapters are unambiguous (sorted by content: so simply first position that differs from adapter before and after)
	uintFragCount sum;
	uintReadLen first_diff_before_a1(0), first_diff_after_a1, first_diff_before_a2, first_diff_after_a2;
	for( auto a1=counts_.size(); a1--; ){
		first_diff_after_a1 = 0;
		if(a1){ // Not last in loop
			while(first_diff_after_a1 < min(length(seqs_.at(0).at(a1)), length(seqs_.at(0).at(a1-1))) && seqs_.at(0).at(a1)[first_diff_after_a1] == seqs_.at(0).at(a1-1)[first_diff_after_a1] ){
				++first_diff_after_a1;
			}
		}

		first_diff_before_a2 = 0;
		for( auto a2=counts_.at(0).size(); a2--; ){ // All vectors i are the same length so we can simply take the length of 0
			first_diff_after_a2 = 0;
			if(a2){ // Not last in loop
				while(first_diff_after_a2 < min(length(seqs_.at(1).at(a2)), length(seqs_.at(1).at(a2-1))) && seqs_.at(1).at(a2)[first_diff_after_a2] == seqs_.at(1).at(a2-1)[first_diff_after_a2] ){
					++first_diff_after_a2;
				}
			}

			sum = 0;
			for( auto pos1 = max(max(first_diff_before_a1,first_diff_after_a1), static_cast<uintReadLen>(counts_.at(a1).at(a2).from())); pos1 < counts_.at(a1).at(a2).to(); ++pos1){
				for( auto pos2 = max(max(first_diff_before_a2,first_diff_after_a2), static_cast<uintReadLen>(counts_.at(a1).at(a2).at(pos1).from())); pos2 < counts_.at(a1).at(a2).at(pos1).to(); ++pos2){
					sum += counts_.at(a1).at(a2).at(pos1).at(pos2);
				}
			}

			count_sum_.at(0).at(a1) += sum;
			count_sum_.at(1).at(a2) += sum;

			first_diff_before_a2 = first_diff_after_a2;
		}

		first_diff_before_a1 = first_diff_after_a1;
	}
}

void AdapterStats::Shrink(){
	for( auto &dim1 : counts_){
		for( auto &adapter_pair : dim1){
			ShrinkVect(adapter_pair);
		}
	}
	for( int i=2; i--; ){
		for( auto &adapter : start_cut_.at(i)){
			adapter.Shrink();
		}
	}
	polya_tail_length_.Shrink();
}

void AdapterStats::PrepareSimulation(){
	for(uintTempSeq seg=2; seg--; ){
		significant_count_.at(seg).clear();
		significant_count_.at(seg).resize(count_sum_.at(seg).size(), 0);

		uintFragCount threshold = ceil(*max_element(count_sum_.at(seg).begin(), count_sum_.at(seg).end()) * kMinFractionOfMaximumForSimulation);

		for(auto i=significant_count_.at(seg).size(); i--; ){
			if(count_sum_.at(seg).at(i) < threshold){
				significant_count_.at(seg).at(i) = 0;
			}
			else{
				significant_count_.at(seg).at(i) = count_sum_.at(seg).at(i);
			}
		}
	}
}
