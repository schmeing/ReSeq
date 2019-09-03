#include "Reference.h"
using reseq::Reference;

#include <algorithm>
using std::min;
using std::max;
using std::sort;
//include <array>
using std::array;
#include <cmath>
using std::log;
using std::pow;
using std::round;
#include <iostream>
using std::cout;
#include <limits>
using std::numeric_limits;
#include <random>
using std::mt19937_64;
using std::uniform_int_distribution;
//include <set>
using std::set;
#include <string>
using std::string;
//include "utilities.h"
using std::pair;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/seq_io.h>
using seqan::CharString;
using seqan::clear;
using seqan::DnaString;
using seqan::Dna5;
using seqan::Dna5String;
using seqan::Exception;
using seqan::FunctorComplement;
using seqan::IupacString;
using seqan::length;
using seqan::ModComplementDna5;
using seqan::ModifiedString;
using seqan::Prefix;
using seqan::prefix;
using seqan::readRecord;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::Size;
//include <seqan/vcf_io.h>
using seqan::VcfFileIn;
using seqan::VcfHeader;
using seqan::VcfRecord;

#include "CMakeConfig.h"
using reseq::utilities::Percent;

bool Reference::CheckVcf() const{
	auto contigs = contigNames(context(vcf_file_));
	uint32_t errors = 0;

	if(length(contigs) != NumberSequences()){
		printErr << "Number of contigs does not match between reference(" << NumberSequences() << ") and variant(" << length(contigs) << ") file." << std::endl;
		++errors;
	}
	for(uint32_t con=0; con < min(length(contigs), NumberSequences()); ++con){
		if(errors < 20 && contigs[con] != ReferenceIdFirstPart(con)){
			printErr << "Contigs at position " << con << " do not match between reference(" << ReferenceId(con) << ") and variant(" << contigs[con] << ") file." << std::endl;
			++errors;
		}
	}

	return !errors;
}

inline bool Reference::ReadFirstVcfRecord(){
	try{
		// Read the first record, so that cur_vcf_record_ is valid
		readRecord(cur_vcf_record_, vcf_file_);
	}
	catch (const Exception& e){
		printErr << "Could not read first vcf record: " << e.what() << std::endl;
		return false;
	}

	return true;
}

bool Reference::ReadVariants(uint32_t end_ref_seq_id, bool positions_only){
	// Read the file record by record until variant is for sequence end_ref_seq_id or later
	uint32_t errors=0;

	uint32_t start_pos, end_pos(0);

	uint16_t pos;
	uint16_t cur_allele, chosen_var;
	vector<uint16_t> allele;
	vector<uint64_t> gt_has_var;
	vector<uint32_t> alt_start_pos;

	if(!positions_only){
		allele.resize(num_alleles_);
		gt_has_var.resize(num_alleles_);
	}

	uint32_t old_ref_id = numeric_limits<uint32_t>::max();

	try{
		while(VcfRecord::INVALID_REFID != cur_vcf_record_.rID && cur_vcf_record_.rID < end_ref_seq_id){
			// Parse record
			if(cur_vcf_record_.rID >= NumberSequences()){
				printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " does not belong to an existing reference sequence." << std::endl;
				if(++errors >= 20){
					break;
				}
				else{
					continue;
				}
			}
			if(cur_vcf_record_.beginPos >= SequenceLength(cur_vcf_record_.rID)){
				printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " starts after the end of the reference sequence." << std::endl;
				if(++errors >= 20){
					break;
				}
				else{
					continue;
				}
			}

			start_pos = cur_vcf_record_.beginPos;
			if(old_ref_id == cur_vcf_record_.rID){
				if(start_pos < end_pos){
					printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " overlaps with a previous variant." << std::endl;
					if(++errors >= 20){
						break;
					}
					else{
						continue;
					}
				}
			}
			else{
				old_ref_id = cur_vcf_record_.rID;
			}

			// Skip variants that contain N's in the reference for full variants, positions only takes them as well
			uint16_t n_count = 0;
			if(!positions_only){
				for(auto var_pos=length(cur_vcf_record_.ref); var_pos--; ){
					if('N' == cur_vcf_record_.ref[var_pos]){
						++n_count;
					}
				}
			}

			end_pos = start_pos + length(cur_vcf_record_.ref);
			if( 0 == n_count && cur_vcf_record_.ref != infix(ReferenceSequence(cur_vcf_record_.rID), start_pos, end_pos) ){
				printErr << "The specified reference in vcf file '" << cur_vcf_record_.ref << "' is not identical with the specified reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << ": '" << infix(ReferenceSequence(cur_vcf_record_.rID), start_pos, end_pos ) << "'." << std::endl;
				if(++errors >= 20){
					break;
				}
				else{
					continue;
				}
			}

			if(positions_only){
				for(auto pos=start_pos; pos<end_pos; ++pos){
					variant_positions_.at(cur_vcf_record_.rID).push_back(pos);
				}
			}
			else if(0 == n_count){
				// Get genotypes
				cur_allele = 0;
				for( auto &genotype : cur_vcf_record_.genotypeInfos ){
					if(cur_allele >= num_alleles_){
						printErr << "Found to many alleles in genotype definition '";
						for( auto &genotype : cur_vcf_record_.genotypeInfos ){
							std::cout << ' ' << genotype;
						}
						std::cout << "'" << std::endl;
					}

					pos=0;
					chosen_var=0;
					while(pos < length(genotype) && ':' != genotype[pos]){
						if('|' == genotype[pos] || '/' == genotype[pos]){
							allele.at(cur_allele++) = chosen_var;
							chosen_var = 0;
						}
						else if('0' <= genotype[pos] && '9' >= genotype[pos]){
							chosen_var *= 10;
							chosen_var += genotype[pos] - 48;
						}
						else{
							printErr << "Unallowed character '" << genotype[pos] << "' in genotype definition '" << genotype << "'" << std::endl;
							if(++errors >= 20){
								break;
							}
							else{
								continue;
							}
						}
						++pos;
					}
					allele.at(cur_allele++) = chosen_var;
				}

				if(cur_allele < num_alleles_){
					printErr << "Could not find enough alleles in genotype definition '";
					for( auto &genotype : cur_vcf_record_.genotypeInfos ){
						std::cout << ' ' << genotype;
					}
					std::cout << "'" << std::endl;
					if(++errors >= 20){
						break;
					}
					else{
						continue;
					}
				}

				// Find variations for the genotypes
				gt_has_var.clear();
				alt_start_pos.clear();
				alt_start_pos.push_back(0);

				chosen_var = 1; // 0 is the reference sequence
				// Split alternative sequences
				for(pos = 0; pos < length(cur_vcf_record_.alt); ++pos){
					if(',' == cur_vcf_record_.alt[pos]){
						alt_start_pos.push_back(pos+1); // Position after the ',' is the start of the next alternative

						// Check which genotypes have this variant
						gt_has_var.push_back(0);
						for(cur_allele=num_alleles_; cur_allele--; ){ // Has to be backwards, because allele 0 is in the rightmost bit
							gt_has_var.back() = gt_has_var.back() << 1;
							if(allele.at(cur_allele) == chosen_var){
								++gt_has_var.back();
							}
						}
						++chosen_var;
					}
				}
				alt_start_pos.push_back(pos+1); // To be in line with the one after the ',' we use here one after the end (which itself is one after the last character)

				// Check which genotypes have the final variant
				gt_has_var.push_back(0);
				for(cur_allele=num_alleles_; cur_allele--; ){ // Has to be backwards, because allele 0 is in the rightmost bit
					gt_has_var.back() = gt_has_var.back() << 1;
					if(allele.at(cur_allele) == chosen_var){
						++gt_has_var.back();
					}
					else if(allele.at(cur_allele) > chosen_var){
						printErr << "Variant number " << allele.at(cur_allele) << " does not exist for sequence id " << cur_vcf_record_.rID << " and position " << cur_vcf_record_.beginPos << std::endl;
						if(++errors >= 20){
							break;
						}
						else{
							continue;
						}
					}
				}
				++chosen_var;

				// Enter variants into vector (split into single reference positions)
				for(pos = 0; pos < length(cur_vcf_record_.ref); ++pos){
					for(uint16_t n_alt = 0; n_alt < gt_has_var.size(); ++n_alt){
						// Check if any genotype has this variant
						if(gt_has_var.at(n_alt)){
							// Enter record into variations vector
							if(pos+1 == length(cur_vcf_record_.ref) && pos+1 < alt_start_pos.at(n_alt+1)-1 - alt_start_pos.at(n_alt)){
								// Insertion
								InsertVariant(cur_vcf_record_.rID, start_pos+pos, infix(cur_vcf_record_.alt, alt_start_pos.at(n_alt)+pos, alt_start_pos.at(n_alt+1)-1), gt_has_var.at(n_alt));
							}
							else if(pos < alt_start_pos.at(n_alt+1)-1 - alt_start_pos.at(n_alt)){
								// Base Mutation
								if(cur_vcf_record_.ref[pos] != cur_vcf_record_.alt[alt_start_pos.at(n_alt)+pos]){ // If variant is identical to reference, we don't need to add it
									InsertVariant(cur_vcf_record_.rID, start_pos+pos, DnaString(cur_vcf_record_.alt[alt_start_pos.at(n_alt)+pos]), gt_has_var.at(n_alt));
								}
							}
							else{
								// Deletion
								InsertVariant(cur_vcf_record_.rID, start_pos+pos, DnaString(""), gt_has_var.at(n_alt));
							}
						}
					}
				}
			}

			// Read next record
			if(atEnd(vcf_file_)){
				cur_vcf_record_.rID = VcfRecord::INVALID_REFID;
				read_variation_for_num_sequences_ = NumberSequences();
			}
			else{
				readRecord(cur_vcf_record_, vcf_file_);

				if( cur_vcf_record_.rID < read_variation_for_num_sequences_ ){
					printErr << "Variant file is not properly position sorted. Found sequence id " << cur_vcf_record_.rID << " after id " << read_variation_for_num_sequences_ << std::endl;
					if(++errors >= 20){
						break;
					}
					else{
						continue;
					}
				}
				else if(cur_vcf_record_.rID == read_variation_for_num_sequences_){
					if(cur_vcf_record_.beginPos < start_pos){
						printErr << "Variant file is not properly position sorted. Found in sequence id " << cur_vcf_record_.rID << " position " << cur_vcf_record_.beginPos << " after position " << start_pos << std::endl;
						if(++errors >= 20){
							break;
						}
						else{
							continue;
						}
					}
				}
				else{
					read_variation_for_num_sequences_ = cur_vcf_record_.rID;
				}
			}
		}
	}
	catch (const Exception& e){
		printErr << "Could not read vcf record: " << e.what() << std::endl;
		++errors;
	}

	if(errors){
		cur_vcf_record_.rID = VcfRecord::INVALID_REFID;
	}

	return !errors;
}

void Reference::CloseVcfFile(){
	if(VcfRecord::INVALID_REFID != cur_vcf_record_.rID){
		printErr << "Variant file is not sorted by position. The following variation record and all after it have been ignored during simulation:\n";
		std::cout << cur_vcf_record_.rID << '\t' << cur_vcf_record_.beginPos+1 << '\t' << cur_vcf_record_.id << '\t' << cur_vcf_record_.ref << '\t' << cur_vcf_record_.alt << '\t' << cur_vcf_record_.qual << '\t' << cur_vcf_record_.filter << '\t' << cur_vcf_record_.info << '\t' << cur_vcf_record_.format;
		for(auto &genotype : cur_vcf_record_.genotypeInfos){
			std::cout << '\t' << genotype;
		}
		std::cout << std::endl;
	}
	clear(contigNames(context(vcf_file_)));
	close(vcf_file_);
}

uint32_t Reference::ForwardSurroundingBlock( Size<decltype(reference_sequences_)>::Type ref_seq_id, uint64_t pos ) const{
	const Dna5String &ref_seq(reference_sequences_[ref_seq_id]);

	uint32_t return_value(0);
	if( pos+surrounding_range_ <= length(ref_seq) ){
		for( auto ref_pos = pos; ref_pos < pos+surrounding_range_; ++ref_pos ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
	}
	else{
		// Surrounding outside of sequence at end
		for( auto ref_pos = pos; ref_pos < length(ref_seq); ++ref_pos ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
		// Circle around the reference sequence and continue at the beginning
		for( auto ref_pos = 0; ref_pos < (pos+surrounding_range_)%length(ref_seq); ++ref_pos ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
	}

	return return_value;
}

uint32_t Reference::ReverseSurroundingBlock( Size<decltype(reference_sequences_)>::Type ref_seq_id, uint64_t pos ) const{
	const ModifiedString<const Dna5String, ModComplementDna5> ref_seq(reference_sequences_[ref_seq_id]);

	uint32_t return_value(0);
	if( pos+surrounding_range_ <= length(ref_seq) ){

		for( auto ref_pos = pos+surrounding_range_; ref_pos-- > pos; ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
	}
	else{
		// Surrounding outside of sequence at end
		// Circle around the reference sequence and start at the beginning
		for( auto ref_pos = pos+surrounding_range_-length(ref_seq); ref_pos--; ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
		// Go back to the end
		for( auto ref_pos = length(ref_seq); ref_pos-- > pos; ){
			return_value = (return_value << 2) + static_cast<uint16_t>(ref_seq[ref_pos]);
		}
	}

	return return_value;
}


inline void Reference::AddFragmentSite(
		std::vector<FragmentSite> &sites,
		uint32_t fragment_length,
		uint32_t gc,
		uint32_t n_count,
		const array<int32_t, num_surrounding_blocks_> &start_sur,
		const array<int32_t, num_surrounding_blocks_> &end_sur
		) const{
	bool valid = n_count <= max_n_in_fragment_site_;
	for(auto block = start_sur.size(); block--; ){
		if( 0 > start_sur.at(block) || 0 > end_sur.at(block) ){
			valid = false;
		}
	}
	if(valid){
		uint16_t gc_perc = Percent(gc, fragment_length-n_count);
		sites.emplace_back(gc_perc, start_sur, end_sur);
	}
	else{
		sites.emplace_back(); // Just a placeholder, so that the indexes in sites match starting positions(with the shift due to min_dist_to_ref_seq_ends_), will be removed after counts have been entered
	}
}

Reference::Reference():
	min_dist_to_ref_seq_ends_(50),
	max_n_in_fragment_site_(50),
	num_alleles_(1) // In case we don't load a variant file, we have a single allele
	{
}

const Prefix<const CharString>::Type Reference::ReferenceIdFirstPart( Size<decltype(reference_ids_)>::Type n ) const{
	Size<decltype(reference_ids_)>::Type pos=0; // declare before the for loop to be able to access it afterwards
	for( ; pos < length(reference_ids_[n]) && ' ' != reference_ids_[n][pos]; ++pos ); // loop until pos is at the first ' ' or at the end of the string
	return prefix(reference_ids_[n], pos);
}

void Reference::ReferenceSequence(
		DnaString &insert_string,
		Size<decltype(reference_sequences_)>::Type seq_id,
		Size<Dna5String>::Type start_pos,
		Size<Dna5String>::Type min_length,
		bool reversed,
		Size<Dna5String>::Type seq_size ) const{
	if( reversed ){
		insert_string = infix(reference_sequences_[seq_id], start_pos-min_length, start_pos);
		reverseComplement(insert_string);
	}
	else{
		insert_string = infix(reference_sequences_[seq_id], start_pos, start_pos+min_length);
	}

	if( seq_size > length(insert_string) ){
		resize(insert_string, seq_size);
	}
}

void Reference::ReferenceSequence(
		DnaString &insert_string,
		Size<decltype(reference_sequences_)>::Type seq_id,
		Size<Dna5String>::Type start_pos,
		Size<Dna5String>::Type min_length,
		bool reversed,
		const vector<Variant> &variants,
		pair<int32_t, uint16_t> first_variant, // {id, posCurrentlyAt}
		uint16_t allele,
		Size<Dna5String>::Type seq_size ) const{
	if( reversed ){
		insert_string = "";
		auto cur_start = start_pos;
		auto cur_var = first_variant.first;
		if(first_variant.second){
			insert_string += ReverseComplementor(prefix(variants.at(cur_var).var_seq_, first_variant.second));
			--cur_var;
			--cur_start;
		}
		for( ; cur_var >= 0 && length(insert_string) < min_length; --cur_var ){
			if(variants.at(cur_var).InAllele(allele)){
				if(cur_start-variants.at(cur_var).position_ > min_length-length(insert_string)){
					// Variant after return sequence
					insert_string += ReverseComplementor(infix(reference_sequences_[seq_id], cur_start+length(insert_string)-min_length, cur_start));
				}
				else{
					// Variant inside return sequence
					insert_string += ReverseComplementor(infix(reference_sequences_[seq_id], variants.at(cur_var).position_+1, cur_start));
					insert_string += ReverseComplementor(variants.at(cur_var).var_seq_);
					cur_start = variants.at(cur_var).position_;
				}
			}
		}
		if(cur_var == -1 && length(insert_string) < min_length){
			// Fill rest after the last variant
			insert_string += ReverseComplementor(infix(reference_sequences_[seq_id], cur_start+length(insert_string)-min_length, cur_start));
		}
	}
	else{
		insert_string = "";
		auto cur_start = start_pos;
		auto cur_var = first_variant.first;
		if(first_variant.second){
			insert_string += suffix(variants.at(cur_var).var_seq_, first_variant.second);
			++cur_var;
			++cur_start;
		}
		for( ; cur_var < variants.size() && length(insert_string) < min_length; ++cur_var ){
			if(variants.at(cur_var).InAllele(allele)){
				if(variants.at(cur_var).position_-cur_start >= min_length-length(insert_string)){
					// Variant after return sequence
					insert_string += infix(reference_sequences_[seq_id], cur_start, cur_start+min_length-length(insert_string));
				}
				else{
					// Variant inside return sequence
					insert_string += infix(reference_sequences_[seq_id], cur_start, variants.at(cur_var).position_);
					insert_string += variants.at(cur_var).var_seq_;
					cur_start = variants.at(cur_var).position_+1;
				}
			}
		}
		if(cur_var == variants.size() && length(insert_string) < min_length){
			// Fill rest after the last variant
			insert_string += infix(reference_sequences_[seq_id], cur_start, cur_start+min_length-length(insert_string));
		}
	}

	if( seq_size > length(insert_string) ){
		resize(insert_string, seq_size);
	}
}

uint64_t Reference::TotalSize() const{
	uint64_t sum = 0;
	for(const auto &seq : reference_sequences_){
		sum += length(seq);
	}
	return sum;
}

void Reference::RefSeqsSortedByNxx(vector<pair<double, uint32_t>> &ref_seqs) const{
	double tot_size = TotalSize();
	tot_size /= 100; // Nxx in percent

	ref_seqs.reserve( NumberSequences() );
	for(uint32_t id=0; id < NumberSequences(); ++id){
		ref_seqs.emplace_back( SequenceLength(id)/tot_size , id );
	}

	sort(ref_seqs.begin(), ref_seqs.end(), std::greater<pair<double, uint32_t>>());

	double sum = 0.0;
	for(uint32_t id=0; id < NumberSequences(); ++id){
		sum += ref_seqs.at(id).first;
		ref_seqs.at(id).first = sum;
	}
}

uint32_t Reference::RefSeqsInNxx(vector<bool> &ref_seqs, double nxx_ref_seqs) const{
	// Get nxx for all reference sequences
	vector<pair<double, uint32_t>> ref_seq_nxx;
	RefSeqsSortedByNxx(ref_seq_nxx);

	// Count how many sequences are needed to get nxx_ref_seqs taking only the largest ones
	uint32_t needed_ref_seqs = 0;
	ref_seqs.resize(NumberSequences(), false);
	while(needed_ref_seqs+1 < ref_seq_nxx.size() && nxx_ref_seqs > ref_seq_nxx.at(needed_ref_seqs).first){
		ref_seqs.at(needed_ref_seqs++) = true;
	}

	ref_seqs.at(needed_ref_seqs++) = true; // Get last one and convert index to size

	return needed_ref_seqs;
}

uint32_t Reference::NumRefSeqBinsInNxx(const std::vector<bool> &ref_seqs, uint64_t max_ref_seq_bin_length) const{
	uint32_t num_bins = 0;
	for(uint32_t ref_id=0; ref_id<ref_seqs.size(); ++ref_id){
		if(ref_seqs.at(ref_id)){
			num_bins += SequenceLength(ref_id) / max_ref_seq_bin_length + 1;
		}
	}
	return num_bins;
}

double Reference::SumBias(
		double &max_bias,
		Size<decltype(reference_sequences_)>::Type ref_seq_id,
		Size<Dna5String>::Type fragment_length,
		double general_bias,
		const Vect<double> &gc_bias,
		const array<vector<double>, num_surrounding_blocks_> &sur_bias) const{
	// Ignores min_dist_to_ref_seq_ends_ and does not account for N's as it is used for simulation only
	double tot(0.0);

	uint32_t n_count(0);
	uint64_t gc = GCContentAbsolut(n_count, ref_seq_id, 0, fragment_length);
	array<uint32_t, num_surrounding_blocks_> start_sur, end_sur;
	array<double, num_surrounding_blocks_> start_bias, end_bias;
	ForwardSurrounding(start_sur, ref_seq_id, 0);
	ReverseSurrounding(end_sur, ref_seq_id, fragment_length-1);

	for(auto block = start_bias.size(); block--; ){
		start_bias.at(block) = sur_bias.at(block).at(start_sur.at(block));
		end_bias.at(block) = sur_bias.at(block).at(end_sur.at(block));
	}

	double bias = Bias(general_bias, gc_bias[Percent(gc, fragment_length)], start_bias, end_bias);
	if(bias > max_bias){
		max_bias = bias;
	}
	tot += bias;

	const Dna5String &ref_seq(reference_sequences_[ref_seq_id]);
	for( uint32_t start_pos=0; start_pos < length(ref_seq)-fragment_length; ){
		UpdateGC( gc, n_count, ref_seq, start_pos, start_pos+fragment_length );
		UpdateReverseSurrounding( end_sur, ref_seq, start_pos+fragment_length );
		UpdateForwardSurrounding( start_sur, ref_seq, ++start_pos );

		for(auto block = start_bias.size(); block--; ){
			start_bias.at(block) = sur_bias.at(block).at(start_sur.at(block));
			end_bias.at(block) = sur_bias.at(block).at(end_sur.at(block));
		}

		double bias = Bias(general_bias, gc_bias[Percent(gc, fragment_length)], start_bias, end_bias);
		if(bias > max_bias){
			max_bias = bias;
		}
		tot += bias;
	}

	return tot;
}

double Reference::SumBias(
		std::vector<utilities::VectorAtomic<uint64_t>> &gc_sites,
		Size<decltype(reference_sequences_)>::Type ref_seq_id,
		Size<Dna5String>::Type fragment_length,
		double general_bias,
		const Vect<double> &gc_bias,
		const array<vector<double>, num_surrounding_blocks_> &sur_bias) const{
	// Uses min_dist_to_ref_seq_ends_ as it is used for stats creation
	double tot(0.0);

	uint32_t n_count(0);
	uint64_t gc = GCContentAbsolut(n_count, ref_seq_id, min_dist_to_ref_seq_ends_, min_dist_to_ref_seq_ends_+fragment_length);
	array<int32_t, num_surrounding_blocks_> start_sur, end_sur;
	array<double, num_surrounding_blocks_> start_bias, end_bias;
	ForwardSurroundingWithN(start_sur, ref_seq_id, min_dist_to_ref_seq_ends_);
	ReverseSurroundingWithN(end_sur, ref_seq_id, min_dist_to_ref_seq_ends_+fragment_length-1);

	bool valid=true;
	for(auto block = start_bias.size(); block--; ){
		if( 0 > start_sur.at(block) ){
			valid = false;
		}
		else{
			start_bias.at(block) = sur_bias.at(block).at(start_sur.at(block));
		}
		if( 0 > end_sur.at(block) ){
			valid = false;
		}
		else{
			end_bias.at(block) = sur_bias.at(block).at(end_sur.at(block));
		}
	}

	if(valid){
		tot += Bias(general_bias, gc_bias[Percent(gc, fragment_length)], start_bias, end_bias);
		++gc_sites.at(Percent(gc, fragment_length));
	}

	const Dna5String &ref_seq(reference_sequences_[ref_seq_id]);
	for( uint32_t start_pos=min_dist_to_ref_seq_ends_; start_pos < seqan::length(ref_seq)-fragment_length-min_dist_to_ref_seq_ends_; ){
		UpdateGC( gc, n_count, ref_seq, start_pos, start_pos+fragment_length );
		UpdateReverseSurroundingWithN( end_sur, ref_seq, start_pos+fragment_length, ref_seq_id );
		UpdateForwardSurroundingWithN( start_sur, ref_seq, ++start_pos, ref_seq_id );

		valid=true;
		for(auto block = start_bias.size(); block--; ){
			if( 0 > start_sur.at(block) ){
				valid = false;
			}
			else{
				start_bias.at(block) = sur_bias.at(block).at(start_sur.at(block));
			}
			if( 0 > end_sur.at(block) ){
				valid = false;
			}
			else{
				end_bias.at(block) = sur_bias.at(block).at(end_sur.at(block));
			}
		}
		if(valid){
			// If 0==fragment_length-n_count surrounding would be invalid, so no check needed
			tot += Bias(general_bias, gc_bias[Percent(gc, fragment_length)], start_bias, end_bias);
			++gc_sites.at(Percent(gc, fragment_length));
		}
	}

	return tot;
}

void Reference::GetFragmentSites( vector<FragmentSite> &sites, Size<decltype(reference_sequences_)>::Type ref_seq_id, Size<Dna5String>::Type fragment_length, uint32_t start, uint32_t end ) const{
	// Uses min_dist_to_ref_seq_ends_ as it is used for stats creation
	sites.clear();

	const Dna5String &ref_seq(reference_sequences_[ref_seq_id]);
	auto start_pos = max(static_cast<uint32_t>(min_dist_to_ref_seq_ends_), start);
	auto end_pos = min(length(ref_seq)-fragment_length-min_dist_to_ref_seq_ends_, static_cast<uint64_t>(end));

	uint32_t n_count(0);
	uint64_t gc = GCContentAbsolut(n_count, ref_seq_id, start_pos, start_pos+fragment_length);
	array<int32_t, num_surrounding_blocks_> start_sur, end_sur;
	ForwardSurroundingWithN(start_sur, ref_seq_id, start_pos);
	ReverseSurroundingWithN(end_sur, ref_seq_id, start_pos+fragment_length-1);

	AddFragmentSite( sites, fragment_length, gc, n_count, start_sur, end_sur );

	for( uint32_t frag_start=start_pos; frag_start < end_pos; ){
		UpdateGC( gc, n_count, ref_seq, frag_start, frag_start+fragment_length );
		UpdateReverseSurroundingWithN( end_sur, ref_seq, frag_start+fragment_length, ref_seq_id );
		UpdateForwardSurroundingWithN( start_sur, ref_seq, ++frag_start, ref_seq_id );

		AddFragmentSite( sites, fragment_length, gc, n_count, start_sur, end_sur );
	}
}

bool Reference::ReadFasta(const char *fasta_file){
	bool success = true;

	SeqFileIn ref;
	seqan::StringSet<seqan::IupacString> tmp_ref_seqs;
	clear(reference_ids_);
	if( !open(ref, fasta_file) ){
		printErr << "Could not open " << fasta_file << " for reading.\n";
		success = false;
	}
	else if(atEnd(ref)){
		printErr << fasta_file << " does not contain any reference sequences.\n";
		success = false;
	}
	else{
		try{
			readRecords(reference_ids_, tmp_ref_seqs, ref);
			printInfo << "Read in " << length(reference_ids_) << " reference sequences.\n";
		}
		catch(const Exception &e){
			printErr << "Could not read record in " << fasta_file << ": " << e.what() << '\n';
			success = false;
		}

	}

	reference_sequences_ = tmp_ref_seqs;

	return success;
}

void Reference::ReplaceN( uint64_t seed ){
	mt19937_64 rgen;
	rgen.seed(seed);
	uniform_int_distribution<> rdis(0, 3);

	for( auto &seq : reference_sequences_){
		for( auto &base : seq ){
			if(base > 3){
				base = rdis(rgen);
			}
		}
	}
}

bool Reference::hasN() const{
	for( auto &seq : reference_sequences_){
		for( auto &base : seq ){
			if(base > 3){
				return true;
			}
		}
	}

	return false;
}

bool Reference::WriteFasta(const char *fasta_file) const{
	bool success = true;

	SeqFileOut ref_out;
	if( !open(ref_out, fasta_file) ){
		printErr << "Could not open " << fasta_file << " for writing.\n";
		success = false;
	}
	else{
		try{
			writeRecords(ref_out, reference_ids_, reference_sequences_);
			printInfo << "Wrote " << length(reference_ids_) << " reference sequences.\n";
		}
		catch(const Exception &e){
			printErr << "Could not write generated sequences to " << fasta_file << ": " << e.what() << '\n';
			success = false;
		}
	}

	return success;
}

bool Reference::PrepareVariantFile(const string &var_file){
	bool success = true;

	try{
		// Open input file
		if( !open(vcf_file_, var_file.c_str()) ){
			printErr << "Could not open vcf file '" << var_file << "'." << std::endl;
			success = false;
		}
		else{
			// Read header
			VcfHeader header;
			readHeader(header, vcf_file_);

			// Check if vcf fits to reference
			if( !CheckVcf() ){
				success = false;
			}
			else{
				if(atEnd(vcf_file_)){
					printErr << "Vcf file '" << var_file << "' has no records." << std::endl;
					success = false;
				}
			}
		}
	}
	catch (const Exception& e){
		printErr << "Could not prepare vcf file '" << var_file << "' for record readin:" << e.what() << std::endl;
		success = false;
	}

	read_variation_for_num_sequences_ = 0;
	cleared_variation_for_num_sequences_ = 0;

	return success;
}

bool Reference::ReadFirstVariants(){
	if( ReadFirstVcfRecord() ){
		// Set number of alleles to correct value
		num_alleles_ = 0;
		for( auto &genotype : cur_vcf_record_.genotypeInfos ){
			uint16_t pos=0;
			while(pos < length(genotype) && ':' != genotype[pos]){
				if('|' == genotype[pos] || '/' == genotype[pos]){
					++num_alleles_; // Additional alleles in a population
				}
				++pos;
			}
			++num_alleles_; // Every population has at least one allele defined
		}

		if(1 == length(cur_vcf_record_.genotypeInfos)){
			if( 1 == num_alleles_ ){
				printInfo << "Reading variants for a single population with a single allele." << std::endl;
			}
			else{
				printInfo << "Reading variants for a single population with " << num_alleles_ << " alleles." << std::endl;
			}
		}
		else{
			printInfo << "Reading variants for " << length(cur_vcf_record_.genotypeInfos) << " populations with a total of " << num_alleles_ << " alleles." << std::endl;
		}

		variants_.resize(NumberSequences());
		return ReadVariants(2);
	}

	return false;
}

bool Reference::ReadFirstVariantPositions(){
	if( ReadFirstVcfRecord() ){
		variant_positions_.resize(NumberSequences());
		return ReadVariantPositions(2);
	}

	return false;
}

void Reference::ClearVariants(uint32_t end_ref_seq_id){
	if(VariantsLoaded()){
		while( cleared_variation_for_num_sequences_ < end_ref_seq_id ){
			variants_.at(cleared_variation_for_num_sequences_).clear();
			variants_.at(cleared_variation_for_num_sequences_++).shrink_to_fit();
		}
	}
}

void Reference::ClearVariantPositions(uint32_t end_ref_seq_id){
	if(VariantPositionsLoaded()){
		while( cleared_variation_for_num_sequences_ < end_ref_seq_id ){
			variant_positions_.at(cleared_variation_for_num_sequences_).clear();
			variant_positions_.at(cleared_variation_for_num_sequences_++).shrink_to_fit();
		}
	}
}

void Reference::ClearAllVariants(){
	if(VariantsLoaded()){
		CloseVcfFile();
		variants_.clear();
		variants_.shrink_to_fit();
	}
}

void Reference::ClearAllVariantPositions(){
	if(VariantPositionsLoaded()){
		CloseVcfFile();
		variant_positions_.clear();
		variant_positions_.shrink_to_fit();
	}
}
