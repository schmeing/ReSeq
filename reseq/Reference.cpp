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
#include <iterator>
using std::distance;
#include <limits>
using std::numeric_limits;
#include <random>
using std::mt19937_64;
using std::uniform_int_distribution;
//include <set>
using std::set;
#include <string>
using std::string;
#include <utility>
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
using seqan::IupacString;
using seqan::length;
using seqan::Prefix;
using seqan::prefix;
using seqan::readRecord;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
//include <seqan/vcf_io.h>
using seqan::VcfFileIn;
using seqan::VcfHeader;
using seqan::VcfRecord;

#include "CMakeConfig.h"
//include "Surrounding.h"
using reseq::Surrounding;
using reseq::SurroundingBias;
//include "utilities.hpp"
using reseq::utilities::ComplementedConstDna5String;
using reseq::utilities::at;
using reseq::utilities::DivideAndCeil;
using reseq::utilities::IsN;
using reseq::utilities::Percent;
using reseq::utilities::ReverseComplementorDna;
using reseq::utilities::SetToMin;


inline void Reference::PushBackExclusionRegion( pair<uintSeqLen, uintSeqLen> region, uintRefSeqId ref_seq, uintSeqLen maximum_fragment_length){
	// If not all fragment length fit between the two regions, merge them
	if(excluded_regions_.at(ref_seq).back().second + maximum_fragment_length <= region.first){
		// Distance is big enough, insert new region
		excluded_regions_.at(ref_seq).push_back(region);
	}
	else{
		// Merge regions
		auto it_first_merged_region = excluded_regions_.at(ref_seq).rbegin() + 1;
		while( it_first_merged_region != excluded_regions_.at(ref_seq).rend() && (it_first_merged_region)->second + maximum_fragment_length > region.first ){
			++it_first_merged_region;
		}
		--it_first_merged_region; // Go back to the element that needed merging

		// Set new values to the first region that still needs to be merged
		SetToMin(it_first_merged_region->first, region.first);
		it_first_merged_region->second = max(excluded_regions_.at(ref_seq).back().second, region.second);

		// Remove all the regions after it that are merged
		excluded_regions_.at(ref_seq).resize( excluded_regions_.at(ref_seq).size() - distance(excluded_regions_.at(ref_seq).rbegin(), it_first_merged_region) );
	}
}

bool Reference::CheckVcf() const{
	auto contigs = contigNames(context(vcf_file_));
	uintErrorCount errors = 0;

	if(length(contigs) != NumberSequences()){
		printErr << "Number of contigs does not match between reference(" << NumberSequences() << ") and variant(" << length(contigs) << ") file." << std::endl;
		++errors;
	}
	for(uintRefSeqId con=0; con < min(static_cast<uintRefSeqId>(length(contigs)), NumberSequences()); ++con){
		if(errors < 20 && at(contigs, con) != ReferenceIdFirstPart(con)){
			printErr << "Contigs at position " << con << " do not match between reference(" << ReferenceId(con) << ") and variant(" << at(contigs, con) << ") file." << std::endl;
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

bool Reference::ReadVariants(uintRefSeqId end_ref_seq_id, bool positions_only){
	// Read the file record by record until variant is for sequence end_ref_seq_id or later
	uintErrorCount errors=0;
	bool tmp_success = true;
	uintSeqLen start_pos, end_pos(0);

	uintSeqLen pos;
	uintAlleleId cur_allele, chosen_var;
	vector<uintAlleleId> allele;
	vector<uintAlleleBitArray> gt_has_var;
	vector<uintSeqLen> alt_start_pos;

	if(!positions_only){
		allele.resize(num_alleles_);
		gt_has_var.resize(num_alleles_);
	}

	uintRefSeqId old_ref_id = numeric_limits<uintRefSeqId>::max();

	try{
		while(VcfRecord::INVALID_REFID != cur_vcf_record_.rID && cur_vcf_record_.rID < end_ref_seq_id){
			// Parse record
			if(cur_vcf_record_.rID >= NumberSequences()){
				printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " does not belong to an existing reference sequence." << std::endl;
				if(++errors >= kMaxErrorsShownPerFile){
					printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
					break;
				}
			}
			else{
				if(cur_vcf_record_.beginPos >= SequenceLength(cur_vcf_record_.rID)){
					printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " starts after the end of the reference sequence." << std::endl;
					if(++errors >= kMaxErrorsShownPerFile){
						printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
						break;
					}
				}
				else{
					start_pos = cur_vcf_record_.beginPos;

					if(old_ref_id == cur_vcf_record_.rID){
						if(start_pos < end_pos){
							printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " overlaps with a previous variant." << std::endl;
							if(++errors >= kMaxErrorsShownPerFile){
								printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
								break;
							}
						}
					}
					else{
						old_ref_id = cur_vcf_record_.rID;
					}

					end_pos = start_pos + length(cur_vcf_record_.ref);
					Dna5String vcf_ref_var_ = cur_vcf_record_.ref;

					if( utilities::HasN(vcf_ref_var_) ){
						printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " has an reference column containing ambiguous bases (e.g. N). Please change or remove them, but make sure the reference file stays consistent with this column." << std::endl;
						if(++errors >= kMaxErrorsShownPerFile){
							printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
							break;
						}
					}
					else{
						if( vcf_ref_var_ != infix(ReferenceSequence(cur_vcf_record_.rID), start_pos, end_pos) ){
							printErr << "The specified reference in vcf file '" << vcf_ref_var_ << "' is not identical with the specified reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << ": '" << infix(ReferenceSequence(cur_vcf_record_.rID), start_pos, end_pos ) << "'." << std::endl;
							if(++errors >= kMaxErrorsShownPerFile){
								printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
								break;
							}
						}
					}

					if(positions_only){
						for(auto pos=start_pos; pos<end_pos; ++pos){
							variant_positions_.at(cur_vcf_record_.rID).push_back(pos);
						}
					}
					else{
						// Get genotypes
						tmp_success = true;
						cur_allele = 0;
						for( auto &genotype : cur_vcf_record_.genotypeInfos ){
							if(cur_allele >= num_alleles_){
								tmp_success = false;
								printErr << "Found to many alleles in genotype definition '";
								for( auto &genotype : cur_vcf_record_.genotypeInfos ){
									std::cout << ' ' << genotype;
								}
								std::cout << "'" << std::endl;
								if(++errors >= kMaxErrorsShownPerFile){
									printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
								}
								break;
							}

							pos=0;
							chosen_var=0;
							while(pos < length(genotype) && ':' != at(genotype, pos)){
								if('|' == at(genotype, pos) || '/' == at(genotype, pos)){
									allele.at(cur_allele++) = chosen_var;
									chosen_var = 0;
								}
								else if('0' <= at(genotype, pos) && '9' >= at(genotype, pos)){
									chosen_var *= 10;
									chosen_var += at(genotype, pos) - 48;
								}
								else{
									tmp_success = false;
									printErr << "Unallowed character '" << at(genotype, pos) << "' in genotype definition '" << genotype << "'" << std::endl;
									if(++errors >= kMaxErrorsShownPerFile){
										printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
										break;
									}
								}
								++pos;
							}

							if(errors >= kMaxErrorsShownPerFile){
								break;
							}

							if(tmp_success){
								allele.at(cur_allele++) = chosen_var;
							}
							else{
								allele.at(cur_allele++) = 0;
							}
						}

						if(errors >= kMaxErrorsShownPerFile){
							break;
						}

						if(cur_allele < num_alleles_){
							tmp_success = false;
							printErr << "Could not find enough alleles in genotype definition '";
							for( auto &genotype : cur_vcf_record_.genotypeInfos ){
								std::cout << ' ' << genotype;
							}
							std::cout << "'" << std::endl;
							if(++errors >= kMaxErrorsShownPerFile){
								printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
								break;
							}
						}

						if(tmp_success){
							// Find variations for the genotypes
							gt_has_var.clear();
							alt_start_pos.clear();
							alt_start_pos.push_back(0);

							if(length(cur_vcf_record_.alt) > numeric_limits<uintSeqLen>::max()){
								printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " has an alternative column of length " << length(cur_vcf_record_.alt) << ", but currently only a maximum of " << numeric_limits<uintSeqLen>::max() << " characters are supported." << std::endl;
								if(++errors >= kMaxErrorsShownPerFile){
									printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
									break;
								}
							}
							else{
								chosen_var = 1; // 0 is the reference sequence
								// Split alternative sequences
								for(pos = 0; pos < length(cur_vcf_record_.alt); ++pos){
									if(',' == at(cur_vcf_record_.alt, pos)){
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
										if(++errors <= kMaxErrorsShownPerFile){
											printErr << "Variant number " << allele.at(cur_allele) << " does not exist for sequence id " << cur_vcf_record_.rID << " and position " << cur_vcf_record_.beginPos << std::endl;
										}
										if(errors >= kMaxErrorsShownPerFile){
											printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
										}
									}
								}
								++chosen_var;

								// Enter variants into vector (split into single reference positions)
								Dna5String inserted_variant;
								for(pos = 0; pos < length(vcf_ref_var_); ++pos){
									for(uintAlleleId n_alt = 0; n_alt < gt_has_var.size(); ++n_alt){
										// Check if any genotype has this variant
										if(gt_has_var.at(n_alt)){
											// Enter record into variations vector
											if(pos+1 == length(vcf_ref_var_) && pos+1 < alt_start_pos.at(n_alt+1)-1 - alt_start_pos.at(n_alt)){
												// Insertion
												inserted_variant = infix(cur_vcf_record_.alt, alt_start_pos.at(n_alt)+pos, alt_start_pos.at(n_alt+1)-1);

											}
											else if(pos < alt_start_pos.at(n_alt+1)-1 - alt_start_pos.at(n_alt)){
												// Base Mutation
												if(at(vcf_ref_var_, pos) != at(cur_vcf_record_.alt, alt_start_pos.at(n_alt)+pos)){
													inserted_variant = at(cur_vcf_record_.alt, alt_start_pos.at(n_alt)+pos);
												}
												else{
													continue; // If variant is identical to reference, we don't need to add it
												}
											}
											else{
												// Deletion
												inserted_variant = Dna5String("");
											}

											if( utilities::HasN(inserted_variant) ){
												if(++errors <= kMaxErrorsShownPerFile){
													printErr << "Variant starting in reference sequence " << cur_vcf_record_.rID << " at position " << start_pos << " has an alternative column containing ambiguous bases (e.g. N). Please change or remove them." << std::endl;
												}
												if(errors >= kMaxErrorsShownPerFile){
													printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
												}
											}
											else{
												InsertVariant(cur_vcf_record_.rID, start_pos+pos, inserted_variant, gt_has_var.at(n_alt));
											}
										}
									}
								}
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
					if(++errors >= kMaxErrorsShownPerFile){
						printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
						break;
					}
					else{
						continue;
					}
				}
				else if(cur_vcf_record_.rID == read_variation_for_num_sequences_){
					if(cur_vcf_record_.beginPos < start_pos){
						printErr << "Variant file is not properly position sorted. Found in sequence id " << cur_vcf_record_.rID << " position " << cur_vcf_record_.beginPos << " after position " << start_pos << std::endl;
						if(++errors >= kMaxErrorsShownPerFile){
							printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
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
		printErr << "Variant file is not sorted by position. The following variation record and all after it have been ignored during simulation:" << std::endl;
		std::cout << cur_vcf_record_.rID << '\t' << cur_vcf_record_.beginPos+1 << '\t' << cur_vcf_record_.id << '\t' << cur_vcf_record_.ref << '\t' << cur_vcf_record_.alt << '\t' << cur_vcf_record_.qual << '\t' << cur_vcf_record_.filter << '\t' << cur_vcf_record_.info << '\t' << cur_vcf_record_.format;
		for(auto &genotype : cur_vcf_record_.genotypeInfos){
			std::cout << '\t' << genotype;
		}
		std::cout << std::endl;
	}
	clear(contigNames(context(vcf_file_)));
	close(vcf_file_);
}

inline void Reference::AddFragmentSite(
		std::vector<FragmentSite> &sites,
		uintSeqLen fragment_length,
		uintSeqLen gc,
		uintSeqLen n_count,
		const Surrounding &start_sur,
		const Surrounding &end_sur
		) const{
	if(n_count <= kMaxNInFragmentSite && start_sur.Valid() && end_sur.Valid()){
		uintPercent gc_perc = Percent(gc, fragment_length-n_count);
		sites.emplace_back(gc_perc, start_sur, end_sur);
	}
	else{
		sites.emplace_back(); // Just a placeholder, so that the indexes in sites match starting positions(with the shift due to kMinDistToContigEnds), will be removed after counts have been entered
	}
}

Reference::Reference():
	num_alleles_(1) // In case we don't load a variant file, we have a single allele
	{
}

const Prefix<const CharString>::Type Reference::ReferenceIdFirstPart( uintRefSeqId n ) const{
	uintSeqLen pos=0; // declare before the for loop to be able to access it afterwards
	for( ; pos < length(ReferenceId(n)) && ' ' != at(ReferenceId(n), pos); ++pos ); // loop until pos is at the first ' ' or at the end of the string
	return prefix(ReferenceId(n), pos);
}

void Reference::ReferenceSequence(
		DnaString &insert_string,
		uintRefSeqId seq_id,
		uintSeqLen start_pos,
		uintSeqLen frag_length,
		bool reversed) const{
	if( reversed ){
		insert_string = infix(ReferenceSequence(seq_id), start_pos-frag_length, start_pos);
		reverseComplement(insert_string);
	}
	else{
		insert_string = infix(ReferenceSequence(seq_id), start_pos, start_pos+frag_length);
	}
}

void Reference::ReferenceSequence(
		DnaString &insert_string,
		uintRefSeqId seq_id,
		uintSeqLen start_pos,
		uintSeqLen frag_length,
		bool reversed,
		const vector<Variant> &variants,
		pair<intVariantId, uintSeqLen> first_variant, // {id, posCurrentlyAt}
		uintAlleleId allele) const{
	if( reversed ){
		insert_string = "";
		auto cur_start = start_pos;
		auto cur_var = first_variant.first;
		if(first_variant.second){
			insert_string += ReverseComplementorDna(prefix(variants.at(cur_var).var_seq_, first_variant.second));
			--cur_var;
			--cur_start;
		}
		for( ; cur_var >= 0 && length(insert_string) < frag_length; --cur_var ){
			if(variants.at(cur_var).InAllele(allele)){
				if(cur_start-variants.at(cur_var).position_ > frag_length-length(insert_string)){
					// Variant after return sequence
					insert_string += ReverseComplementorDna(infix(ReferenceSequence(seq_id), cur_start+length(insert_string)-frag_length, cur_start));
				}
				else{
					// Variant inside return sequence
					insert_string += ReverseComplementorDna(infix(ReferenceSequence(seq_id), variants.at(cur_var).position_+1, cur_start));
					insert_string += ReverseComplementorDna(variants.at(cur_var).var_seq_);
					cur_start = variants.at(cur_var).position_;
				}
			}
		}
		if(cur_var == -1 && length(insert_string) < frag_length){
			// Fill rest after the last variant
			insert_string += ReverseComplementorDna(infix(ReferenceSequence(seq_id), cur_start+length(insert_string)-frag_length, cur_start));
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
		for( ; cur_var < variants.size() && length(insert_string) < frag_length; ++cur_var ){
			if(variants.at(cur_var).InAllele(allele)){
				if(variants.at(cur_var).position_-cur_start >= frag_length-length(insert_string)){
					// Variant after return sequence
					insert_string += infix(ReferenceSequence(seq_id), cur_start, cur_start+frag_length-length(insert_string));
				}
				else{
					// Variant inside return sequence
					insert_string += infix(ReferenceSequence(seq_id), cur_start, variants.at(cur_var).position_);
					insert_string += variants.at(cur_var).var_seq_;
					cur_start = variants.at(cur_var).position_+1;
				}
			}
		}
		if(cur_var == variants.size() && length(insert_string) < frag_length){
			// Fill rest after the last variant
			insert_string += infix(ReferenceSequence(seq_id), cur_start, cur_start+frag_length-length(insert_string));
		}
	}

	if(length(insert_string) > frag_length){ // This can happen if we insert a variant that is longer than the end of the sequence as we do not check for the length of variants
		resize(insert_string, frag_length);
	}
}

reseq::uintRefLenCalc Reference::TotalSize() const{
	uintRefLenCalc sum = 0;
	for(const auto &seq : reference_sequences_){
		sum += length(seq);
	}
	return sum;
}

void Reference::RefSeqsSortedByNxx(vector<pair<double, uintRefSeqId>> &ref_seqs) const{
	double tot_size = TotalSize();
	tot_size /= 100; // Nxx in percent

	ref_seqs.reserve( NumberSequences() );
	for(uintRefSeqId id=0; id < NumberSequences(); ++id){
		ref_seqs.emplace_back( SequenceLength(id)/tot_size , id );
	}

	sort(ref_seqs.begin(), ref_seqs.end(), std::greater<pair<double, uintRefSeqId>>());

	double sum = 0.0;
	for(uintRefSeqId id=0; id < NumberSequences(); ++id){
		sum += ref_seqs.at(id).first;
		ref_seqs.at(id).first = sum;
	}
}

reseq::uintRefSeqId Reference::RefSeqsInNxx(vector<bool> &ref_seqs, double nxx_ref_seqs) const{
	// Get nxx for all reference sequences
	vector<pair<double, uintRefSeqId>> ref_seq_nxx;
	RefSeqsSortedByNxx(ref_seq_nxx);

	// Count how many sequences are needed to get nxx_ref_seqs taking only the largest ones
	uintRefSeqId needed_ref_seqs = 0;
	ref_seqs.resize(NumberSequences(), false);
	while(needed_ref_seqs+1 < ref_seq_nxx.size() && nxx_ref_seqs > ref_seq_nxx.at(needed_ref_seqs).first){
		ref_seqs.at(needed_ref_seqs++) = true;
	}

	ref_seqs.at(needed_ref_seqs++) = true; // Get last one and convert index to size

	return needed_ref_seqs;
}

reseq::uintRefSeqBin Reference::NumRefSeqBinsInNxx(const std::vector<bool> &ref_seqs, uintSeqLen max_ref_seq_bin_length) const{
	uintRefSeqBin num_bins = 0;
	for(uintRefSeqId ref_id=0; ref_id<ref_seqs.size(); ++ref_id){
		if(ref_seqs.at(ref_id)){
			num_bins += DivideAndCeil(SequenceLength(ref_id),  max_ref_seq_bin_length);
		}
	}
	return num_bins;
}

double Reference::SumBias(
		double &max_bias,
		uintRefSeqId ref_seq_id,
		uintSeqLen fragment_length,
		double general_bias,
		const Vect<double> &gc_bias,
		const SurroundingBias &sur_bias) const{
	// Ignores exclusion regions and does not account for N's as it is used for simulation only
	double tot(0.0);

	uintSeqLen n_count(0);
	uintSeqLen gc = GCContentAbsolut(n_count, ref_seq_id, 0, fragment_length);

	Surrounding start_sur, end_sur;
	ForwardSurrounding(start_sur, ref_seq_id, 0);
	ReverseSurrounding(end_sur, ref_seq_id, fragment_length-1);

	double bias = Bias(general_bias, gc_bias[Percent(gc, fragment_length)], sur_bias.Bias(start_sur), sur_bias.Bias(end_sur));
	if(bias > max_bias){
		max_bias = bias;
	}
	tot += bias;

	const Dna5String &ref_seq(ReferenceSequence(ref_seq_id));
	for( uintSeqLen start_pos=0; start_pos < length(ref_seq)-fragment_length; ){
		UpdateGC( gc, n_count, ref_seq, start_pos, start_pos+fragment_length );
		end_sur.UpdateReverse( ref_seq, start_pos+fragment_length );
		start_sur.UpdateForward( ref_seq, ++start_pos );

		double bias = Bias(general_bias, gc_bias[Percent(gc, fragment_length)], sur_bias.Bias(start_sur), sur_bias.Bias(end_sur));
		if(bias > max_bias){
			max_bias = bias;
		}
		tot += bias;
	}

	return tot;
}

double Reference::SumBias(
		std::vector<utilities::VectorAtomic<uintFragCount>> &gc_sites,
		uintRefSeqId ref_seq_id,
		uintSeqLen fragment_length,
		double general_bias,
		const Vect<double> &gc_bias,
		const SurroundingBias &sur_bias) const{
	// Uses exclusion regions as it is used for stats creation
	double tot(0.0);

	uintSeqLen n_count(0);
	uintSeqLen start_pos = excluded_regions_.at(ref_seq_id).front().second;
	uintSeqLen gc = GCContentAbsolut(n_count, ref_seq_id, start_pos, start_pos+fragment_length);

	Surrounding start_sur, end_sur;
	ForwardSurroundingWithN(start_sur, ref_seq_id, start_pos);
	ReverseSurroundingWithN(end_sur, ref_seq_id, start_pos+fragment_length-1);

	if(start_sur.Valid() && end_sur.Valid()){
		tot += Bias(general_bias, gc_bias[Percent(gc, fragment_length)], sur_bias.Bias(start_sur), sur_bias.Bias(end_sur));
		++gc_sites.at(Percent(gc, fragment_length));
	}

	const Dna5String &ref_seq(ReferenceSequence(ref_seq_id));
	uintSeqLen next_exclusion_region(1);
	for( ; ++start_pos < length(ref_seq); ){
		if(start_pos == excluded_regions_.at(ref_seq_id).at(next_exclusion_region).first-fragment_length+1){
			// We hit an exclusion region, go to end of it
			start_pos = excluded_regions_.at(ref_seq_id).at(next_exclusion_region++).second;

			if(start_pos >= length(ref_seq)){
				break;
			}

			gc = GCContentAbsolut(n_count, ref_seq_id, start_pos, start_pos+fragment_length);
			ForwardSurroundingWithN(start_sur, ref_seq_id, start_pos);
			ReverseSurroundingWithN(end_sur, ref_seq_id, start_pos+fragment_length-1);
		}
		else{
			UpdateGC( gc, n_count, ref_seq, start_pos-1, start_pos+fragment_length-1 );
			end_sur.UpdateReverseWithN( ref_seq, start_pos+fragment_length-1);
			start_sur.UpdateForwardWithN( ref_seq, start_pos );
		}

		if(start_sur.Valid() && end_sur.Valid()){
			// If 0==fragment_length-n_count surrounding would be invalid, so no check needed
			tot += Bias(general_bias, gc_bias[Percent(gc, fragment_length)], sur_bias.Bias(start_sur), sur_bias.Bias(end_sur));
			++gc_sites.at(Percent(gc, fragment_length));
		}
	}

	return tot;
}

void Reference::GetFragmentSites( vector<FragmentSite> &sites, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintSeqLen start, uintSeqLen end ) const{
	// Uses exclusion regions as it is used for stats creation
	sites.clear();

	const Dna5String &ref_seq(ReferenceSequence(ref_seq_id));
	auto start_pos = max(excluded_regions_.at(ref_seq_id).front().second, start);
	auto end_pos = min(excluded_regions_.at(ref_seq_id).back().first-fragment_length+1, end);

	uintSeqLen n_count(0);
	uintSeqLen gc = GCContentAbsolut(n_count, ref_seq_id, start_pos, start_pos+fragment_length);

	Surrounding start_sur, end_sur;
	ForwardSurroundingWithN(start_sur, ref_seq_id, start_pos);
	ReverseSurroundingWithN(end_sur, ref_seq_id, start_pos+fragment_length-1);

	AddFragmentSite( sites, fragment_length, gc, n_count, start_sur, end_sur );

	auto next_exclusion_region = NextExclusionRegion(ref_seq_id, start_pos);

	for( uintSeqLen frag_start=start_pos; ++frag_start < end_pos; ){
		if(frag_start == excluded_regions_.at(ref_seq_id).at(next_exclusion_region).first-fragment_length+1){
			// We hit an exclusion region, go to end of it
			frag_start = excluded_regions_.at(ref_seq_id).at(next_exclusion_region++).second;

			if(frag_start >= end_pos){
				break;
			}

			// Add placeholder for the fragment length as corrected positions only handle the excluded regions themselves
			sites.resize( sites.size()+fragment_length-1 );

			// Set values to end of excluded region
			gc = GCContentAbsolut(n_count, ref_seq_id, frag_start, frag_start+fragment_length);
			ForwardSurroundingWithN(start_sur, ref_seq_id, frag_start);
			ReverseSurroundingWithN(end_sur, ref_seq_id, frag_start+fragment_length-1);
		}
		else{
			UpdateGC( gc, n_count, ref_seq, frag_start-1, frag_start+fragment_length-1 );
			end_sur.UpdateReverseWithN( ref_seq, frag_start+fragment_length-1 );
			start_sur.UpdateForwardWithN( ref_seq, frag_start );
		}

		AddFragmentSite( sites, fragment_length, gc, n_count, start_sur, end_sur );
	}
}

bool Reference::ReadFasta(const char *fasta_file){
	bool success = true;

	SeqFileIn ref;
	seqan::StringSet<seqan::IupacString> tmp_ref_seqs;
	clear(reference_ids_);
	if( !open(ref, fasta_file) ){
		printErr << "Could not open " << fasta_file << " for reading." << std::endl;
		success = false;
	}
	else if(atEnd(ref)){
		printErr << fasta_file << " does not contain any reference sequences." << std::endl;
		success = false;
	}
	else{
		try{
			readRecords(reference_ids_, tmp_ref_seqs, ref);
			printInfo << "Read in " << length(reference_ids_) << " reference sequences." << std::endl;
		}
		catch(const Exception &e){
			printErr << "Could not read record in " << fasta_file << ": " << e.what() << std::endl;
			success = false;
		}
	}

	if(!success){
		return false;
	}

	if( length(tmp_ref_seqs) > numeric_limits<uintRefSeqId>::max() ){
		printErr << "Reference has  " << length(tmp_ref_seqs) << " sequence entries. Currently a maximum of " << numeric_limits<uintRefSeqId>::max() << " is supported." << std::endl;
		return false;
	}

	uintErrorCount errors=0;
	for( auto ref_seq_id = length(tmp_ref_seqs); ref_seq_id--;  ){
		if( length(at(tmp_ref_seqs, ref_seq_id)) > numeric_limits<uintSeqLen>::max() ){
			printErr << "Reference sequence " << at(reference_ids_, ref_seq_id) << " is " << length(at(tmp_ref_seqs, ref_seq_id)) << " bases long. Currently a maximum of " << numeric_limits<uintSeqLen>::max() << " is supported." << std::endl;

			if(++errors >= kMaxErrorsShownPerFile){
				printErr << "Maximum number of errors reached. Additional errors are not shown for this reference." << std::endl;
				return false;
			}
		}
	}

	if(errors){
		return false;
	}

	reference_sequences_ = tmp_ref_seqs;

	return true;
}

void Reference::ReplaceN( uintSeed seed ){
	mt19937_64 rgen;
	rgen.seed(seed);
	uniform_int_distribution<> rdis(0, 3);

	for( auto &seq : reference_sequences_){
		for( auto &base : seq ){
			if( IsN(base) ){
				base = rdis(rgen);
			}
		}
	}
}

bool Reference::HasN() const{
	for( auto &seq : reference_sequences_){
		utilities::HasN(seq);
	}

	return false;
}

bool Reference::WriteFasta(const char *fasta_file) const{
	bool success = true;

	SeqFileOut ref_out;
	if( !open(ref_out, fasta_file) ){
		printErr << "Could not open " << fasta_file << " for writing." << std::endl;
		success = false;
	}
	else{
		try{
			writeRecords(ref_out, reference_ids_, reference_sequences_);
			printInfo << "Wrote " << length(reference_ids_) << " reference sequences." << std::endl;
		}
		catch(const Exception &e){
			printErr << "Could not write generated sequences to " << fasta_file << ": " << e.what() << std::endl;
			success = false;
		}
	}

	return success;
}

void Reference::PrepareExclusionRegions(){
	excluded_regions_.resize( NumberSequences() );
	sum_of_excluded_bases_.resize( NumberSequences() );
	obtained_exclusion_regions_for_num_sequences_ = 0;
	cleared_exclusion_regions_for_num_sequences_ = 0;
}

void Reference::ObtainExclusionRegions( uintRefSeqId end_ref_seq_id, uintSeqLen maximum_fragment_length ){
	for(uintRefSeqId ref_seq = obtained_exclusion_regions_for_num_sequences_; ref_seq < end_ref_seq_id && ref_seq < NumberSequences(); ++ref_seq){
		excluded_regions_.at(ref_seq).clear();

		if(SequenceLength(ref_seq) < maximum_fragment_length + 2*kMinDistToContigEnds){
			excluded_regions_.at(ref_seq).emplace_back(0, SequenceLength(ref_seq)); // Sequence is too short and completely excluded
		}
		else{
			// Add first region of ref seq
			uintSeqLen pos=0;
			for(; pos < SequenceLength(ref_seq) && IsN(at(ReferenceSequence(ref_seq), pos)); ++pos); // Count starting N's
			excluded_regions_.at(ref_seq).emplace_back(0,pos+kMinDistToContigEnds);

			// Add regions in the middle of ref seq
			uintSeqLen counted_ns(0);
			for(; pos < SequenceLength(ref_seq); ++pos){
				if( IsN(at(ReferenceSequence(ref_seq), pos)) ){
					++counted_ns;
				}
				else{
					if( kMinNToSplitContigs <= counted_ns ){
						if( kMinDistToContigEnds < pos-counted_ns){
							PushBackExclusionRegion({pos-counted_ns - kMinDistToContigEnds, pos + kMinDistToContigEnds}, ref_seq, maximum_fragment_length);
						}
						else{
							PushBackExclusionRegion({0, pos + kMinDistToContigEnds}, ref_seq, maximum_fragment_length);
						}
					}

					counted_ns = 0;
				}

			}

			// Add last region of ref seq
			PushBackExclusionRegion({SequenceLength(ref_seq) - counted_ns - kMinDistToContigEnds, SequenceLength(ref_seq)}, ref_seq, maximum_fragment_length);
		}

		// Sum excluded regions
		sum_of_excluded_bases_.at(ref_seq) = 0;
		for(auto &region : excluded_regions_.at(ref_seq) ){
			sum_of_excluded_bases_.at(ref_seq) += region.second - region.first;
		}

		++obtained_exclusion_regions_for_num_sequences_;
	}
}

bool Reference::FragmentExcluded( uintSeqLen &last_region_id, uintRefSeqId &last_ref_seq, uintRefSeqId ref_seq_id, uintSeqLen fragment_start, uintSeqLen fragment_end ) const{
	// Check if fragment overlaps with an exclusion region
	if(last_ref_seq != ref_seq_id){
		if( last_ref_seq < ref_seq_id ){
			last_region_id = 0;
		}
		else{
			last_region_id = excluded_regions_.at(ref_seq_id).size()-1;
		}

		last_ref_seq = ref_seq_id;
	}

	// Find a region that starts before the fragment end
	while( excluded_regions_.at(ref_seq_id).at(last_region_id).first >= fragment_end ){ // Cannot reach negative last_region_id, because first exclusion region starts at ref seq start
		--last_region_id;
	}

	// Find a region that ends after the fragment start
	while( excluded_regions_.at(ref_seq_id).at(last_region_id).second <= fragment_start ){ // Cannot reach over size of excluded_regions_, because last exclusion region end at ref seq end
		++last_region_id;
	}

	// If an overlapping region exists last_region_id is now pointing at it
	return excluded_regions_.at(ref_seq_id).at(last_region_id).first < fragment_end; // excluded_regions_.at(last_region_id).second > fragment_start is already valid
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
			while(pos < length(genotype) && ':' != at(genotype, pos)){
				if('|' == at(genotype, pos) || '/' == at(genotype, pos)){
					++num_alleles_; // Additional alleles in a population
				}
				++pos;
			}
			++num_alleles_; // Every population has at least one allele defined
		}

		if(num_alleles_ > Variant::kMaxAlleles){
			printErr << "Currently only " << Variant::kMaxAlleles << " alleles are supported, but file has " << num_alleles_ << '.' << std::endl;
		}
		else{
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

void Reference::ClearVariants(uintRefSeqId end_ref_seq_id){
	if(VariantsLoaded()){
		while( cleared_variation_for_num_sequences_ < end_ref_seq_id ){
			variants_.at(cleared_variation_for_num_sequences_).clear();
			variants_.at(cleared_variation_for_num_sequences_++).shrink_to_fit();
		}
	}
}

void Reference::ClearVariantPositions(uintRefSeqId end_ref_seq_id){
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
