#include "SimulatorTest.h"
using reseq::SimulatorTest;

#include <algorithm>
using std::max;
//include <array>
using std::array;
#include <bitset>
using std::bitset;
#include <string>
using std::string;
//include <vector>
using std::vector;

//include <seqan/seq_io.h>
using seqan::DnaString;
using seqan::length;

//include "utilities.hpp"
using reseq::utilities::ReverseComplementorDna;
using reseq::utilities::at;
using reseq::utilities::Percent;

void SimulatorTest::Register(){
	// Guarantees that library is included
}

void SimulatorTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new Simulator ) << "Could not allocate memory for FragmentDuplicationStats object\n";
}

void SimulatorTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void SimulatorTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void SimulatorTest::ChooseAlleles(vector<uintAlleleId> &chosen_allele_ids, vector<bool> &reverse_selection, uintAlleleId non_zero_strands, uintAlleleId possible_strands) const{
	EXPECT_FALSE( non_zero_strands <= possible_strands/2 ) << "Not all alleles were select in the test. Something went wrong";

	// More than half the possible alleles need to be drawn: Draw inverse selection
	chosen_allele_ids.clear();
	reverse_selection.clear();
	reverse_selection.resize(possible_strands, true);
	while( chosen_allele_ids.size() < possible_strands-non_zero_strands ){
		double random_value = 0.5;
		test_->SelectAllele(chosen_allele_ids, reverse_selection, possible_strands, random_value);
	}

	EXPECT_EQ(possible_strands-non_zero_strands, chosen_allele_ids.size());

	test_->ReverseSelection(chosen_allele_ids, reverse_selection, possible_strands);

	EXPECT_EQ(possible_strands, chosen_allele_ids.size());
}

void SimulatorTest::TestCoverageConversion(){
	// CoveragePropLostFromAdapters
	DataStats stats(NULL);

	stats.read_lengths_by_fragment_length_.at(1)[200][150] = 10;
	stats.read_lengths_by_fragment_length_.at(0)[150][150] = 10;
	stats.read_lengths_by_fragment_length_.at(1)[100][150] = 5; // 5*50
	stats.read_lengths_by_fragment_length_.at(1)[50][150] = 10; // (10-5)*100 [-5 due to next line]
	stats.non_mapped_read_lengths_by_fragment_length_.at(1)[50][150] = 5; // 5*150
	stats.read_lengths_by_fragment_length_.at(0)[0][150] = 10; // 10*150

	double adapter_part = test_->CoveragePropLostFromAdapters(stats);
	EXPECT_DOUBLE_EQ( 3000.0/(45*150), adapter_part );

	// CoverageToNumberPairs
	double coverage = 100;
	uintRefLenCalc total_ref_size = 50000;
	double average_read_length = 150;

	uintFragCount total_pairs = test_->CoverageToNumberPairs(coverage, total_ref_size, average_read_length, adapter_part);
	EXPECT_EQ( 30000, total_pairs );

	// NumberPairsToCoverage
	EXPECT_NEAR( coverage, test_->NumberPairsToCoverage(total_pairs, total_ref_size, average_read_length, adapter_part), 0.01 );
}

void SimulatorTest::TestSelectAllele(){
	// Test 2 alleles
	vector<uintAlleleId> chosen_allele_ids;
	vector<bool> reverse_selection;

	reverse_selection.resize(2, true);
	while( chosen_allele_ids.size() < 2 ){
		test_->SelectAllele(chosen_allele_ids, reverse_selection, 2, 0.5);
	}

	EXPECT_EQ(1, chosen_allele_ids.at(0));
	EXPECT_EQ(0, chosen_allele_ids.at(1));

	// Test 4 alleles
	chosen_allele_ids.clear();
	reverse_selection.clear();
	reverse_selection.resize(4, true);
	while( chosen_allele_ids.size() < 4 ){
		test_->SelectAllele(chosen_allele_ids, reverse_selection, 4, 0.5);
	}

	EXPECT_EQ(2, chosen_allele_ids.at(0));
	EXPECT_EQ(1, chosen_allele_ids.at(1));
	EXPECT_EQ(3, chosen_allele_ids.at(2));
	EXPECT_EQ(0, chosen_allele_ids.at(3));
}

void SimulatorTest::TestVariationInInnerLoopOfSimulateFromGivenBlock(
		Simulator::VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		uintSeqLen cur_start_position,
		uintSeqLen frag_length_from,
		uintSeqLen frag_length_to,
		uintAlleleId num_valid_alleles,
		array<vector<intVariantId>, 2> unhandled_variant_id,
		array<vector<uintSeqLen>, 2> unhandled_bases_in_variant,
		array<vector<intSeqShift>, 2> gc_mod,
		array<vector<intSeqShift>, 2> end_pos_shift,
		array<uintSeqLen, 2> modified_start_pos,
		array<Reference *, 2> comp_ref){
	// Preparation
	DataStats stats(NULL);
	stats.read_lengths_.at(0)[100];
	stats.read_lengths_.at(1)[100];
	stats.errors_.PrepareSimulation();
	Simulator::SimPair sim_reads;
	Surrounding surrounding_end, comp_surrounding;

	auto num_tests = 0;

	uintPercent gc_perc;
	uintSeqLen cur_end_position;

	vector<uintAlleleId> possible_alleles, chosen_allele_ids;
	vector<bool> reverse_selection;

	double probability_chosen = 1.0;

	// Inner loop
	auto frag_len_start = frag_length_from;
	test_->GetPossibleAlleles(possible_alleles, species_reference_, bias_mod, cur_start_position, ref_seq_id);

	for( auto fragment_length=frag_len_start; fragment_length < frag_length_to; ++fragment_length){
		if(test_->ProbabilityAboveThreshold(probability_chosen, ref_seq_id, fragment_length)){
			auto non_zero_strands = stats.FragmentDistribution().DrawNumberNonZeroStrands( possible_alleles.size(), test_->NonZeroThreshold(ref_seq_id, fragment_length), probability_chosen );
			EXPECT_EQ(2*possible_alleles.size(), non_zero_strands);

			if( non_zero_strands ){
				ChooseAlleles(chosen_allele_ids, reverse_selection, non_zero_strands, 2*possible_alleles.size());

				for(auto chosen_id : chosen_allele_ids){
					auto allele = possible_alleles.at(chosen_id/2);
					bool strand = chosen_id%2;

					test_->PrepareBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, species_reference_, cur_start_position, fragment_length, allele );

					comp_ref.at(allele)->ReverseSurrounding( comp_surrounding, ref_seq_id, modified_start_pos.at(allele)+fragment_length-1 );
					EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_end_.at(allele).sur_.at(0)) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Allele: " << allele << std::endl;
					EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_end_.at(allele).sur_.at(1)) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Allele: " << allele << std::endl;
					EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_end_.at(allele).sur_.at(2)) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Allele: " << allele << std::endl;

					auto result = fragment_length - frag_length_from;
					EXPECT_EQ( unhandled_variant_id.at(allele).at(result), bias_mod.unhandled_variant_id_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
					EXPECT_EQ( unhandled_bases_in_variant.at(allele).at(result), bias_mod.unhandled_bases_in_variant_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
					EXPECT_EQ( gc_mod.at(allele).at(result), bias_mod.gc_mod_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
					EXPECT_EQ( end_pos_shift.at(allele).at(result), bias_mod.end_pos_shift_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;

					cur_end_position = cur_start_position + fragment_length + bias_mod.end_pos_shift_.at(allele);
					if( cur_end_position <= species_reference_.SequenceLength(ref_seq_id)){
						// Determine how many read pairs are generated for this strand and allele at this position with this fragment_length
						gc_perc = test_->GetGCPercent( bias_mod, ref_seq_id, species_reference_, cur_end_position, fragment_length, allele );
						EXPECT_EQ( Percent(comp_ref.at(allele)->GCContentAbsolut( ref_seq_id, modified_start_pos.at(allele), modified_start_pos.at(allele)+fragment_length ), fragment_length), gc_perc );

						test_->GetOrgSeq(sim_reads, strand, allele, fragment_length, cur_start_position, cur_end_position, ref_seq_id, species_reference_, stats, bias_mod);

						EXPECT_TRUE( infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id), modified_start_pos.at(allele), modified_start_pos.at(allele)+fragment_length) == prefix(sim_reads.at(strand).org_seq_, fragment_length) ) << prefix(sim_reads.at(strand).org_seq_, fragment_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Allele: " << allele << std::endl;
						EXPECT_TRUE( ReverseComplementorDna(infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id), modified_start_pos.at(allele), modified_start_pos.at(allele)+fragment_length)) == prefix(sim_reads.at(!strand).org_seq_, fragment_length) ) << prefix(sim_reads.at(!allele).org_seq_, fragment_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << fragment_length << " Allele: " << allele << std::endl;

						++num_tests;
					}
				}
			}
		}
	}

	EXPECT_EQ( num_valid_alleles*2*max(unhandled_variant_id.at(0).size(),unhandled_variant_id.at(1).size()), num_tests );
}

void SimulatorTest::TestVariationInSimulateFromGivenBlock(){
	species_reference_.num_alleles_ = 2;

	// samtools faidx ../test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:1001-1020
	// GTTGCGAGATTTGGACGGAC
	species_reference_.variants_.clear();
	species_reference_.variants_.resize(1);
	species_reference_.variants_.at(0).clear();
	std::array<uintAlleleBitArray, reseq::Reference::Variant::kMaxAlleles/64> present_in_alleles = {2}; //bitwise [0:no, 1:yes]
	species_reference_.variants_.at(0).emplace_back(455, "", present_in_alleles);
	species_reference_.variants_.at(0).emplace_back(1002, "TAC", present_in_alleles);
	species_reference_.variants_.at(0).emplace_back(1004, "TGA", present_in_alleles);
	species_reference_.variants_.at(0).emplace_back(1008, "", present_in_alleles);
	species_reference_.variants_.at(0).emplace_back(1011, "C", present_in_alleles);
	species_reference_.variants_.at(0).emplace_back(1012, "", present_in_alleles);
	// GTT--GC--GAGATTTGGACGGAC
	// GTTACGCGAGAG-TTC-GACGGAC

	Reference test_ref;
	resize( test_ref.reference_sequences_, 1 );
	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 455);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 456, 1003);
	at(test_ref.reference_sequences_, 0) += "AC";
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1003, 1004);
	at(test_ref.reference_sequences_, 0) += "TGA";
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1005, 1008);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1009, 1012);
	at(at(test_ref.reference_sequences_, 0), 1013) = 'C';
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1013, 2000);
	Surrounding comp_surrounding;

	test_->coverage_groups_.resize(species_reference_.NumberSequences(), 0);
	test_->non_zero_thresholds_.resize(1);
	test_->non_zero_thresholds_.at(0).resize(100, {0.75, 0.31640625});

	Simulator::VariantBiasVarModifiers bias_mod(1, 2);
	uintRefSeqId ref_seq_id = 0;
	Surrounding surrounding_start;

	// Position 1003
	uintSeqLen cur_start_position = 1003;
	bias_mod.first_variant_id_ = 2; // Insertion at position 2 has already been processed

	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1004 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 13, 2,
			{{ {2,3,3,3,3,4,4,4,5,6,6,6}, {2,2,2,3,3,3,3,4,4,5,6,6} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0,0}, {0,2,1,0,0,0,0,0,0,0,0,0} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0,0}, {0,-1,0,0,0,0,0,0,0,1,0,0} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0,0}, {0,0,-1,-2,-2,-2,-2,-1,-1,-1,0,0} }},
			{cur_start_position, 1004}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1004 + 0
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, ++cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(0)) << " (Result)";
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(1)) << " (Result)";
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(2)) << " (Result)";

	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1005 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 12, 2,
			{{ {3,3,3,3,4,4,4,5,6,6,6}, {2,2,3,3,3,3,4,4,5,6,6} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0}, {2,1,0,0,0,0,0,0,0,0,0} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0}, {-1,0,0,0,0,0,0,0,1,0,0} }},
			{{ {0,0,0,0,0,0,0,0,0,0,0}, {0,-1,-2,-2,-2,-2,-1,-1,-1,0,0} }},
			{cur_start_position, 1005}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(1, bias_mod.start_variant_pos_);

	// Position 1004 + 1
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1006 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 11, 1,
			{{ {}, {3,3,3,3,3,4,4,5,6,6} }},
			{{ {}, {0,0,0,0,0,0,0,0,0,0} }},
			{{ {}, {0,0,0,0,0,0,0,1,0,0} }},
			{{ {}, {0,-1,-1,-1,-1,0,0,0,1,1} }},
			{0, 1006}, {NULL, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(2, bias_mod.start_variant_pos_);

	// Position 1004 + 2
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1007 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 10, 1,
			{{ {}, {3,3,3,3,4,4,5,6,6} }},
			{{ {}, {0,0,0,0,0,0,0,0,0} }},
			{{ {}, {-1,-1,-1,-1,-1,-1,0,-1,-1} }},
			{{ {}, {0,0,0,0,1,1,1,2,2} }},
			{0, 1007}, {NULL, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(3, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1008
	cur_start_position = 1008;
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 8, 1,
			{{ {4,4,4,5,6,6,6}, {} }},
			{{ {0,0,0,0,0,0,0}, {} }},
			{{ {0,0,0,0,0,0,0}, {} }},
			{{ {0,0,0,0,0,0,0}, {} }},
			{cur_start_position, 0}, {&species_reference_, NULL});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(4, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1011
	cur_start_position = 1011;
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, 1, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1013 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 5, 2,
			{{ {5,6,6,6}, {5,6,6,6} }},
			{{ {0,0,0,0}, {0,0,0,0} }},
			{{ {0,0,0,0}, {1,0,0,0} }},
			{{ {0,0,0,0}, {0,1,1,1} }},
			{cur_start_position, 1013}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(5, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);
}

namespace reseq{
	TEST_F(SimulatorTest, BasicFunctonality){
		CreateTestObject();

		TestCoverageConversion();
		TestSelectAllele();
	}

	TEST_F(SimulatorTest, Variants){
		CreateTestObject();
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		LoadReference(test_dir+"ecoli-GCF_000005845.2_ASM584v2_genomic.fa");

		TestVariationInSimulateFromGivenBlock();
	}
}
