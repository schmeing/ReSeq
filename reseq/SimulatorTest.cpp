#include "SimulatorTest.h"
using reseq::SimulatorTest;

#include <array>
using std::array;
#include <bitset>
using std::bitset;
#include <string>
using std::string;
//include <vector>
using std::vector;

//include <seqan/seq_io.h>
using seqan::DnaString;

//include "utilities.hpp"
using reseq::utilities::ReverseComplementorDna;
using reseq::utilities::at;

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

void SimulatorTest::TestVariationInInnerLoopOfSimulateFromGivenBlock(
		Simulator::VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		uintSeqLen cur_start_position,
		uintSeqLen frag_length_from,
		uintSeqLen frag_length_to,
		const vector<array<intVariantId, 2>> unhandled_variant_id,
		const vector<array<uintSeqLen, 2>> unhandled_bases_in_variant,
		const vector<array<intSeqShift, 2>> gc_mod,
		const vector<array<intSeqShift, 2>> end_pos_shift,
		array<uintSeqLen, 2> modified_start_pos,
		array<const Reference *, 2> comp_ref){
	DataStats stats(NULL);
	stats.read_lengths_.at(0)[100];
	stats.read_lengths_.at(1)[100];
	stats.fragment_distribution_.insert_lengths_[frag_length_to-1];
	stats.errors_.PrepareSimulation();
	Simulator::SimPair sim_reads;
	Surrounding surrounding_end, comp_surrounding;

	uintSeqLen cur_end_position = cur_start_position+frag_length_from-1;
	species_reference_.ReverseSurrounding( surrounding_end, ref_seq_id, cur_end_position-1 ); // -1 because it is shifted first thing in the loop

	for(auto frag_length = frag_length_from; frag_length < frag_length_to; ++frag_length){
		surrounding_end.UpdateReverse( species_reference_.ReferenceSequence(ref_seq_id), cur_end_position++ );

		do{ //while( bias_mod.fragment_length_extension_ )
			test_->UpdateBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_end_position, surrounding_end );

			auto result = frag_length - frag_length_from;
			EXPECT_EQ( unhandled_variant_id.at(result).at(0), bias_mod.unhandled_variant_id_.at(0) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( unhandled_variant_id.at(result).at(1), bias_mod.unhandled_variant_id_.at(1) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( unhandled_bases_in_variant.at(result).at(0), bias_mod.unhandled_bases_in_variant_.at(0) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( unhandled_bases_in_variant.at(result).at(1), bias_mod.unhandled_bases_in_variant_.at(1) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( gc_mod.at(result).at(0), bias_mod.gc_mod_.at(0) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( gc_mod.at(result).at(1), bias_mod.gc_mod_.at(1) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( end_pos_shift.at(result).at(0), bias_mod.end_pos_shift_.at(0) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ( end_pos_shift.at(result).at(1), bias_mod.end_pos_shift_.at(1) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;

			comp_ref.at(0)->ReverseSurrounding( comp_surrounding, ref_seq_id, modified_start_pos.at(0)+frag_length-1 );
			EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_end_.at(0).sur_.at(0)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_end_.at(0).sur_.at(1)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_end_.at(0).sur_.at(2)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			comp_ref.at(1)->ReverseSurrounding( comp_surrounding, ref_seq_id, modified_start_pos.at(1)+frag_length-1 );
			EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_end_.at(1).sur_.at(0)) << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_end_.at(1).sur_.at(0)) << " (Result)" << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_end_.at(1).sur_.at(1)) << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_end_.at(1).sur_.at(1)) << " (Result)" << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_end_.at(1).sur_.at(2)) << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_end_.at(1).sur_.at(2)) << " (Result)" << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;

			for(uintAlleleId allele=0; allele < 2; ++allele){
				test_->GetOrgSeq(sim_reads, allele, allele, frag_length, cur_start_position, cur_end_position, ref_seq_id, species_reference_, stats, bias_mod);

				EXPECT_TRUE( infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id), modified_start_pos.at(allele), modified_start_pos.at(allele)+frag_length) == prefix(sim_reads.org_seq_.at(allele), frag_length) ) << prefix(sim_reads.org_seq_.at(allele), frag_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
				EXPECT_TRUE( ReverseComplementorDna(infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id), modified_start_pos.at(allele), modified_start_pos.at(allele)+frag_length)) == prefix(sim_reads.org_seq_.at(!allele), frag_length) ) << prefix(sim_reads.org_seq_.at(!allele), frag_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
			}

			test_->CheckForFragmentLengthExtension( bias_mod, frag_length, frag_length_to, species_reference_, stats );
			EXPECT_EQ(0, bias_mod.fragment_length_extension_);
		} while(bias_mod.fragment_length_extension_);
	}
}

void SimulatorTest::TestVariationInInnerLoopOfSimulateFromGivenBlock(
		Simulator::VariantBiasVarModifiers &bias_mod,
		uintRefSeqId ref_seq_id,
		uintSeqLen cur_start_position,
		uintSeqLen frag_length_from,
		uintSeqLen frag_length_to,
		uintAlleleId allele,
		const vector<intVariantId> unhandled_variant_id,
		const vector<uintSeqLen> unhandled_bases_in_variant,
		const vector<intSeqShift> gc_mod,
		const vector<intSeqShift> end_pos_shift,
		uintSeqLen modified_start_pos,
		const Reference &comp_ref){
	DataStats stats(NULL);
	stats.read_lengths_.at(0)[100];
	stats.read_lengths_.at(1)[100];
	stats.errors_.PrepareSimulation();
	Simulator::SimPair sim_reads;
	Surrounding surrounding_end, comp_surrounding;

	uintSeqLen cur_end_position = cur_start_position+frag_length_from-1;
	species_reference_.ReverseSurrounding( surrounding_end, ref_seq_id, cur_end_position-1 );

	for(auto frag_length = frag_length_from; frag_length < frag_length_to; ++frag_length){
		surrounding_end.UpdateReverse( species_reference_.ReferenceSequence(ref_seq_id), cur_end_position++ );
		test_->UpdateBiasModForCurrentFragmentLength( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_end_position, surrounding_end );

		auto result = frag_length - frag_length_from;
		EXPECT_EQ( unhandled_variant_id.at(result), bias_mod.unhandled_variant_id_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << " Variant position: " << bias_mod.start_variant_pos_ << std::endl;
		EXPECT_EQ( unhandled_bases_in_variant.at(result), bias_mod.unhandled_bases_in_variant_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << " Variant position: " << bias_mod.start_variant_pos_ << std::endl;
		EXPECT_EQ( gc_mod.at(result), bias_mod.gc_mod_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << " Variant position: " << bias_mod.start_variant_pos_ << std::endl;
		EXPECT_EQ( end_pos_shift.at(result), bias_mod.end_pos_shift_.at(allele) ) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << " Variant position: " << bias_mod.start_variant_pos_ << std::endl;

		comp_ref.ReverseSurrounding( comp_surrounding, ref_seq_id, modified_start_pos+frag_length-1 );
		EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_end_.at(allele).sur_.at(0)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
		EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_end_.at(allele).sur_.at(1)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
		EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_end_.at(allele).sur_.at(2)) << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;

		test_->GetOrgSeq(sim_reads, allele, allele, frag_length, cur_start_position, cur_end_position, ref_seq_id, species_reference_, stats, bias_mod);

		EXPECT_TRUE( infix(comp_ref.ReferenceSequence(ref_seq_id), modified_start_pos, modified_start_pos+frag_length) == prefix(sim_reads.org_seq_.at(allele), frag_length) ) << prefix(sim_reads.org_seq_.at(allele), frag_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
		EXPECT_TRUE( ReverseComplementorDna(infix(comp_ref.ReferenceSequence(ref_seq_id), modified_start_pos, modified_start_pos+frag_length)) == prefix(sim_reads.org_seq_.at(!allele), frag_length) ) << prefix(sim_reads.org_seq_.at(!allele), frag_length) << std::endl << "Start position: " << cur_start_position << " Fragment length: " << frag_length << std::endl;
	}
}

void SimulatorTest::TestVariationInSimulateFromGivenBlock(){
	species_reference_.num_alleles_ = 2;

	// samtools faidx ../test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:1001-1020
	// GTTGCGAGATTTGGACGGAC
	species_reference_.variants_.clear();
	species_reference_.variants_.resize(1);
	species_reference_.variants_.at(0).clear();
	uintAlleleBitArray present_in_alleles = 2; //[0:no, 1:yes]
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

	Simulator::VariantBiasVarModifiers bias_mod(1, 2);
	uintRefSeqId ref_seq_id = 0;
	Surrounding surrounding_start;

	// Position 1003
	uintSeqLen cur_start_position = 1003;
	bias_mod.first_variant_id_ = 2; // Insertion at position 2 has already been processed

	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2));
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1004 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 13,
			{{2,2},{3,2},{3,2},{3,3},{3,3},{4,3},{4,3},{4,4},{5,4},{6,5},{6,6},{6,6}},
			{{0,0},{0,2},{0,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
			{{0,0},{0,-1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,1},{0,0},{0,0}},
			{{0,0},{0,0},{0,-1},{0,-2},{0,-2},{0,-2},{0,-2},{0,-1},{0,-1},{0,-1},{0,0},{0,0}},
			{cur_start_position, 1004}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1004 + 0
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, ++cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(0).sur_.at(0)) << " (Result)";
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(0).sur_.at(1)) << " (Result)";
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(0).sur_.at(2)) << " (Result)";
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1005 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 12,
			{{3,2},{3,2},{3,3},{3,3},{4,3},{4,3},{4,4},{5,4},{6,5},{6,6},{6,6}},
			{{0,2},{0,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
			{{0,-1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,1},{0,0},{0,0}},
			{{0,0},{0,-1},{0,-2},{0,-2},{0,-2},{0,-2},{0,-1},{0,-1},{0,-1},{0,0},{0,0}},
			{cur_start_position, 1005}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(1, bias_mod.start_variant_pos_);

	// Position 1004 + 1
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2));
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1006 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 11, 1,
			{3,3,3,3,3,4,4,5,6,6},
			{0,0,0,0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0,1,0,0},
			{0,-1,-1,-1,-1,0,0,0,1,1},
			1006, test_ref);

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(2, bias_mod.first_variant_id_);
	EXPECT_EQ(2, bias_mod.start_variant_pos_);

	// Position 1004 + 2
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2));
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1007 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 10, 1,
			{3,3,3,3,4,4,5,6,6},
			{0,0,0,0,0,0,0,0,0},
			{-1,-1,-1,-1,-1,-1,0,-1,-1},
			{0,0,0,0,1,1,1,2,2},
			1007, test_ref);

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(3, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1008
	cur_start_position = 1008;
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2));

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 8, 0,
			{4,4,4,5,6,6,6},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			cur_start_position, species_reference_);

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(4, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);

	// Position 1011
	cur_start_position = 1011;
	species_reference_.ForwardSurrounding( surrounding_start, ref_seq_id, cur_start_position );
	test_->PrepareBiasModForCurrentStartPos( bias_mod, ref_seq_id, species_reference_, cur_start_position, cur_start_position, surrounding_start ); // cur_end_position = cur_start_position -> Fragment length = 1
	EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.mod_surrounding_start_.at(0).sur_.at(0));
	EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.mod_surrounding_start_.at(0).sur_.at(1));
	EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.mod_surrounding_start_.at(0).sur_.at(2));
	test_ref.ForwardSurrounding( comp_surrounding, ref_seq_id, 1013 );
	EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(0)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(1)) << " (Result)";
	EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl << bitset<20>(bias_mod.mod_surrounding_start_.at(1).sur_.at(2)) << " (Result)";

	TestVariationInInnerLoopOfSimulateFromGivenBlock( bias_mod, ref_seq_id, cur_start_position, 1, 5,
			{{5,5},{6,6},{6,6},{6,6}},
			{{0,0},{0,0},{0,0},{0,0}},
			{{0,1},{0,0},{0,0},{0,0}},
			{{0,0},{0,1},{0,1},{0,1}},
			{cur_start_position, 1013}, {&species_reference_, &test_ref});

	test_->CheckForInsertedBasesToStartFrom( bias_mod, 0, cur_start_position, species_reference_ );
	EXPECT_EQ(5, bias_mod.first_variant_id_);
	EXPECT_EQ(0, bias_mod.start_variant_pos_);
}

namespace reseq{
	TEST_F(SimulatorTest, BasicFunctonality){
		CreateTestObject();

		TestCoverageConversion();
	}

	TEST_F(SimulatorTest, Variants){
		CreateTestObject();
		LoadReference("ecoli-GCF_000005845.2_ASM584v2_genomic.fa");

		TestVariationInSimulateFromGivenBlock();
	}
}
