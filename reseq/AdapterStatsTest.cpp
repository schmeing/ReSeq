#include "AdapterStatsTest.h"
using reseq::AdapterStatsTest;

#include <string>
using std::string;

#include <seqan/sequence.h>
using seqan::CharString;

#include "utilities.hpp"
using reseq::utilities::GetReSeqDir;

void AdapterStatsTest::Register(){
	// Guarantees that library is included
}

void AdapterStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new AdapterStats() ) << "Could not allocate memory for DataStats object\n";
}

void AdapterStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void AdapterStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void AdapterStatsTest::LoadAdapters(){
	string adapter_dir;
	ASSERT_TRUE( GetReSeqDir(adapter_dir, "adapters", "TruSeq_single.fa") );
	EXPECT_TRUE( test_->LoadAdapters( (adapter_dir+"TruSeq_single.fa").c_str(), (adapter_dir+"TruSeq_single.mat").c_str() ) ) << "Problem loading the default adapters";
	test_->PrepareAdapters(33, 101);
}

void AdapterStatsTest::TestLoading(const AdapterStats &test){
	// The tests below were done manually and are not updated yet with verification commands
	ASSERT_EQ(48, test.seqs_.at(0).size() );
	EXPECT_EQ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", string(toCString(static_cast<CharString>(test.seqs_.at(0).at(0)))) );
	EXPECT_EQ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG", string(toCString(static_cast<CharString>(test.seqs_.at(0).at(23)))) );
	ASSERT_EQ(48, test.names_.at(0).size() );
	EXPECT_EQ("TruSeq_Adapter_Index_5_f", test.names_.at(0).at(0) );
	EXPECT_EQ("TruSeq_Adapter_Index_3_f", test.names_.at(0).at(23) );

	ASSERT_EQ(2, test.seqs_.at(1).size() );
	EXPECT_EQ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", string(toCString(static_cast<CharString>(test.seqs_.at(1).at(1)))) );
	ASSERT_EQ(2, test.names_.at(1).size() );
	EXPECT_EQ("TruSeq_Universal_Adapter_r", test.names_.at(1).at(1) );

	ASSERT_EQ(48, test.combinations_.size() );
	ASSERT_EQ(2, test.combinations_.at(0).size() );
	EXPECT_TRUE( test.combinations_.at(0).at(0) );
	ASSERT_EQ(2, test.combinations_.at(23).size() );
	EXPECT_TRUE( test.combinations_.at(23).at(0) );
}

void AdapterStatsTest::TestSumming(AdapterStats &test){
	test.counts_.resize(test.seqs_.at(0).size());
	for( auto &dim2 : test.counts_ ){
		dim2.resize(test.seqs_.at(1).size());
	}
	for( auto seg=2; seg--; ){
		test.start_cut_.at(seg).resize(test.seqs_.at(seg).size()); // Is the size SumCounts uses to resize count_sum_
	}

	// seqtk seq ../adapters/TruSeq_single.fa | awk '(NR%2==0 && NR>2)' | sort
	// cat <(seqtk seq ../adapters/TruSeq_single.fa | awk 'NR==2') <(seqtk seq -r ../adapters/TruSeq_single.fa | awk 'NR==2') | sort
	test.counts_.at(7).at(0)[35][0] = 1; // ambiguous
	test.counts_.at(7).at(0)[35][1] = 2; // ambiguous
	test.counts_.at(7).at(0)[36][0] = 4; // ambiguous
	test.counts_.at(7).at(0)[36][1] = 8;
	test.counts_.at(7).at(1)[35][0] = 16; // ambiguous
	test.counts_.at(7).at(1)[35][1] = 32; // ambiguous
	test.counts_.at(7).at(1)[36][0] = 64; // ambiguous
	test.counts_.at(7).at(1)[36][1] = 128;
	test.counts_.at(8).at(0)[34][0] = 256; // ambiguous
	test.counts_.at(8).at(0)[34][1] = 512; // ambiguous
	test.counts_.at(8).at(0)[35][0] = 1024; // ambiguous
	test.counts_.at(8).at(0)[35][1] = 2048;
	test.counts_.at(8).at(1)[34][0] = 4096; // ambiguous
	test.counts_.at(8).at(1)[34][1] = 8192; // ambiguous
	test.counts_.at(8).at(1)[35][0] = 16384; // ambiguous
	test.counts_.at(8).at(1)[35][1] = 32768;

	test.SumCounts();

	EXPECT_EQ(136, test.count_sum_.at(0).at(7));
	EXPECT_EQ(34816, test.count_sum_.at(0).at(8));
	EXPECT_EQ(2056, test.count_sum_.at(1).at(0));
	EXPECT_EQ(32896, test.count_sum_.at(1).at(1));

	test.PrepareSimulation();
	EXPECT_EQ(0, test.significant_count_.at(0).at(7));
	EXPECT_EQ(34816, test.significant_count_.at(0).at(8));
	EXPECT_EQ(0, test.significant_count_.at(1).at(0));
	EXPECT_EQ(32896, test.significant_count_.at(1).at(1));
}

void AdapterStatsTest::TestAdapters(const AdapterStats &test, const char *context){
	TestLoading(test);

	// The tests below were done manually and are not updated yet with verification commands
	for( uintAdapterId i=0; i < test.counts_.size(); ++i  ){
		for( uintAdapterId j=0; j < test.counts_.at(i).size(); ++j  ){
			if( 8==i && 1==j ){
				EXPECT_EQ(59, test.counts_.at(i).at(j).size()) << "Combination " << test.names_.at(0).at(i) << ' ' << test.names_.at(1).at(j) << " is not detected properly for " << context;
				EXPECT_EQ(2, test.counts_.at(i).at(j)[100][100]) << "Position 0 not correctly detected for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[42][42]) << "Position 58 not correctly detected for " << context;

			}
			else if( 0==i && 1==j ){
				// Ambiguous adapter is always at first index it can be (so 0). Trying to redistribute them properly just messes stuff up, so they are simply ignored during summing.
				EXPECT_EQ(16, test.counts_.at(i).at(j).size()) << "Combination " << test.names_.at(0).at(i) << ' ' << test.names_.at(1).at(j) << " is not detected properly for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[20].size()) << "Size at position 80 not correctly detected for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[20][20]) << "Position 80 not correctly detected for " << context; // Adapter with mapping quality 8 still used for adapters
				EXPECT_EQ(1, test.counts_.at(i).at(j)[19].size()) << "Size at position 81 not correctly detected for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[19][19]) << "Position 81 not correctly detected for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[5].size()) << "Size at position 95 not correctly detected for " << context;
				EXPECT_EQ(1, test.counts_.at(i).at(j)[5][5]) << "Position 95 not correctly detected for " << context;
			}
			else{
				EXPECT_EQ(0, test.counts_.at(i).at(j).size()) << "Combination " << test.names_.at(0).at(i) << ' ' << test.names_.at(1).at(j) << " is detected without being in the data for " << context;
			}
		}

		if( 8==i ){
			EXPECT_EQ(2, test.start_cut_.at(0).at(i).size()) << test.names_.at(0).at(i) << " is not detected for " << context;
			EXPECT_EQ(1, test.start_cut_.at(0).at(i)[0]) << " 0-cut is not detected for " << context;
			EXPECT_EQ(1, test.start_cut_.at(0).at(i)[1]) << " 1-cut is not detected for " << context;
		}
		else{
			EXPECT_EQ(0, test.start_cut_.at(0).at(i).size()) << test.names_.at(0).at(i) << " is detected without being in the data for " << context;
		}
	}
	EXPECT_EQ(2, test.start_cut_.at(1).at(1).size()) << test.names_.at(1).at(0) << " is not detected for " << context;
	EXPECT_EQ(1, test.start_cut_.at(1).at(1)[0]) << " 0-cut is not detected for " << context;
	EXPECT_EQ(1, test.start_cut_.at(1).at(1)[1]) << " 1-cut is not detected for " << context;

	// Ambiguous adapter ignored, so 3 instead of 4
	EXPECT_EQ(3, test.count_sum_.at(0).at(8)) << "Distribution did not work properly for " << test.names_.at(0).at(8) << " for " << context;
	EXPECT_EQ(3, test.count_sum_.at(1).at(1)) << "Distribution did not work properly for " << test.names_.at(1).at(0) << " for " << context;

	EXPECT_EQ(13, test.polya_tail_length_.size() ) << "for " << context;
	EXPECT_EQ(1, test.polya_tail_length_[1] ) << "for " << context;
	EXPECT_EQ(1, test.polya_tail_length_[2] ) << "for " << context;
	EXPECT_EQ(1, test.polya_tail_length_[9] ) << "for " << context;
	EXPECT_EQ(1, test.polya_tail_length_[13] ) << "for " << context;
	EXPECT_EQ(4, SumVect(test.polya_tail_length_) ) << "for " << context;

	EXPECT_EQ(50, test.overrun_bases_.at(0) ) << "for " << context;
	EXPECT_EQ(38, test.overrun_bases_.at(1) ) << "for " << context;
	EXPECT_EQ(31, test.overrun_bases_.at(2) ) << "for " << context;
	EXPECT_EQ(14, test.overrun_bases_.at(3) ) << "for " << context;
	EXPECT_EQ(0, test.overrun_bases_.at(4) ) << "for " << context;
}

void AdapterStatsTest::TestNexteraAdapters(const AdapterStats &test){
	// Adapters are internally sorted by sequence
	// seqtk seq ../adapters/Nextera_XTv2.fa | awk '(NR%2==1){name=$0}(NR%2==0){print $0, name}' | head -n24 | sort | awk '{print NR-1, $0}'
	ASSERT_EQ(48, test.names_.at(0).size() );
	EXPECT_EQ("Index1Read_N701_f", test.names_.at(0).at(21) );
	EXPECT_EQ("Index1Read_N729_f", test.names_.at(0).at(9) );
	EXPECT_EQ("Index1Read_N701_r", test.names_.at(0).at(41) );
	EXPECT_EQ("Index1Read_N729_r", test.names_.at(0).at(46) );

	// seqtk seq -r ../adapters/Nextera_XTv2.fa | awk '(NR%2==1){name=$0}(NR%2==0){print $0, name}' | tail -n16 | sort | awk '{print NR-1, $0}'
	ASSERT_EQ(32, test.names_.at(1).size() );
	EXPECT_EQ("Index2Read_S502_f", test.names_.at(1).at(7) );
	EXPECT_EQ("Index2Read_S522_f", test.names_.at(1).at(14) );
	EXPECT_EQ("Index2Read_S502_r", test.names_.at(1).at(20) );
	EXPECT_EQ("Index2Read_S522_r", test.names_.at(1).at(30) );

	ASSERT_EQ(48, test.combinations_.size() );
	for(auto i=test.combinations_.size(); i--; ){
		ASSERT_EQ(32, test.combinations_.at(i).size() );
		for(auto j=test.combinations_.at(i).size(); j--; ){
			EXPECT_TRUE( test.combinations_.at(i).at(j) );
		}
	}

	// awk 's=index($0,"CCGAGCCCACGAGAC") {print length($0)-s+1}' ecoli-S1L001-adapters-R1.fq
	// awk 's=index($0,"GACGCTGCCGACGA") {print length($0)-s+1}' ecoli-S1L001-adapters-R2.fq
	EXPECT_EQ(1, test.counts_.at(38).at(28)[103][103]);
	EXPECT_EQ(1, test.counts_.at(38).at(28)[85][85]);
	EXPECT_EQ(1, test.counts_.at(38).at(28)[43][43]);

	// None of the three starts at the beginning
	EXPECT_EQ(0, test.start_cut_.at(0).at(38)[0]);
	EXPECT_EQ(0, test.start_cut_.at(1).at(28)[0]);

	EXPECT_EQ(3, test.count_sum_.at(0).at(38)) << "Distribution did not work properly for " << test.names_.at(0).at(38) << '\n';
	EXPECT_EQ(3, test.count_sum_.at(1).at(28)) << "Distribution did not work properly for " << test.names_.at(1).at(28) << '\n';
}

namespace reseq{
	TEST_F(AdapterStatsTest, LoadingAndSumming){
		CreateTestObject();

		LoadAdapters();
		TestLoading(*test_);

		TestSumming(*test_);
	}
}
