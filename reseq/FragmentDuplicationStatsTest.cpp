#include "FragmentDuplicationStatsTest.h"
using reseq::FragmentDuplicationStatsTest;

#include <string>
using std::string;
#include <vector>
using std::vector;

void FragmentDuplicationStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new FragmentDuplicationStats ) << "Could not allocate memory for FragmentDuplicationStats object\n";
}

void FragmentDuplicationStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void FragmentDuplicationStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void FragmentDuplicationStatsTest::TestSrr490124Equality(const FragmentDuplicationStats &test, const char *context){
	EXPECT_EQ(1, test.duplication_number_.size()) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.duplication_number_[1]) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
}

void FragmentDuplicationStatsTest::TestDuplicates(const FragmentDuplicationStats &test){
	// For the different strands (numbers have to be summed):
	// samtools view -q 10 -f 67 ecoli-duplicates.bam | awk '(0 != substr($1,12,1) && ($4>$8 || $4==$8 && int($2%256/128))){len=0; for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e || "D"==e){len+=num}; num=0}}; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $8 "-" $4+len-1)}' | seqtk seq | awk '(0==NR%2){len=length($0); print int(gsub(/[GC]/,"",$0)*100/len+0.5), len}' | sort -k1,1n -k2,2n | uniq -c
	// samtools view -q 10 -f 131 ecoli-duplicates.bam | awk '(0 != substr($1,12,1) && ($4>$8 || $4==$8 && int($2%256/128))){len=0; for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e || "D"==e){len+=num}; num=0}}; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $8 "-" $4+len-1)}' | seqtk seq | awk '(0==NR%2){len=length($0); print int(gsub(/[GC]/,"",$0)*100/len+0.5), len}' | sort -k1,1n -k2,2n | uniq -c
	EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[36][50].size() ) << "duplication_number_by_gc_insert_length_ not shrunken correctly\n";
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[33][99][2] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[36][50][3] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[36][150][2] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[39][151][1] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[39][151][2] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[44][201][2] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[45][301][1] );
    EXPECT_EQ( 1, test.duplication_number_by_gc_insert_length_[51][151][1] );

	TestVectEquality({1,{3,4,1}}, test.duplication_number_, "duplicates test", "duplication_number_", " not correct for ");
}

void FragmentDuplicationStatsTest::TestCrossDuplicates(const FragmentDuplicationStats &test){
	EXPECT_EQ(0, test.duplication_number_by_gc_insert_length_.size() ) << "Cross duplicates do not count to duplications anymore\n";
	EXPECT_EQ(0, test.duplication_number_.size() ) << "Cross duplicates do not count to duplications anymore\n";
}

namespace reseq{
	TEST_F(FragmentDuplicationStatsTest, DispersionCalculation){
		CreateTestObject();

		test_->CalculateDispersionPlot();
		EXPECT_EQ(0, test_->mean_list_.size()) << "Empty vector returns wrong results\n";
		EXPECT_EQ(0, test_->dispersion_list_.size()) << "Empty vector returns wrong results\n";

		// Not valid: No count of 2 (mean>variance)
		test_->duplication_number_by_gc_insert_length_[1][1][0] = 3;
		test_->duplication_number_by_gc_insert_length_[1][1][1] = 2;
		// Valid
		test_->duplication_number_by_gc_insert_length_[1][3][0] = 17;
		test_->duplication_number_by_gc_insert_length_[1][3][1] = 2;
		test_->duplication_number_by_gc_insert_length_[1][3][2] = 1;
		// Not valid: mean>variance
		test_->duplication_number_by_gc_insert_length_[1][5][0] = 2;
		test_->duplication_number_by_gc_insert_length_[1][5][1] = 2;
		test_->duplication_number_by_gc_insert_length_[1][5][2] = 1;

		test_->CalculateDispersionPlot();
		EXPECT_EQ(1, test_->dispersion_list_.size()) << "Shrinkage of vector did not work properly or wrong values are inserted\n";
		EXPECT_DOUBLE_EQ(0.2, test_->mean_list_[0]) << "Mean wrong\n";
		EXPECT_DOUBLE_EQ(2.0/3, test_->dispersion_list_[0]) << "Dispersion wrong\n";
	}
}
