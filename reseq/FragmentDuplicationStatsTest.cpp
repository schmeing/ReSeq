#include "FragmentDuplicationStatsTest.h"
using reseq::FragmentDuplicationStatsTest;

#include <string>
using std::string;
#include <vector>
using std::vector;

void FragmentDuplicationStatsTest::Register(){
	// Guarantees that library is included
}

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
	EXPECT_EQ(8, test.duplication_number_.size()) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.duplication_number_[1]) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
	EXPECT_EQ(8, test.duplication_number_[8]) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
}

void FragmentDuplicationStatsTest::TestDuplicates(const FragmentDuplicationStats &test){
	// For the different strands (numbers have to be summed):
	TestVectEquality({1,{3,8,3}}, test.duplication_number_, "duplicates test", "duplication_number_", " not correct for ");
}

void FragmentDuplicationStatsTest::TestCrossDuplicates(const FragmentDuplicationStats &test){
	EXPECT_EQ(0, test.duplication_number_.size() ) << "Cross duplicates do not count to duplications anymore\n";
}
