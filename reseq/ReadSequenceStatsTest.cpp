#include "ReadSequenceStatsTest.h"
using reseq::ReadSequenceStatsTest;

void ReadSequenceStatsTest::Register() {
    // Guarantees that library is included
}

void ReadSequenceStatsTest::CreateTestObject() {
    test_ = std::make_unique<ReadSequenceStats>();
    ASSERT_TRUE(test_) << "Could not allocate memory for ReadSequenceStats object\n";
}

void ReadSequenceStatsTest::DeleteTestObject() {
    test_.reset();
}

void ReadSequenceStatsTest::TearDown() {
    BasicTestClass::TearDown();
    DeleteTestObject();
}

namespace reseq {

TEST_F(ReadSequenceStatsTest, PrepareAndFinalize) {
    CreateTestObject();
    test_->PrepareAccumulators(60, 100);
    test_->IncrementMappingQuality(0, 30);
    test_->IncrementMappingQuality(0, 30);
    test_->IncrementMappingQuality(1, 20);
    test_->Finalize();
    EXPECT_EQ(2, test_->ProperPairMappingQuality()[30]) << "ProperPairMappingQuality[30] should be 2";
    EXPECT_EQ(1, test_->ImproperPairMappingQuality()[20]) << "ImproperPairMappingQuality[20] should be 1";
}

TEST_F(ReadSequenceStatsTest, ShrinkPreservesData) {
    CreateTestObject();
    test_->PrepareAccumulators(60, 100);
    test_->IncrementMappingQuality(0, 30);
    test_->Finalize();
    auto size_before = test_->ProperPairMappingQuality().size();
    ASSERT_GT(size_before, 0) << "No data to test shrink";
    test_->Shrink();
    EXPECT_GT(test_->ProperPairMappingQuality().size(), 0) << "Data lost after Shrink";
    EXPECT_EQ(1, test_->ProperPairMappingQuality()[30]) << "ProperPairMappingQuality[30] changed after Shrink";
}

TEST_F(ReadSequenceStatsTest, EmptyFinalize) {
    CreateTestObject();
    test_->PrepareAccumulators(60, 100);
    test_->Finalize();
    EXPECT_TRUE(test_->ProperPairMappingQuality().empty())
        << "ProperPairMappingQuality should be empty when no data added";
    EXPECT_TRUE(test_->ImproperPairMappingQuality().empty())
        << "ImproperPairMappingQuality should be empty when no data added";
    EXPECT_TRUE(test_->SingleReadMappingQuality().empty())
        << "SingleReadMappingQuality should be empty when no data added";
}

} // namespace reseq
