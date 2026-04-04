#include "FragmentDuplicationStatsTest.h"
using reseq::FragmentDuplicationStatsTest;

#include <string>
using std::string;
#include <vector>
using std::vector;

void FragmentDuplicationStatsTest::Register() {
    // Guarantees that library is included
}

void FragmentDuplicationStatsTest::CreateTestObject() {
    test_ = std::make_unique<FragmentDuplicationStats>();
    ASSERT_TRUE(test_) << "Could not allocate memory for FragmentDuplicationStats object\n";
}

void FragmentDuplicationStatsTest::DeleteTestObject() {
    test_.reset();
}

void FragmentDuplicationStatsTest::TearDown() {
    BasicTestClass::TearDown();
    DeleteTestObject();
}

void FragmentDuplicationStatsTest::TestSrr490124Equality(const FragmentDuplicationStats& test, const char* context) {
    EXPECT_EQ(8, test.duplication_number_.size())
        << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.duplication_number_[1]) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
    EXPECT_EQ(8, test.duplication_number_[8]) << "SRR490124-4pairs duplication_number_ wrong for " << context << '\n';
}

void FragmentDuplicationStatsTest::TestDuplicates(const FragmentDuplicationStats& test) {
    // For the different strands (numbers have to be summed):
    TestVectEquality({1, {3, 8, 3}}, test.duplication_number_, "duplicates test", "duplication_number_",
                     " not correct for ");
}

void FragmentDuplicationStatsTest::TestCrossDuplicates(const FragmentDuplicationStats& test) {
    EXPECT_EQ(0, test.duplication_number_.size()) << "Cross duplicates do not count to duplications anymore\n";
}

namespace reseq {

TEST_F(FragmentDuplicationStatsTest, AddDuplication) {
    CreateTestObject();
    test_->PrepareTmpDuplicationVector();
    test_->AddDuplication(1);
    test_->AddDuplication(1);
    test_->AddDuplication(5);
    test_->AddDuplication(FragmentDuplicationStats::kMaxDuplication + 10);
    test_->FinalizeDuplicationVector();
    // FinalizeDuplicationVector multiplies counts by dup number (#Sites -> #Fragments)
    EXPECT_EQ(2, test_->DuplicationNumber()[1]) << "Duplication count for dup=1 incorrect (2 sites * 1)";
    EXPECT_EQ(5, test_->DuplicationNumber()[5]) << "Duplication count for dup=5 incorrect (1 site * 5)";
    EXPECT_EQ(FragmentDuplicationStats::kMaxDuplication + 1,
              test_->DuplicationNumber()[FragmentDuplicationStats::kMaxDuplication + 1])
        << "Overflow bin count incorrect (1 site * (kMaxDuplication+1))";
}

TEST_F(FragmentDuplicationStatsTest, FinalizeDuplicationVector) {
    CreateTestObject();
    test_->PrepareTmpDuplicationVector();
    test_->AddDuplication(3);
    test_->AddDuplication(3);
    test_->AddDuplication(3);
    test_->FinalizeDuplicationVector();
    const auto& dn = test_->DuplicationNumber();
    EXPECT_FALSE(dn.empty()) << "DuplicationNumber empty after finalization";
    EXPECT_EQ(9, dn[3]) << "DuplicationNumber[3] should be 9 (3 sites * 3)";
}

TEST_F(FragmentDuplicationStatsTest, EmptyFinalize) {
    CreateTestObject();
    test_->PrepareTmpDuplicationVector();
    test_->FinalizeDuplicationVector();
    EXPECT_TRUE(test_->DuplicationNumber().empty()) << "DuplicationNumber should be empty when no data added";
}

} // namespace reseq
