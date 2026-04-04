#ifndef COVERAGESTATSTEST_H
#define COVERAGESTATSTEST_H
#include "CoverageStats.h"

#include <cstdint>
#include <memory>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq {
class CoverageStatsTest : public BasicTestClassWithReference {
  public:
    static void Register();

  protected:
    std::unique_ptr<CoverageStats> test_;

    void CreateTestObject();
    void DeleteTestObject();

    virtual void TearDown();

    void TestNonSystematicErrorRate();
    CoverageStats::CoverageBlock* BootstrapFirstBlock(CoverageStats& cs, uintRefSeqId seq_id, uintSeqLen start_pos);
    void TestDequeLifecycle();
    void TestFindByIndex();
    void TestCleanupRecyclesIndices();
    void TestConcurrentCoverageIncrement();

  public:
    CoverageStatsTest() {}

    static void TestSrr490124Equality(const CoverageStats& test, const char* context);
    static void TestDuplicates(const CoverageStats& test);
    static void TestVariants(const CoverageStats& test);
    static void TestCrossDuplicates(const CoverageStats& test);
    static void TestCoverage(const CoverageStats& test);
};
} // namespace reseq

#endif // COVERAGESTATSTEST_H
