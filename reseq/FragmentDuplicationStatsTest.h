#ifndef FRAGMENTDUPLICATIONSTATSTEST_H
#define FRAGMENTDUPLICATIONSTATSTEST_H
#include "FragmentDuplicationStats.h"

#include <cstdint>
#include <memory>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq {
class FragmentDuplicationStatsTest : public BasicTestClassWithReference {
  public:
    static void Register();

  protected:
    std::unique_ptr<FragmentDuplicationStats> test_;

    void CreateTestObject();
    void DeleteTestObject();

    virtual void TearDown();

  public:
    FragmentDuplicationStatsTest() {}

    static void TestSrr490124Equality(const FragmentDuplicationStats& test, const char* context);
    static void TestDuplicates(const FragmentDuplicationStats& test);
    static void TestCrossDuplicates(const FragmentDuplicationStats& test);
};
} // namespace reseq

#endif // FRAGMENTDUPLICATIONSTATSTEST_H
