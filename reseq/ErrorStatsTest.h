#ifndef ERRORSTATSTEST_H
#define ERRORSTATSTEST_H
#include "ErrorStats.h"

#include <cstdint>
#include <memory>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq {
class ErrorStatsTest : public BasicTestClass {
  protected:
    std::unique_ptr<ErrorStats> test_;

    void CreateTestObject();
    void DeleteTestObject();

    virtual void TearDown();

  public:
    ErrorStatsTest() {}

    static void TestSrr490124Equality(const ErrorStats& test, const char* context);
    static void TestDuplicates(const ErrorStats& test);
    static void TestVariants(const ErrorStats& test);
    static void TestAdapters(const ErrorStats& test, const char* context, bool bwa = false);
};
} // namespace reseq

#endif // ERRORSTATSTEST_H
