#ifndef READSEQUENCESTATSTEST_H
#define READSEQUENCESTATSTEST_H

#include "ReadSequenceStats.h"

#include <cstdint>
#include <memory>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq {
class ReadSequenceStatsTest : public BasicTestClass {
  public:
    static void Register();

  protected:
    std::unique_ptr<ReadSequenceStats> test_;

    void CreateTestObject();
    void DeleteTestObject();

    virtual void TearDown();

  public:
    ReadSequenceStatsTest() {}
};
} // namespace reseq

#endif // READSEQUENCESTATSTEST_H
