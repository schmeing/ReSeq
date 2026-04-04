#include <atomic>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

namespace reseq {
std::atomic<uint16_t> kVerbosityLevel{2};
bool kNoDebugOutput = false;
} // namespace reseq
#include "logging.hpp"

#include "AdapterStatsTest.h"
#include "CoverageStatsTest.h"
#include "DataStatsTest.h"
#include "ErrorStatsTest.h"
#include "FragmentDistributionStatsTest.h"
#include "FragmentDuplicationStatsTest.h"
#include "ProbabilityEstimatesTest.h"
#include "QualityStatsTest.h"
#include "ReadSequenceStatsTest.h"
#include "ReferenceTest.h"
#include "RegressionTest.h"
#include "SeqQualityStatsTest.h"
#include "SimulatorTest.h"
#include "SurroundingTest.h"
#include "TileStatsTest.h"
#include "utilitiesTest.h"
#include "VectTest.h"

#include "BasicTestClass.hpp"

// Definitions so that referencing a const static is valid
#include <seqan/sequence.h>
seqan::FunctorComplement<seqan::Dna5> reseq::utilities::Complement::Dna5;
seqan::FunctorComplement<seqan::Dna> reseq::utilities::Complement::Dna;

int main(int argc, char** argv) {
    std::string test_dir;
    if (!reseq::BasicTestClass::GetTestDir(test_dir)) {
        std::cerr << "ERROR: Cannot find test data directory. Aborting." << std::endl;
        return 1;
    }

    // Register all test classes (mirrors main.cpp registration order)
    reseq::AdapterStatsTest::Register();
    reseq::DataStatsTest::Register();
    reseq::FragmentDistributionStatsTest::Register(1);
    reseq::FragmentDuplicationStatsTest::Register();
    reseq::ProbabilityEstimatesTest::Register();
    reseq::ReadSequenceStatsTest::Register();
    reseq::ReferenceTest::Register();
    reseq::SeqQualityStatsTest::Register();
    reseq::SimulatorTest::Register();
    reseq::SurroundingTest::Register();
    reseq::TileStatsTest::Register();
    reseq::VectTest::Register();
    reseq::utilitiesTest::Register();
    reseq::RegressionTest::Register();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
