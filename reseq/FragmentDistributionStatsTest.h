#ifndef FRAGMENTDISTRIBUIONSTATSTEST_H
#define FRAGMENTDISTRIBUIONSTATSTEST_H
#include "FragmentDistributionStats.h"

#include <array>
#include <mutex>
#include <random>
#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class FragmentDistributionStatsTest : public BasicTestClassWithReference{
	public: 
		static void Register();

	protected:
		FragmentDistributionStats *test_;

		void CreateTestObject(const Reference *ref);
		void DeleteTestObject();

		static void TestOutskirtContent(const FragmentDistributionStats &test, uint16_t template_segment, uint32_t at_pos, uint64_t cont_a, uint64_t cont_c, uint64_t cont_g, uint64_t cont_t, const char * context );

		virtual void TearDown();

		void CreateCoverageDataHelper(uint8_t gc_perc, const std::array<int32_t, Reference::num_surrounding_blocks_> &start_sur, const std::array<int32_t, Reference::num_surrounding_blocks_> &end_sur, uint32_t start_pos, uint32_t fragment_length, std::mt19937_64 &rgen);
		void CreateCoverageData(uint32_t fragment_length);

	public:
		FragmentDistributionStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const FragmentDistributionStats &test, const char *context);
		static void TestDuplicates(const FragmentDistributionStats &test);
		static void TestCrossDuplicates(const FragmentDistributionStats &test);
		static void TestCoverage(const FragmentDistributionStats &test);
		static void TestAdapters(const FragmentDistributionStats &test);
		static void BiasCalculationThread(FragmentDistributionStats &test, const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, std::mutex &print_mutex);
	};
}

#endif // FRAGMENTDISTRIBUIONSTATSTEST_H
