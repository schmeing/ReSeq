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
		static void Register(uintNumThreads num_threads);

	protected:
		FragmentDistributionStats *test_;

		void CreateTestObject(const Reference *ref);
		void DeleteTestObject();

		static void TestOutskirtContent(const FragmentDistributionStats &test, uintTempSeq template_segment, uintSeqLen at_pos, uintNucCount cont_a, uintNucCount cont_c, uintNucCount cont_g, uintNucCount cont_t, const char * context );

		virtual void TearDown();

		void CreateCoverageDataHelper(uintPercent gc_perc, uintSeqLen start_pos, uintSeqLen fragment_length, std::mt19937_64 &rgen);
		void CreateCoverageData(uintSeqLen fragment_length);

		template<size_t N> void CheckDrawnCounts(double bias, std::array<double, N> thresholds, const Surrounding &start_sur, const Surrounding &end_sur);

	public:
		FragmentDistributionStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const FragmentDistributionStats &test, const char *context);
		static void TestDuplicates(const FragmentDistributionStats &test);
		static void TestCrossDuplicates(const FragmentDistributionStats &test);
		static void TestCoverage(const FragmentDistributionStats &test);
		static void TestAdapters(const FragmentDistributionStats &test, const char *context, bool bwa=false);
		static void BiasCalculationThread(FragmentDistributionStats &test, const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, std::mutex &print_mutex);
		void TestBiasCalculationVectorsPreprocessing();
		void TestBiasCalculationVectorsNormalizations();
		void TestBiasCalculationVectorsSpline();
		void TestBiasCalculationVectorsLikelihoods();
		void TestBiasCalculation();
		void TestUniformBias();
		void TestDrawCounts();
		void TestRefSeqSplitting();
		void TestRefBinProcessing();
	};
}

#endif // FRAGMENTDISTRIBUIONSTATSTEST_H
