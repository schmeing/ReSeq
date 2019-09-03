#ifndef COVERAGESTATSTEST_H
#define COVERAGESTATSTEST_H
#include "CoverageStats.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class CoverageStatsTest : public BasicTestClass{
	protected:
		CoverageStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

	public:
		CoverageStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const CoverageStats &test, const char *context);
		static void TestDuplicates(const CoverageStats &test);
		static void TestVariants(const CoverageStats &test);
		static void TestCrossDuplicates(const CoverageStats &test);
		static void TestCoverage(const CoverageStats &test);
	};
}

#endif // COVERAGESTATSTEST_H
