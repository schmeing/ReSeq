#ifndef FRAGMENTDUPLICATIONSTATSTEST_H
#define FRAGMENTDUPLICATIONSTATSTEST_H
#include "FragmentDuplicationStats.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class FragmentDuplicationStatsTest : public BasicTestClassWithReference{
	protected:
		FragmentDuplicationStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

	public:
		FragmentDuplicationStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const FragmentDuplicationStats &test, const char *context);
		static void TestDuplicates(const FragmentDuplicationStats &test);
		static void TestCrossDuplicates(const FragmentDuplicationStats &test);
	};
}

#endif // FRAGMENTDUPLICATIONSTATSTEST_H
