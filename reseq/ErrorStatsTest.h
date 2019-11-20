#ifndef ERRORSTATSTEST_H
#define ERRORSTATSTEST_H
#include "ErrorStats.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class ErrorStatsTest : public BasicTestClass{
	protected:
		ErrorStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

	public:
		ErrorStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const ErrorStats &test, const char *context);
		static void TestDuplicates(const ErrorStats &test);
		static void TestVariants(const ErrorStats &test);
		static void TestAdapters(const ErrorStats &test, const char *context, bool bwa=false);
	};
}

#endif // ERRORSTATSTEST_H
