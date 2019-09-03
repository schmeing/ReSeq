#ifndef QUALITYSTATSTEST_H
#define QUALITYSTATSTEST_H
#include "QualityStats.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class QualityStatsTest : public BasicTestClass{
	protected:
		QualityStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

		template<typename T> static void TestSizeAtIndex_0_1_98_99(
				T& vector,
				uint64_t size0,
				uint64_t size1,
				uint64_t size98,
				uint64_t size99,
				const char *context,
				const char *err_msg
				){
			EXPECT_EQ(size0, vector[0].size()) << err_msg << context << '\n';
			EXPECT_EQ(size1, vector[1].size()) << err_msg << context << '\n';
			EXPECT_EQ(size98, vector[98].size()) << err_msg << context << '\n';
			EXPECT_EQ(size99, vector[99].size()) << err_msg << context << '\n';
		}

		template<typename T> static void TestValueAtIndex_0_1_98_99(
				T& vector,
				uint64_t val0,
				uint64_t val1,
				uint64_t val98,
				uint64_t val99,
				const char *context,
				const char *err_msg
				){
			EXPECT_EQ(val0, vector[0]) << err_msg << context << '\n';
			EXPECT_EQ(val1, vector[1]) << err_msg << context << '\n';
			EXPECT_EQ(val98, vector[98]) << err_msg << context << '\n';
			EXPECT_EQ(val99, vector[99]) << err_msg << context << '\n';
		}

	public:
		QualityStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const QualityStats &test, const char *context);
		static void TestTiles(const QualityStats &test);
		static void TestDuplicates(const QualityStats &test);
		static void TestVariants(const QualityStats &test);
		static void TestCoverage(const QualityStats &test);
		static void TestAdapters(const QualityStats &test);
	};
}

#endif // QUALITYSTATSTEST_H
