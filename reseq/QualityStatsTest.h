#ifndef QUALITYSTATSTEST_H
#define QUALITYSTATSTEST_H
#include "QualityStats.h"

#include <cmath>
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
				const Vect<T> &vector,
				size_t size0,
				size_t size1,
				size_t size98,
				size_t size99,
				const char *context,
				const char *err_msg
				){
			EXPECT_EQ(size0, vector[0].size()) << err_msg << context << '\n';
			EXPECT_EQ(size1, vector[1].size()) << err_msg << context << '\n';
			EXPECT_EQ(size98, vector[98].size()) << err_msg << context << '\n';
			EXPECT_EQ(size99, vector[99].size()) << err_msg << context << '\n';
		}

		template<typename T, typename U> static void TestValueAtIndex_0_1_98_99(
				const Vect<T> &vector,
				U val0,
				U val1,
				U val98,
				U val99,
				const char *context,
				const char *err_msg
				){
			EXPECT_EQ(val0, vector[0]) << err_msg << context << '\n';
			EXPECT_EQ(val1, vector[1]) << err_msg << context << '\n';
			EXPECT_EQ(val98, vector[98]) << err_msg << context << '\n';
			EXPECT_EQ(val99, vector[99]) << err_msg << context << '\n';
		}

		template<typename U> static void TestRoundedValueAtIndex_0_1_98_99(
				const Vect<double> &vector,
				U val0,
				U val1,
				U val98,
				U val99,
				const char *context,
				const char *err_msg
				){
			EXPECT_EQ(val0, static_cast<uint64_t>(std::round(vector[0]))) << err_msg << context << '\n';
			EXPECT_EQ(val1, static_cast<uint64_t>(std::round(vector[1]))) << err_msg << context << '\n';
			EXPECT_EQ(val98, static_cast<uint64_t>(std::round(vector[98]))) << err_msg << context << '\n';
			EXPECT_EQ(val99, static_cast<uint64_t>(std::round(vector[99]))) << err_msg << context << '\n';
		}

		void TestRemoveSystematic();

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
		static void TestAdapters(const QualityStats &test, const char *context, bool bwa=false);
	};
}

#endif // QUALITYSTATSTEST_H
