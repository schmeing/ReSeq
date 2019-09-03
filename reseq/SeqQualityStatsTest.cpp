#include "SeqQualityStats.hpp"

#include <stdint.h>

#include "gtest/gtest.h"

namespace reseq{
	TEST(SeqQualityStatsTest, EmptyClass){
		SeqQualityStats<uint32_t> test;
		test.Calculate();
		EXPECT_EQ(0, test.mean_);
		EXPECT_EQ(0, test.minimum_);
		EXPECT_EQ(0, test.first_quartile_);
		EXPECT_EQ(0, test.median_);
		EXPECT_EQ(0, test.third_quartile_);
		EXPECT_EQ(0, test.maximum_);
		EXPECT_TRUE(0 == test) << "Equality operator for class 'SeqQualityStats' had wrong result for empty qualities_.\n";
	}

	TEST(SeqQualityStatsTest, Values5){
		SeqQualityStats<uint16_t> test;
		test[5] = 4;
		test[10] = 5;
		test[12] = 1;
		test[15] = 4;
		test[34] = 5;
		test.Calculate();
		EXPECT_EQ(16, test.mean_);
		EXPECT_EQ(5, test.minimum_);
		EXPECT_EQ(10, test.first_quartile_);
		EXPECT_EQ(12, test.median_);
		EXPECT_EQ(34, test.third_quartile_);
		EXPECT_EQ(34, test.maximum_);
		EXPECT_TRUE(35 == test) << "Equality operator for class 'SeqQualityStats' had wrong result.\n";
	}

	TEST(SeqQualityStatsTest, Value1){
		SeqQualityStats<uint64_t> test;
		test[7] = 1;
		test.Calculate();
		EXPECT_EQ(7, test.mean_);
		EXPECT_EQ(7, test.minimum_);
		EXPECT_EQ(7, test.first_quartile_);
		EXPECT_EQ(7, test.median_);
		EXPECT_EQ(7, test.third_quartile_);
		EXPECT_EQ(7, test.maximum_);
	}

	TEST(SeqQualityStatsTest, Value6){
		SeqQualityStats<uint32_t> test;
		test[7] = 1;
		test[15] = 1;
		test[36] = 1;
		test[39] = 1;
		test[40] = 1;
		test[41] = 1;
		test.Calculate();
		EXPECT_EQ(30, test.mean_);
		EXPECT_EQ(7, test.minimum_);
		EXPECT_EQ(15, test.first_quartile_);
		EXPECT_EQ(39, test.median_);
		EXPECT_EQ(40, test.third_quartile_);
		EXPECT_EQ(41, test.maximum_);
	}

	TEST(SeqQualityStatsTest, Value8){
		SeqQualityStats<uint32_t> test;
		test[7] = 1;
		test[15] = 1;
		test[36] = 1;
		test[39] = 1;
		test[40] = 1;
		test[41] = 1;
		test[6] = 1;
		test[42] = 1;
		test.Calculate();
		EXPECT_EQ(28, test.mean_);
		EXPECT_EQ(6, test.minimum_);
		EXPECT_EQ(7, test.first_quartile_);
		EXPECT_EQ(39, test.median_);
		EXPECT_EQ(41, test.third_quartile_);
		EXPECT_EQ(42, test.maximum_);
	}

	TEST(SeqQualityStatsTest, Value9){
		SeqQualityStats<uint32_t> test;
		test[7] = 1;
		test[15] = 1;
		test[36] = 1;
		test[39] = 1;
		test[40] = 1;
		test[41] = 1;
		test[6] = 1;
		test[42] = 1;
		test[43] = 1;
		test.Calculate();
		EXPECT_EQ(30, test.mean_);
		EXPECT_EQ(6, test.minimum_);
		EXPECT_EQ(7, test.first_quartile_);
		EXPECT_EQ(39, test.median_);
		EXPECT_EQ(42, test.third_quartile_);
		EXPECT_EQ(43, test.maximum_);
	}

	TEST(SeqQualityStatsTest, Value11){
		SeqQualityStats<uint32_t> test;
		test[7] = 1;
		test[15] = 1;
		test[36] = 1;
		test[39] = 1;
		test[40] = 1;
		test[41] = 1;
		test[6] = 1;
		test[42] = 1;
		test[43] = 1;
		test[47] = 1;
		test[49] = 1;
		test.Calculate();
		EXPECT_EQ(33, test.mean_);
		EXPECT_EQ(6, test.minimum_);
		EXPECT_EQ(15, test.first_quartile_);
		EXPECT_EQ(40, test.median_);
		EXPECT_EQ(43, test.third_quartile_);
		EXPECT_EQ(49, test.maximum_);
	}

	TEST(SeqQualityStatsTest, Probability){
		SeqQualityStats<uint32_t> test;
		test[10] = 11;
		test[20] = 20;
		test[30] = 30;
		test[40] = 39;
		test.CalculateProbabilityMean();
		EXPECT_EQ(19, test.probability_mean_);
	}
}
