#include "utilitiesTest.h"
using reseq::utilitiesTest;

#include <stdint.h>

#include "gtest/gtest.h"

void utilitiesTest::Register(){
	// Guarantees that library is included
}

namespace reseq{
	namespace utilities{
		TEST(utilitiesTest, Divide){
			EXPECT_EQ(0, Divide(0, 2)) << "'Divide' function from utilities does not return proper zero\n";
			EXPECT_EQ(1, Divide(17438564308265206, 17438564308265206)) << "'Divide' function from utilities does not return proper 1\n";
			EXPECT_EQ(11, Divide(73500, 7000)) << "'Divide' function from utilities does not round up properly\n";
			EXPECT_EQ(10, Divide(73400, 7000)) << "'Divide' function from utilities does not round down properly\n";
		}

		TEST(utilitiesTest, Percent){
			EXPECT_EQ(0, Percent(0, 2)) << "'Percent' function from utilities does not return proper zero\n";
			EXPECT_EQ(100, Percent(17438564308265206, 17438564308265206)) << "'Percent' function from utilities does not return proper 100\n";
			EXPECT_EQ(11, Percent(735, 7000)) << "'Percent' function from utilities does not round up properly\n";
			EXPECT_EQ(10, Percent(734, 7000)) << "'Percent' function from utilities does not round down properly\n";
		}

		TEST(utilitiesTest, MeanWithRoundingToFirst){
			EXPECT_EQ(2, MeanWithRoundingToFirst(2, 2)) << "'MeanWithRoundingToFirst' function from utilities does not return proper value without rounding\n";
			EXPECT_EQ(3, MeanWithRoundingToFirst(3, 2)) << "'MeanWithRoundingToFirst' function from utilities does not round up properly\n";
			EXPECT_EQ(2, MeanWithRoundingToFirst(2, 3)) << "'MeanWithRoundingToFirst' function from utilities does not round down properly\n";
		}

		TEST(utilitiesTest, SetToMinMax){
			uint64_t test = 5;

			SetToMin(test, 6);
			EXPECT_EQ(5, test);

			SetToMin(test, 5);
			EXPECT_EQ(5, test);

			SetToMin(test, 4);
			EXPECT_EQ(4, test);

			test = 5;

			SetToMax(test, 4);
			EXPECT_EQ(5, test);

			SetToMax(test, 5);
			EXPECT_EQ(5, test);

			SetToMax(test, 6);
			EXPECT_EQ(6, test);
		}
	}
}
