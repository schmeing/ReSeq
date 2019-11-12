#include "utilitiesTest.h"
using reseq::utilitiesTest;
using reseq::utilities::VectorAtomic;

#include <stdint.h>

#include "gtest/gtest.h"

//include <seqan/modifier.h>
using seqan::Dna5StringReverseComplement;
//include <seqan/sequence.h>
using seqan::Dna5String;

using seqan::reverseComplement;

void utilitiesTest::Register(){
	// Guarantees that library is included
}

namespace reseq{
	namespace utilities{
		TEST(utilitiesTest, VectorAtomic){
			VectorAtomic<unsigned int> test;
			EXPECT_EQ(100, test += 100);
			EXPECT_EQ(80, test += -20);
		}

		TEST(utilitiesTest, Divide){
			EXPECT_EQ(0, Divide(0, 2)) << "'Divide' function from utilities does not return proper zero\n";
			EXPECT_EQ(1, Divide(17438564308265206, 17438564308265206)) << "'Divide' function from utilities does not return proper 1\n";
			EXPECT_EQ(11, Divide(73500, 7000)) << "'Divide' function from utilities does not round up properly\n";
			EXPECT_EQ(10, Divide(73400, 7000)) << "'Divide' function from utilities does not round down properly\n";
		}

		TEST(utilitiesTest, HasN){
			Dna5String test = "ACGTgcat";
			Dna5StringReverseComplement test2(test);

			EXPECT_FALSE( HasN(test) );
			EXPECT_FALSE( HasN(test2) );

			test = "ACGTNgcat";
			EXPECT_TRUE( HasN(test) );
			EXPECT_TRUE( HasN(test2) );

			test = "ACGTgcatn";
			EXPECT_TRUE( HasN(test) );
			EXPECT_TRUE( HasN(test2) );

			test = "ACGTgcatS";
			EXPECT_TRUE( HasN(test) );
			EXPECT_TRUE( HasN(test2) );

			test = "ACGTgcatw";
			EXPECT_TRUE( HasN(test) );
			EXPECT_TRUE( HasN(test2) );
		}

		TEST(utilitiesTest, Percent){
			EXPECT_EQ(0, Percent(0, 2)) << "'Percent' function from utilities does not return proper zero\n";
			EXPECT_EQ(100, Percent(17438564308265206, 17438564308265206)) << "'Percent' function from utilities does not return proper 100\n";
			EXPECT_EQ(11, Percent(735, 7000)) << "'Percent' function from utilities does not round up properly\n";
			EXPECT_EQ(10, Percent(734, 7000)) << "'Percent' function from utilities does not round down properly\n";

			EXPECT_EQ(50, SafePercent(734, 0)) << "'SafePercent' function from utilities does not catch the zero denominator case properly\n";
		}

		TEST(utilitiesTest, MeanWithRoundingToFirst){
			EXPECT_EQ(2, MeanWithRoundingToFirst(2, 2)) << "'MeanWithRoundingToFirst' function from utilities does not return proper value without rounding\n";
			EXPECT_EQ(3, MeanWithRoundingToFirst(3, 2)) << "'MeanWithRoundingToFirst' function from utilities does not round up properly\n";
			EXPECT_EQ(2, MeanWithRoundingToFirst(2, 3)) << "'MeanWithRoundingToFirst' function from utilities does not round down properly\n";
		}

		TEST(utilitiesTest, IsGCN){
			Dna5String test = "ACGTNgcatnSw";

			EXPECT_FALSE( IsGC(test, 0) );
			EXPECT_TRUE( IsGC(test, 1) );
			EXPECT_TRUE( IsGC(test, 2) );
			EXPECT_FALSE( IsGC(test, 3) );
			EXPECT_FALSE( IsGC(test, 4) );
			EXPECT_TRUE( IsGC(test, 5) );
			EXPECT_TRUE( IsGC(test, 6) );
			EXPECT_FALSE( IsGC(test, 7) );
			EXPECT_FALSE( IsGC(test, 8) );
			EXPECT_FALSE( IsGC(test, 9) );
			EXPECT_FALSE( IsGC(test, 10) );
			EXPECT_FALSE( IsGC(test, 11) );

			EXPECT_FALSE( IsN(test, 0) );
			EXPECT_FALSE( IsN(test, 1) );
			EXPECT_FALSE( IsN(test, 2) );
			EXPECT_FALSE( IsN(test, 3) );
			EXPECT_TRUE( IsN(test, 4) );
			EXPECT_FALSE( IsN(test, 5) );
			EXPECT_FALSE( IsN(test, 6) );
			EXPECT_FALSE( IsN(test, 7) );
			EXPECT_FALSE( IsN(test, 8) );
			EXPECT_TRUE( IsN(test, 9) );
			EXPECT_TRUE( IsN(test, 10) );
			EXPECT_TRUE( IsN(test, 11) );

			Dna5StringReverseComplement test2(test);

			EXPECT_FALSE( IsGC(test2, 11) );
			EXPECT_TRUE( IsGC(test2, 10) );
			EXPECT_TRUE( IsGC(test2, 9) );
			EXPECT_FALSE( IsGC(test2, 8) );
			EXPECT_FALSE( IsGC(test2, 7) );
			EXPECT_TRUE( IsGC(test2, 6) );
			EXPECT_TRUE( IsGC(test2, 5) );
			EXPECT_FALSE( IsGC(test2, 4) );
			EXPECT_FALSE( IsGC(test2, 3) );
			EXPECT_FALSE( IsGC(test2, 2) );
			EXPECT_FALSE( IsGC(test2, 1) );
			EXPECT_FALSE( IsGC(test2, 0) );

			EXPECT_FALSE( IsN(test2, 11) );
			EXPECT_FALSE( IsN(test2, 10) );
			EXPECT_FALSE( IsN(test2, 9) );
			EXPECT_FALSE( IsN(test2, 8) );
			EXPECT_TRUE( IsN(test2, 7) );
			EXPECT_FALSE( IsN(test2, 6) );
			EXPECT_FALSE( IsN(test2, 5) );
			EXPECT_FALSE( IsN(test2, 4) );
			EXPECT_FALSE( IsN(test2, 3) );
			EXPECT_TRUE( IsN(test2, 2) );
			EXPECT_TRUE( IsN(test2, 1) );
			EXPECT_TRUE( IsN(test2, 0) );

			reverseComplement(test);
			EXPECT_FALSE( IsGC(test, 11) );
			EXPECT_TRUE( IsGC(test, 10) );
			EXPECT_TRUE( IsGC(test, 9) );
			EXPECT_FALSE( IsGC(test, 8) );
			EXPECT_FALSE( IsGC(test, 7) );
			EXPECT_TRUE( IsGC(test, 6) );
			EXPECT_TRUE( IsGC(test, 5) );
			EXPECT_FALSE( IsGC(test, 4) );
			EXPECT_FALSE( IsGC(test, 3) );
			EXPECT_FALSE( IsGC(test, 2) );
			EXPECT_FALSE( IsGC(test, 1) );
			EXPECT_FALSE( IsGC(test, 0) );

			EXPECT_FALSE( IsN(test, 11) );
			EXPECT_FALSE( IsN(test, 10) );
			EXPECT_FALSE( IsN(test, 9) );
			EXPECT_FALSE( IsN(test, 8) );
			EXPECT_TRUE( IsN(test, 7) );
			EXPECT_FALSE( IsN(test, 6) );
			EXPECT_FALSE( IsN(test, 5) );
			EXPECT_FALSE( IsN(test, 4) );
			EXPECT_FALSE( IsN(test, 3) );
			EXPECT_TRUE( IsN(test, 2) );
			EXPECT_TRUE( IsN(test, 1) );
			EXPECT_TRUE( IsN(test, 0) );

			EXPECT_FALSE( IsGC(test2, 0) );
			EXPECT_TRUE( IsGC(test2, 1) );
			EXPECT_TRUE( IsGC(test2, 2) );
			EXPECT_FALSE( IsGC(test2, 3) );
			EXPECT_FALSE( IsGC(test2, 4) );
			EXPECT_TRUE( IsGC(test2, 5) );
			EXPECT_TRUE( IsGC(test2, 6) );
			EXPECT_FALSE( IsGC(test2, 7) );
			EXPECT_FALSE( IsGC(test2, 8) );
			EXPECT_FALSE( IsGC(test2, 9) );
			EXPECT_FALSE( IsGC(test2, 10) );
			EXPECT_FALSE( IsGC(test2, 11) );

			EXPECT_FALSE( IsN(test2, 0) );
			EXPECT_FALSE( IsN(test2, 1) );
			EXPECT_FALSE( IsN(test2, 2) );
			EXPECT_FALSE( IsN(test2, 3) );
			EXPECT_TRUE( IsN(test2, 4) );
			EXPECT_FALSE( IsN(test2, 5) );
			EXPECT_FALSE( IsN(test2, 6) );
			EXPECT_FALSE( IsN(test2, 7) );
			EXPECT_FALSE( IsN(test2, 8) );
			EXPECT_TRUE( IsN(test2, 9) );
			EXPECT_TRUE( IsN(test2, 10) );
			EXPECT_TRUE( IsN(test2, 11) );
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

		TEST(utilitiesTest, Sign){
			EXPECT_EQ(1, Sign(321421));
			EXPECT_EQ(0, Sign(0));
			EXPECT_EQ(-1, Sign(-14363));
		}
	}
}
