#include "utilitiesTest.h"
using reseq::utilitiesTest;
using reseq::utilities::VectorAtomic;

//include <array>
using std::array;
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

		TEST(utilitiesTest, DominantBase){
			Dna5String seq = "CAGATTTTGGAANAGTNN";

			DominantBase test, test2;
			test.Set(seq,0);
			EXPECT_EQ(1, static_cast<uintBaseCall>(test.Get()));
			test.Update(at(seq,0), seq, 0);
			EXPECT_EQ(1, static_cast<uintBaseCall>(test.Get()));
			test.Clear();
			test.Set(seq,1);
			EXPECT_EQ(1, static_cast<uintBaseCall>(test.Get()));
			test2 = test;
			EXPECT_EQ(1, static_cast<uintBaseCall>(test2.Get()));
			test.Clear();
			test.Set(seq,2);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test.Get()));
			test2.Update(at(seq,1), seq, 1);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test2.Get()));
			test2 = test;
			EXPECT_EQ(0, static_cast<uintBaseCall>(test2.Get()));

			array<uintBaseCall, 15> correct_results = {2,0,0,3,3,3,3,3,2,0,0,0,0,0,3};
			for(auto pos=3; pos < length(seq); ++pos){
				test.Clear();
				test.Set(seq,pos);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test.Get())) << "Position " << pos;
				test2.Update(at(seq,pos-1), seq, pos-1);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test2.Get())) << "Position " << pos;
			}
			test2.Update(at(seq,length(seq)-1), seq, length(seq)-1); // Should run without crashing, output doesn't matter as it won't be used

			DominantBaseWithMemory test3, test4;
			test4.Update(at(seq,0));
			EXPECT_EQ(1, static_cast<uintBaseCall>(test4.Get()));
			test3.Set(seq,0);
			EXPECT_EQ(1, static_cast<uintBaseCall>(test3.Get()));
			test3.Clear();
			test3.Set(seq,1);
			EXPECT_EQ(1, static_cast<uintBaseCall>(test3.Get()));
			test4.Update(at(seq,1));
			EXPECT_EQ(1, static_cast<uintBaseCall>(test4.Get()));
			test3.Clear();
			test3.Set(seq,2);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test3.Get()));
			test4.Update(at(seq,2));
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));
			test4.Clear();
			test4 = test3;
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));

			for(auto pos=3; pos < length(seq); ++pos){
				test3.Clear();
				test3.Set(seq,pos);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test3.Get())) << "Position " << pos;
				test4.Update(at(seq,pos));
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test4.Get())) << "Position " << pos;
			}

			Dna5StringReverseComplement seq2(seq); //"NNACTNTTCCAAAATCTG"

			test.Clear();
			test.Set(seq2,0);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test.Get()));
			test.Update(at(seq2,0), seq2, 0);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test.Get()));
			test.Clear();
			test.Set(seq2,1);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test.Get()));
			test2 = test;
			EXPECT_EQ(0, static_cast<uintBaseCall>(test2.Get()));
			test.Clear();
			test.Set(seq2,2);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test.Get()));
			test2.Update(at(seq2,1), seq2, 1);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test2.Get()));
			test2 = test;
			EXPECT_EQ(0, static_cast<uintBaseCall>(test2.Get()));

			correct_results = {0,1,3,3,3,3,3,1,1,0,0,0,0,0,3};
			for(auto pos=3; pos < length(seq2); ++pos){
				test.Clear();
				test.Set(seq2,pos);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test.Get())) << "Position " << pos;
				test2.Update(at(seq2,pos-1), seq2, pos-1);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test2.Get())) << "Position " << pos;
			}
			test2.Update(at(seq2,length(seq2)-1), seq2, length(seq2)-1); // Should run without crashing, output doesn't matter as it won't be used

			test4.Clear();
			test4.Update(at(seq2,0));
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));
			test3.Clear();
			test3.Set(seq2,0);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test3.Get()));
			test3.Clear();
			test3.Set(seq2,1);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test3.Get()));
			test4.Update(at(seq2,1));
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));
			test3.Clear();
			test3.Set(seq2,2);
			EXPECT_EQ(0, static_cast<uintBaseCall>(test3.Get()));
			test4.Update(at(seq2,2));
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));
			test4.Clear();
			test4 = test3;
			EXPECT_EQ(0, static_cast<uintBaseCall>(test4.Get()));

			for(auto pos=3; pos < length(seq2); ++pos){
				test3.Clear();
				test3.Set(seq2,pos);
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test3.Get())) << "Position " << pos;
				test4.Update(at(seq2,pos));
				EXPECT_EQ(correct_results.at(pos-3), static_cast<uintBaseCall>(test4.Get())) << "Position " << pos;
			}

		}

		TEST(utilitiesTest, Divide){
			EXPECT_EQ(0, Divide(0, 2)) << "'Divide' function from utilities does not return proper zero\n";
			EXPECT_EQ(1, Divide(17438564308265206, 17438564308265206)) << "'Divide' function from utilities does not return proper 1\n";
			EXPECT_EQ(11, Divide(73500, 7000)) << "'Divide' function from utilities does not round up properly\n";
			EXPECT_EQ(10, Divide(73400, 7000)) << "'Divide' function from utilities does not round down properly\n";

			EXPECT_EQ(0, DivideAndCeil(0, 2));
			EXPECT_EQ(1, DivideAndCeil(17438564308265206, 17438564308265206));
			EXPECT_EQ(11, DivideAndCeil(70001, 7000));
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
