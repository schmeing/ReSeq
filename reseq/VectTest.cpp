#include "VectTest.h"
using reseq::VectTest;

#include <stdint.h>
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "gtest/gtest.h"

void VectTest::Register(){
	// Guarantees that library is included
}

namespace reseq{
	TEST(VectTest, BasicFunctionality){
		Vect<uint32_t> test;
		EXPECT_EQ(0, test.vec_.first) << "Offset not properly initialized in class 'Vect'.\n";
		EXPECT_TRUE( test.vec_.second.empty() ) << "Std::vector not properly initialized in class 'Vect'.\n";
		EXPECT_TRUE(0 == test) << "Equality operator for class 'Vect' had wrong result for empty vec_.\n";
		EXPECT_NO_THROW( test.Shrink() ) << "Crashed during shrinking of empty object\n";

		test[5] = 15;
		EXPECT_EQ(6, test.vec_.second.size()) << "Error assuring correct size at assessing element in class 'Vect'.\n";
		EXPECT_EQ(15, test.vec_.second[5]) << "Error setting element in class 'Vect'.\n";

		test.SetOffset(6);
		EXPECT_EQ(1, test.vec_.second.size()) << "Error assuring correct size at setting offset in class 'Vect'.\n";
		EXPECT_EQ(5, test.vec_.first) << "Error setting highest possible offset in class 'Vect'.\n";
		EXPECT_TRUE(1 == test) << "Equality operator for class 'Vect' had wrong result for non-empty vec_.\n";

		test[7] = 25;
		EXPECT_EQ(3, test.vec_.second.size()) << "Error assuring correct size at assessing element after setting offset in class 'Vect'.\n";
		EXPECT_EQ(25, test.vec_.second[2]) << "Error setting element after setting offset in class 'Vect'.\n";
		EXPECT_EQ(0, test[3]) << "Error initializing entries lower than offset in class 'Vect'.\n";
		EXPECT_EQ(0, test[9]) << "Error initializing entries higher than original size in class 'Vect'.\n";

		test.Shrink();
		EXPECT_EQ(3, test.size()) << "Size not correct after shrinking for class 'Vect'.\n";
		EXPECT_EQ(5, test.vec_.first) << "Offset not correct after shrinking for class 'Vect'.\n";
		EXPECT_EQ(25, test[7]) << "Data not intact after shrinking for class 'Vect'.\n";

		Vect<uint32_t> test2;
		test2[5] = 5;
		test2[6] = 1;
		test2 += test;
		EXPECT_EQ(20, test2[5]) << "operator+= not working properly for class 'Vect'.\n";
		EXPECT_EQ(1, test2[6]) << "operator+= not working properly for class 'Vect'.\n";
		EXPECT_EQ(25, test2[7]) << "operator+= not working properly for class 'Vect'.\n";

		test2 /= 2;
		EXPECT_EQ(10, test2[5]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(0, test2[6]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(12, test2[7]) << "operator/= not working properly for class 'Vect'.\n";

		EXPECT_EQ(12, test2.Max()) << "Max function not working properly for class 'Vect'.\n";

		test2 /= test2[7];
		EXPECT_EQ(0, test2[5]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(0, test2[6]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(1, test2[7]) << "operator/= not working properly for class 'Vect'.\n";

		test2 = test;
		test2 /= test2[5];
		EXPECT_EQ(1, test2[5]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(0, test2[6]) << "operator/= not working properly for class 'Vect'.\n";
		EXPECT_EQ(1, test2[7]) << "operator/= not working properly for class 'Vect'.\n";

		test.Set(3,7,10);
		EXPECT_EQ(10, test[3]) << "Set function did not work properly.\n";
		EXPECT_EQ(10, test[6]) << "Set function did not work properly.\n";
		EXPECT_EQ(25, test[7]) << "Set function did not work properly.\n";
		EXPECT_EQ(3, test.vec_.first) << "Offset not correct after set function.\n";

		test.Shift(2);
		EXPECT_EQ(10, test[5]) << "Positive shift did not work properly.\n";
		EXPECT_EQ(5, test.vec_.first) << "Offset not correct after positive shift.\n";
		test.Shift(-1);
		EXPECT_EQ(10, test[4]) << "Negative shift did not work properly.\n";
		EXPECT_EQ(4, test.vec_.first) << "Offset not correct after negative shift.\n";
		test.Shift(-100);
		EXPECT_EQ(10, test[0]) << "Blocking to big negative shifts did not work properly.\n";
		EXPECT_EQ(25, test[4]) << "Blocking to big negative shifts did not work properly.\n";
	}

	TEST(VectTest, CopyAndClear){
		// Prepare test Vect
		Vect<uint32_t> test;
		test.vec_.first = 3;
		test.vec_.second.push_back(5);
		test.vec_.second.push_back(0);
		test.vec_.second.push_back(3);
		// Copy it
		Vect<uint32_t> test2(test);
		Vect<uint32_t> test3;
		test3 = test;
		EXPECT_EQ(test.size(), test3.size()) << "Size of copy is not identical for class 'Vect'.\n";
		for( decltype( test.vec_.first ) i = test.from(); i < test.to(); ++i ){
			EXPECT_EQ(test[i], test3[i]) << "Copy is not identical for class 'Vect'.\n";
		}
		EXPECT_EQ(test.from(), test3.from()) << "Offset of copy is not identical for class 'Vect'.\n";
		// Clear the original
		test.Clear();
		EXPECT_TRUE( test.empty() ) << "Std:vector not properly cleared for class 'Vect'.\n";
		EXPECT_EQ(0, test.vec_.first) << "Offset not properly reset for class 'Vect'.\n";
		EXPECT_FALSE( test2.empty() ) << "Reference instead of deep copy was created in copy constructor for class 'Vect'.\n";
		EXPECT_FALSE( test3.empty() ) << "Reference instead of deep copy was created in operator= for class 'Vect'.\n";
		EXPECT_EQ(test2.size(), test3.size()) << "Size of equal copies is not identical for class 'Vect'.\n";
		for( decltype( test.vec_.first ) i = test2.from(); i < test2.to(); ++i ){
			EXPECT_EQ(test2[i], test3[i]) << "Equal copies are not identical for class 'Vect'.\n";
		}
		EXPECT_EQ(test2.from(), test3.from()) << "Offset of equal copies is not identical for class 'Vect'.\n";

		// Add 0 at the beginning in test2 to change offset but no values
		test2[0] = 0;
		EXPECT_TRUE( test2 == test3 ) << "Comparison operator not working properly for class 'Vect'.\n";
		// Add 1 at the beginning in test2 to change values
		test2[0] = 1;
		EXPECT_FALSE( test2 == test3 ) << "Comparison operator not working properly for class 'Vect'.\n";
	}

	TEST(VectTest, Offset){
		Vect<uint16_t> test; // command
		test.SetOffset(10); // command
		test[15] = 3; // command
		EXPECT_EQ(6, test.size()) << "Setting offset before assigning the first value does not result in correct size for class 'Vect'.\n";
		EXPECT_EQ(10, test.from()) << "Setting offset before assigning the first value does not result in correct offset for class 'Vect'.\n";

		test[15] = 0; // command
		test.Shrink(); // command
		EXPECT_EQ(0, test.size()) << "Shrinking all zero vector does not result in correct size for class 'Vect'.\n";
		EXPECT_EQ(10, test.from()) << "Shrinking all zero vector does not result in correct offset for class 'Vect'.\n";
	}

	TEST(VectTest, MultiDimFunctions){
		// ShrinkVect
		Vect< Vect< Vect< Vect<uint16_t> > > > test_4d;
		test_4d[5][3][7][1] = 3;
		test_4d[5][3][2][1] = 0;
		test_4d[7][2][1][0] = 0;
		string setup_error = "Test setup failed\n";
		EXPECT_EQ(8, test_4d.size()) << setup_error;
		EXPECT_EQ(4, test_4d[5].size()) << setup_error;
		EXPECT_EQ(8, test_4d[5][3].size()) << setup_error;
		EXPECT_EQ(2, test_4d[5][3][7].size()) << setup_error;
		EXPECT_EQ(3, test_4d[7].size()) << setup_error;
		EXPECT_EQ(2, test_4d[7][2].size()) << setup_error;
		EXPECT_EQ(1, test_4d[7][2][1].size()) << setup_error;

		ShrinkVect(test_4d);
		string test_error = "Multidimensional shrinking did not work properly\n";
		EXPECT_EQ(1, test_4d.size()) << test_error;
		EXPECT_EQ(1, test_4d[5].size()) << test_error;
		EXPECT_EQ(1, test_4d[5][3].size()) << test_error;
		EXPECT_EQ(1, test_4d[5][3][7].size()) << test_error;

		// SumVect
		test_4d[5][3][2][2] = 7;
		EXPECT_EQ(10, SumVect(test_4d)) << "Multidimensional summing did not work properly\n";

		vector<vector<vector<double>>> stest_3d;
		stest_3d.resize(3);
		for(auto &sub2d: stest_3d){
			sub2d.resize(3);
			for(auto &sub: sub2d){
				sub.resize(3, 0.0);
			}
		}
		stest_3d[0][0][0] = 2.4;
		stest_3d[0][1][2] = 2.4;
		stest_3d[2][2][2] = 2.4;
		EXPECT_DOUBLE_EQ(7.2, SumVectD(stest_3d)) << "Multidimensional summing did not work properly for std and double\n";

		// GetLimitsSecondDimension
		Vect< Vect<uint64_t> > test_2d;
		test_2d[3][2] = 1;
		test_2d[3][6] = 1;
		test_2d[5][3] = 1;
		test_2d[5][7] = 1;
		ShrinkVect(test_2d);

		decltype(test_2d.from()) from, to;
		GetLimitsSecondDimension(test_2d, from, to);
		EXPECT_EQ(2, from) << "GetLimitsSecondDimension did not work properly\n";
		EXPECT_EQ(8, to) << "GetLimitsSecondDimension did not work properly\n";

		test_2d[7][3] = 1;
		test_2d[7].Shrink();
		GetLimitsSecondDimension(test_2d, from, to);
		EXPECT_EQ(2, from) << "GetLimitsSecondDimension did not properly ignore empty bins\n";
		EXPECT_EQ(8, to) << "GetLimitsSecondDimension did not properly ignore empty bins\n";

		// GetIndicesFirstDimension
		vector<uint64_t> indices;
		GetIndicesFirstDimension(test_2d, indices);
		EXPECT_EQ(3,indices.size()) << "GetIndicesFirstDimension did not work properly\n";
		EXPECT_EQ(3,indices[0]) << "GetIndicesFirstDimension did not work properly\n";
		EXPECT_EQ(5,indices[1]) << "GetIndicesFirstDimension did not work properly\n";
		EXPECT_EQ(7,indices[2]) << "GetIndicesFirstDimension did not work properly\n";

		// GetIndicesSecondDimension
		GetIndicesSecondDimension(test_2d, indices);
		EXPECT_EQ(4,indices.size()) << "GetIndicesSecondDimension did not work properly\n";
		EXPECT_EQ(2,indices[0]) << "GetIndicesSecondDimension did not work properly\n";
		EXPECT_EQ(3,indices[1]) << "GetIndicesSecondDimension did not work properly\n";
		EXPECT_EQ(6,indices[2]) << "GetIndicesSecondDimension did not work properly\n";
		EXPECT_EQ(7,indices[3]) << "GetIndicesSecondDimension did not work properly\n";

		// ClearSecondDimension
		auto capacity = test_2d[3].capacity();

		ClearSecondDimension(test_2d);
		EXPECT_EQ(5,test_2d.size()) << "ClearSecondDimension did not work properly\n";
		EXPECT_EQ(3,test_2d.from()) << "ClearSecondDimension did not work properly\n";
		EXPECT_EQ(0,test_2d[3].size()) << "ClearSecondDimension did not work properly\n";
		EXPECT_EQ(0,test_2d[5].size()) << "ClearSecondDimension did not work properly\n";
		EXPECT_EQ(capacity,test_2d[3].capacity()) << "ClearSecondDimension changed capacity of Vect\n";

		// Shift2d
		test_2d.Clear();
		test_2d[5][3] = 3;
		ShrinkVect(test_2d);
		Shift2d(test_2d, 3, -1);
		EXPECT_EQ(3, test_2d[8][2]) << "Multidimensional shifting did not work properly\n";

		// SetVect2d
		SetVect2d(test_2d, 5, 9, 1, 4, (uint64_t)10);
		test_error = "Twodimensional setting did not work properly\n";
		EXPECT_EQ(10, test_2d[8][2]) << test_error;
		EXPECT_EQ(10, test_2d[5][3]) << test_error;
		EXPECT_EQ(0, test_2d[4][3]) << test_error;
		EXPECT_EQ(0, test_2d[5][4]) << test_error;
		EXPECT_EQ(120, SumVect(test_2d)) << test_error;


	}
}
