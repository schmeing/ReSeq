#ifndef BASICTESTCLASS_HPP
#define BASICTESTCLASS_HPP

#include <algorithm>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "reportingUtils.hpp"

#include "CMakeConfig.h"
#include "Reference.h"
#include "SeqQualityStats.hpp"

namespace reseq{
	class BasicTestClass : public ::testing::Test{
	private:
		const uint16_t reduced_verbosity_ = 2;
		uint16_t real_verbosity_;

	protected:
		inline void ReduceVerbosity(){ kVerbosityLevel = reduced_verbosity_; };
		inline void RestoreVerbosity(){ kVerbosityLevel = real_verbosity_; };

		virtual void SetUp(){
			ReduceVerbosity();
		}
		virtual void TearDown(){
			RestoreVerbosity();
		}

		template<typename T> static void TestVectorEquality(
				const std::vector<T> &vector1,
				const std::vector<T> &vector2,
				const char *context,
				const std::string &err_msg_part1,
				const std::string &err_msg_part2
				){
			EXPECT_EQ(vector1.size(), vector2.size()) << err_msg_part1 << ".size()" << err_msg_part2 << context << '\n';
			for( uint64_t i=0; i < std::min(vector1.size(), vector2.size()); ++i ){
				EXPECT_EQ(vector1.at(i), vector2.at(i)) << err_msg_part1 << '[' << i << ']' << err_msg_part2 << context << '\n';
			}
		}

		template<typename T> static void TestVectEquality(
				const T &vector1,
				const T &vector2,
				const char *context,
				const std::string &err_msg_part1,
				const std::string &err_msg_part2
				){
			EXPECT_EQ(vector1.size(), vector2.size()) << err_msg_part1 << ".size()" << err_msg_part2 << context << '\n';
			for( auto i=vector1.from(); i < vector1.to(); ++i ){
				EXPECT_EQ(vector1[i], vector2[i]) << err_msg_part1 << '[' << i << ']' << err_msg_part2 << context << '\n';
			}
		}

		template<typename T> static void TestDoubleVectEquality(
				const T &vector1,
				const T &vector2,
				const char *context,
				const std::string &err_msg_part1,
				const std::string &err_msg_part2,
				double precision = 0.0
				){
			EXPECT_EQ(vector1.size(), vector2.size()) << err_msg_part1 << ".size()" << err_msg_part2 << context << '\n';
			for( auto i=vector1.from(); i < vector1.to(); ++i ){
				if( precision ){
					EXPECT_NEAR(vector1[i], vector2[i], precision) << err_msg_part1 << '[' << i << ']' << err_msg_part2 << context << '\n';
				}
				else{
					EXPECT_DOUBLE_EQ(vector1[i], vector2[i]) << err_msg_part1 << '[' << i << ']' << err_msg_part2 << context << '\n';
				}
			}
		}

		template<typename T> static void TestVectEquality2d(
				const T &vector1,
				const T &vector2,
				const char *context,
				const std::string &err_msg_part1,
				const std::string &err_msg_part2
				){
			EXPECT_EQ(vector1.size(), vector2.size()) << err_msg_part1 << ".size()" << err_msg_part2 << context << '\n';
			for( auto i=vector1.from(); i < vector1.to(); ++i ){
				EXPECT_EQ(vector1[i].size(), vector2[i].size()) << err_msg_part1 << '[' << i << ']' << ".size()" << err_msg_part2 << context << '\n';
				for( auto j=vector1[i].from(); j < vector1[i].to(); ++j ){
					EXPECT_EQ(vector1[i][j], vector2[i][j]) << err_msg_part1 << '[' << i << ']' << '[' << j << ']' << err_msg_part2 << context << '\n';
				}
			}
		}

		template<typename T> static void TestSeqQualityStats(
				const SeqQualityStats<T> &stats,
				uint16_t size,
				unsigned char mean,
				unsigned char minimum,
				unsigned char median,
				unsigned char maximum,
				const char *context,
				const std::string &err_msg
				){
			EXPECT_EQ(size, stats.size()) << "size " << err_msg << context << '\n';
			EXPECT_EQ(mean, stats.mean_) << "mean " << err_msg << context << '\n';
			EXPECT_EQ(minimum, stats.minimum_) << "minimum " << err_msg << context << '\n';
			EXPECT_EQ(median, stats.median_) << "median" << err_msg << context << '\n';
			EXPECT_EQ(maximum, stats.maximum_) << "maximum " << err_msg << context << '\n';
		}

	public:
		BasicTestClass():
			real_verbosity_(kVerbosityLevel)
			{
		}
	};

	class BasicTestClassWithReference : public BasicTestClass{
	protected:
		Reference species_reference_;

		void LoadReference(const std::string &ref_file){
			ASSERT_TRUE( species_reference_.ReadFasta( (std::string(PROJECT_SOURCE_DIR)+"/test/"+ref_file).c_str() ) );
		}
	};
}

#endif // REFERENCETEST_H
