#include "ProbabilityEstimatesTest.h"
using reseq::ProbabilityEstimatesTest;
using reseq::ProbabilityEstimatesSubClasses::DataStorage;
using reseq::ProbabilityEstimatesSubClasses::LogArrayCalc;

#include <algorithm>
using std::max;
using std::min;
//include <array>
using std::array;
#include <limits>
using std::numeric_limits;
#include <mutex>
using std::mutex;
#include <random>
using std::mt19937_64;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
//include <string>
using std::string;
//include <utility>
using std::pair;
//include <vector>
using std::vector;

//include "utilities.hpp"
using reseq::utilities::Divide;

void ProbabilityEstimatesTest::Register(){
	// Guarantees that library is included
}

inline double ProbabilityEstimatesTest::CalculateAbsolutPrecision(uintMatrixCount correct_margin, double marginal_sum, double precision){
	// Precision is defined to be always >1 so this has to be reverted first
	double sum_factor = correct_margin > marginal_sum ? precision : 1.0/precision;
	// We are interested in the absolut difference and not absolut value
	sum_factor = fabs( sum_factor - 1);
	// Precision is not evaluated for correct_margings equal to 0, because it would be infinity, as we still check for this cases everything considered a perfect fit (<epsilon) is accepted here
	return max(kEpsilon, marginal_sum*sum_factor);
}

void ProbabilityEstimatesTest::SetUpDataQual(DataStats &stats, uintBaseCall base, const vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position){
	stats.qualities_.base_quality_for_error_rate_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(0);
	stats.qualities_.base_quality_for_preceding_quality_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(1);
	stats.qualities_.preceding_quality_for_error_rate_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(2);
	stats.qualities_.base_quality_stats_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margin_quality_position;
	stats.qualities_.error_rate_for_position_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(4);
	stats.qualities_.preceding_quality_for_position_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(5);
	stats.qualities_.base_quality_for_sequence_quality_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(6);
	stats.qualities_.sequence_quality_for_error_rate_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(7);
	stats.qualities_.preceding_quality_for_sequence_quality_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(8);
	stats.qualities_.sequence_quality_for_position_per_tile_reference_.at(kTemplateSegment).at(base)[0] = margins.at(9);
}

void ProbabilityEstimatesTest::SetUpDataBaseCall(DataStats &stats, const vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position){
	stats.errors_.called_bases_by_base_quality_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(0);
	stats.errors_.called_bases_by_position_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(1);
	stats.qualities_.base_quality_stats_per_tile_per_error_reference_.at(kTemplateSegment).at(0).at(0)[0] = margin_quality_position;
	stats.errors_.called_bases_by_error_num_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(3);
	stats.errors_.error_num_by_quality_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(4);
	stats.errors_.error_num_by_position_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(5);
	stats.errors_.called_bases_by_error_rate_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(6);
	stats.qualities_.base_quality_for_error_rate_per_tile_per_error_reference_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(7);
	stats.qualities_.error_rate_for_position_per_tile_per_error_reference_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(8);
	stats.errors_.error_num_by_error_rate_per_tile_.at(kTemplateSegment).at(0).at(0)[0] = margins.at(9);
}

void ProbabilityEstimatesTest::SetUpDataDomError(DataStats &stats, const vector<Vect<Vect<uintMatrixCount>>> &margins){
	stats.coverage_.dominant_errors_by_distance_.at(0).at(0).at(0) = margins.at(0);
	stats.coverage_.dominant_errors_by_gc_.at(0).at(0).at(0) = margins.at(1);
	stats.coverage_.gc_by_distance_de_.at(0).at(0).at(0) = margins.at(2);
	stats.coverage_.dominant_errors_by_start_rates_.at(0).at(0).at(0) = margins.at(3);
	stats.coverage_.start_rates_by_distance_de_.at(0).at(0).at(0) = margins.at(4);
	stats.coverage_.start_rates_by_gc_de_.at(0).at(0).at(0) = margins.at(5);
}

void ProbabilityEstimatesTest::SetUpDataErrorRate(DataStats &stats, const vector<Vect<Vect<uintMatrixCount>>> &margins){
	stats.coverage_.error_rates_by_distance_.at(0).at(0) = margins.at(0);
	stats.coverage_.error_rates_by_gc_.at(0).at(0) = margins.at(1);
	stats.coverage_.gc_by_distance_er_.at(0).at(0) = margins.at(2);
	stats.coverage_.error_rates_by_start_rates_.at(0).at(0) = margins.at(3);
	stats.coverage_.start_rates_by_distance_er_.at(0).at(0) = margins.at(4);
	stats.coverage_.start_rates_by_gc_er_.at(0).at(0) = margins.at(5);
}

void ProbabilityEstimatesTest::SetUp(){
	BasicTestClass::SetUp();

	ResetProbabilityEstimates();
}

void ProbabilityEstimatesTest::TestLogArrayCalcNormalize(){
	LogArrayCalc<3> test, test_normalized;

	// Set dimensions
	for( uintMarginId n=3; n--;){
		test.dim_size_.at(n) = 3;
		test_normalized.dim_size_.at(n) = 3;
	}

	// Randomly fill test
	mt19937_64 rgen;
	rgen.seed(20180315);
	uniform_real_distribution<double> dist(0.5, 1.5);

	for( uintMarginId n=3;n--;){
		test.dim2_.at(n).resize(9);
		for( uintMatrixIndex p1=3; p1--; ){
			for( uintMatrixIndex p2=3; p2--; ){
				if( 1 == p1 && 1 == p2 ){
					test.dim2_.at(n).at(4) = 0.0; // Introduce zeros to see if the normalization can cope with that
				}
				else{
					test.dim2_.at(n).at(p1*3+p2) = dist(rgen);
				}

				// Introduce factors that cancel each other out that are close to the numeric limit to see whether that is a problem for the normalization
				if(2==n){
					test.dim2_.at(n).at(p1*3+p2) *= numeric_limits<double>::max()/2.0;
				}
				else if(1==n){
					test.dim2_.at(n).at(p1*3+p2) *= 2.0/numeric_limits<double>::max();
				}
			}
		}
	}

	// Copy values and renormalize new LogArray
	for( uintMarginId n=3;n--;){
		test_normalized.dim2_.at(n).resize(9);
		for( uintMatrixIndex p1=3; p1--; ){
			for( uintMatrixIndex p2=3; p2--; ){
				test_normalized.dim2_.at(n).at(p1*3+p2) = test.dim2_.at(n).at(p1*3+p2);
			}
		}
	}
	test_normalized.Normalize();

	// Now check if values got as small as they are supposed to be (Nothing can be higher than the max of the real distribution 1.5)
	for( uintMarginId n=3;n--;){
		for( uintMatrixIndex p1=3; p1--; ){
			for( uintMatrixIndex p2=3; p2--; ){
				EXPECT_TRUE( 1.6 > test_normalized.dim2_.at(n).at(p1*3+p2) ); // Chose 1.6, that is good enough and very generous towards precision issues
			}
		}
	}

	// Now check if the results are still the same
	for( uintMatrixIndex p1=3; p1--; ){
		for( uintMatrixIndex p2=3; p2--; ){
			for( uintMatrixIndex p3=3; p3--; ){
				// Double precision cannot be expected after the logs used
				EXPECT_FLOAT_EQ( test.dim2_.at(2).at(p3*3+p2)*test.dim2_.at(1).at(p3*3+p1)*test.dim2_.at(0).at(p2*3+p1),
						test_normalized.dim2_.at(2).at(p3*3+p2)*test_normalized.dim2_.at(1).at(p3*3+p1)*test_normalized.dim2_.at(0).at(p2*3+p1) );
			}
		}
	}
}

void ProbabilityEstimatesTest::TestDataStorageReduceAndExpand(){
	// Prepare test data
	DataStorage<4> test;
	for(uintMarginId n=3; n--; ){
		test.dim_size_.at(n) = 3+2*n;
	}
	test.dim_size_.at(3) = 1;

	test.data_.at(0).resize(5*3, 0.0); // 1,0
	test.data_.at(1).resize(7*3, 0.0); // 2,0
	test.data_.at(2).resize(3, 0.0); // 3,0
	test.data_.at(3).resize(7*5, 0.0); // 2,1
	test.data_.at(4).resize(5, 0.0); // 3,1
	test.data_.at(5).resize(7, 0.0); // 3,2

	double value;
	uintMatrixIndex red_i, red_j, red_k;
	for(auto i=test.dim_size_.at(0); i--; ){
		for(auto j=test.dim_size_.at(1); j--; ){
			for(auto k=test.dim_size_.at(2); k--; ){
				red_i = i;
				if( 0==j ){
					red_j = 0;
				}
				else if(4 == j){
					red_j = 2;
				}
				else{
					red_j = 1;
				}
				if( 0==k ){
					red_k = 0;
				}
				else if(3 < k){
					red_k = 2;
				}
				else{
					red_k = 1;
				}

				if( !((0==red_i && 0==red_j) || (1==red_i && 1==red_k) || (2==red_j && 2==red_k)) ){
					value = (red_i+1)*(red_j+4)*(red_k+7);

					test.data_.at(0).at(j*test.dim_size_.at(0)+i) += value;
					test.data_.at(1).at(k*test.dim_size_.at(0)+i) += value;
					test.data_.at(2).at(i) += value;
					test.data_.at(3).at(k*test.dim_size_.at(1)+j) += value;
					test.data_.at(4).at(j) += value;
					test.data_.at(5).at(k) += value;
				}
			}
		}
	}

	// Test Reduce function
	DataStorage<4> reduced_data;
	std::array<std::vector<uintMatrixIndex>, 4> dim_indices_reduced, dim_indices_count;

	test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 3);

	// Check for correctness
	for(uintMarginId n=3; n--; ){
		EXPECT_EQ(3+2*n, dim_indices_reduced.at(n).size()) << "dim_indices_reduced[" << n << "].size() wrong";
		for(auto i=dim_indices_reduced.at(n).size(); i--; ){
			// Index range is from 0 to 2
			EXPECT_TRUE(dim_indices_reduced.at(n).at(i) < 3) << "dim_indices_reduced[" << n << "][" << i << "] invalid index";
		}
		for(auto i=dim_indices_reduced.at(n).size(); --i; ){
			// First index has to be different then all other indeces
			EXPECT_TRUE(dim_indices_reduced.at(n).at(i) != dim_indices_reduced.at(n).at(0)) << "dim_indices_reduced[" << n << "][" << i << "] is equal to first index";
		}
		for(auto i=dim_indices_reduced.at(n).size()-1; i--; ){
			if(i < 4){
				// Last index has to be different then all other indeces
				EXPECT_TRUE(dim_indices_reduced.at(n).at(i) != dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-1)) << "dim_indices_reduced[" << n << "][" << i << "] is equal to last index";
			}
			else{
				// Special case for dimension 3, where the last index exists three times
				EXPECT_TRUE(dim_indices_reduced.at(n).at(i) == dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-1)) << "dim_indices_reduced[" << n << "][" << i << "] is not equal to last index";
			}
		}
		for(auto i=dim_indices_reduced.at(n).size()-(2==n?4:2); --i; ){
			// Middle indeces have to be equal
			EXPECT_TRUE(dim_indices_reduced.at(n).at(i) == dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-(2==n?4:2))) << "dim_indices_reduced[" << n << "][" << i << "] is not equal to middle indeces";
		}

		EXPECT_EQ(3, dim_indices_count.at(n).size()) << "dim_indices_count[" << n << "].size() wrong";
		EXPECT_EQ(1, dim_indices_count.at(n).at( dim_indices_reduced.at(n).at(0) )) << "dim_indices_count[" << n << "][" << dim_indices_reduced.at(n).at(0) << "] wrong";
		EXPECT_EQ((0==n?1:3), dim_indices_count.at(n).at( dim_indices_reduced.at(n).at(1) )) << "dim_indices_count[" << n << "][" << dim_indices_reduced.at(n).at(1) << "] wrong";
		EXPECT_EQ((2==n?3:1), dim_indices_count.at(n).at( dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-1) )) << "dim_indices_count[" << n << "][" << dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-1) << "] wrong";
	}

	uintMarginId dim_a(4), dim_b(3);
	for(uintMarginId n=6; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = 3;
		}

		for(auto i=dim_indices_reduced.at(dim_a).size(); i--; ){
			for(auto j=dim_indices_reduced.at(dim_b).size(); j--; ){
				EXPECT_DOUBLE_EQ( reduced_data.data_.at(n).at( dim_indices_reduced.at(dim_a).at(i)*dim_indices_count.at(dim_b).size() + dim_indices_reduced.at(dim_b).at(j) ),
						test.data_.at(n).at( i*dim_indices_reduced.at(dim_b).size() + j) * dim_indices_count.at(dim_a).at( dim_indices_reduced.at(dim_a).at(i) ) * dim_indices_count.at(dim_b).at( dim_indices_reduced.at(dim_b).at(j) ) );
			}
		}
	}

	// Test Expand function
	std::array<std::vector<uintMatrixIndex>, 4> expansion_indices;
	std::array<std::vector<uintMatrixIndex>, 4> expansion_count;
	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uintMarginId n=3; n--; ){
		EXPECT_EQ((2==n?6:3), dim_indices_count.at(n).size()) << "dim_indices_count[" << n << "].size() wrong";
		if(0==n){
			for( auto count : dim_indices_count.at(n)){
				EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after first expansion";
			}
			for( auto ind=dim_indices_reduced.at(n).size(); ind--; ){
				EXPECT_EQ(ind, dim_indices_reduced.at(n).at(ind)) << "dim_indices_reduced[" << n << "] did not go back to original";
			}
		}
		else if(2==n){
			bool count_two_left(true);
			for( auto count : dim_indices_count.at(n)){
				EXPECT_TRUE(1 == count || (count_two_left && 2 == count)) << "dim_indices_count[" << n << "] has wrong counts after first expansion";
				if(2 == count){
					count_two_left = false;
				}
			}
			std::vector<uintMatrixIndex> tmp_counts = dim_indices_count.at(n);
			for( auto ind : dim_indices_reduced.at(n)){
				EXPECT_TRUE( tmp_counts.at(ind)-- )  << "dim_indices_count[" << n << "] does not fit with dim_indices_reduced after first expansion";
			}
		}

		EXPECT_EQ( 1, expansion_count.at(n).at(expansion_indices.at(n).at(dim_indices_reduced.at(n).at(0))) ) << "expansion_count[" << n << "] wrong";
		EXPECT_EQ( (0==n?1:3), expansion_count.at(n).at(expansion_indices.at(n).at(dim_indices_reduced.at(n).at(1))) ) << "expansion_count[" << n << "] wrong";
		EXPECT_EQ( (2==n?3:1), expansion_count.at(n).at(expansion_indices.at(n).at(dim_indices_reduced.at(n).at(dim_indices_reduced.at(n).size()-1))) ) << "expansion_count[" << n << "] wrong";
	}

	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uintMarginId n=3; n--; ){
		EXPECT_EQ((2==n?6:3+2*n), dim_indices_count.at(n).size()) << "dim_indices_count[" << n << "].size() wrong";
		if(2!=n){
			for( auto count : dim_indices_count.at(n)){
				EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after second expansion";
			}
			for( auto ind=dim_indices_reduced.at(n).size(); ind--; ){
				EXPECT_EQ(ind, dim_indices_reduced.at(n).at(ind)) << "dim_indices_reduced[" << n << "] did not go back to original";
			}
		}
	}

	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uintMarginId n=3; n--; ){
		EXPECT_EQ(3+2*n, dim_indices_count.at(n).size()) << "dim_indices_count[" << n << "].size() wrong";
		for( auto count : dim_indices_count.at(n)){
			EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after third expansion";
		}
		for( auto ind=dim_indices_reduced.at(n).size(); ind--; ){
			EXPECT_EQ(ind, dim_indices_reduced.at(n).at(ind)) << "dim_indices_reduced[" << n << "] did not go back to original";
		}
	}

	EXPECT_FALSE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
}

void ProbabilityEstimatesTest::TestDataStorageReduceAndExpand2(){
	// Prepare test data
	DataStorage<4> test;
	test.dim_size_.at(0) = 1;
	test.dim_size_.at(1) = 1;
	test.dim_size_.at(2) = 5;
	test.dim_size_.at(3) = 1;
	test.data_.at(0).resize(1, 0.0); // 1,0
	test.data_.at(1).resize(5, 0.0); // 2,0
	test.data_.at(2).resize(1, 0.0); // 3,0
	test.data_.at(3).resize(5, 0.0); // 2,1
	test.data_.at(4).resize(1, 0.0); // 3,1
	test.data_.at(5).resize(5, 0.0); // 3,2

	double value;

	for(auto k=test.dim_size_.at(2); k--; ){
		if(k < 2){
			value=200+k;
		}
		else{
			value=96+2*k;
		}

		test.data_.at(0).at(0) += value;
		test.data_.at(1).at(k) += value;
		test.data_.at(2).at(0) += value;
		test.data_.at(3).at(k) += value;
		test.data_.at(4).at(0) += value;
		test.data_.at(5).at(k) += value;
	}

	// Test Reduce function
	std::array<std::vector<uintMatrixIndex>, 4> dim_indices_reduced, dim_indices_count;

	for( uint16_t test_case=0; test_case<2; ++test_case){
		for( uintMarginId n = 4; n--; ){
			dim_indices_count.at(n).clear();
			dim_indices_reduced.at(n).clear();
		}

		if(0 == test_case){
			DataStorage<4> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 2);
		}
		else{
			// Test Expand function
			DataStorage<4> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 1);

			std::array<std::vector<uintMatrixIndex>, 4> expansion_indices, expansion_count;

			EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";

			for( uintMarginId n = 3; n--; ){
				EXPECT_EQ( 1, expansion_count.at(n).size() ) << "expansion_count[" << n << "].size() wrong";
				EXPECT_EQ( (2==n?5:1), expansion_count.at(n).at(0) ) << "expansion_count[" << n << "][0] wrong";
				EXPECT_EQ( (2==n?2:1), expansion_indices.at(n).size() ) << "expansion_indices[" << n << "].size() wrong";
				for(uintMatrixIndex i=0; i < expansion_indices.at(n).size(); ++i ){
					EXPECT_EQ( 0, expansion_indices.at(n).at(i) ) << "expansion_count[" << n << "][" << i << "] wrong";
				}
			}
		}

		EXPECT_EQ(1, dim_indices_count.at(0).at(0)) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(1, dim_indices_count.at(1).at(0)) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(0, dim_indices_reduced.at(0).at(0)) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(0, dim_indices_reduced.at(1).at(0)) << "Error in simple dimension in case " << test_case;

		if(2 == dim_indices_count.at(2).at(0)){
			EXPECT_EQ(3, dim_indices_count.at(2).at(1)) << "dim_indices_count[2] wrong: " << dim_indices_count.at(2).at(0) << ' ' << dim_indices_count.at(2).at(1) << " in case " << test_case;
			for(uintMatrixIndex i=0; i<2; ++i){
				EXPECT_EQ(0, dim_indices_reduced.at(2).at(i)) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
			for(uintMatrixIndex i=2; i<5; ++i){
				EXPECT_EQ(1, dim_indices_reduced.at(2).at(i)) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
		}
		else{
			EXPECT_EQ(3, dim_indices_count.at(2).at(0)) << "dim_indices_count[2] wrong: " << dim_indices_count.at(2).at(0) << ' ' << dim_indices_count.at(2).at(1) << " in case " << test_case;
			EXPECT_EQ(2, dim_indices_count.at(2).at(1)) << "dim_indices_count[2] wrong: " << dim_indices_count.at(2).at(0) << ' ' << dim_indices_count.at(2).at(1) << " in case " << test_case;
			for(uintMatrixIndex i=0; i<2; ++i){
				EXPECT_EQ(1, dim_indices_reduced.at(2).at(i)) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
			for(uintMatrixIndex i=2; i<5; ++i){
				EXPECT_EQ(2, dim_indices_reduced.at(2).at(i)) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
		}
	}

	for( uint16_t test_case=0; test_case<2; ++test_case){
		for( uintMarginId n = 3; n--; ){
			dim_indices_count.at(n).clear();
			dim_indices_reduced.at(n).clear();
		}

		if(0 == test_case){
			DataStorage<4> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 4);
		}
		else{
			// Test Expand function
			DataStorage<4> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 1);

			std::array<std::vector<uintMatrixIndex>, 4> expansion_indices, expansion_count;

			EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 2) ) << "Expansion returns wrong value";

			for( uintMarginId n = 3; n--; ){
				EXPECT_EQ( 1, expansion_count.at(n).size() ) << "expansion_count[" << n << "].size() wrong";
				EXPECT_EQ( (2==n?5:1), expansion_count.at(n).at(0) ) << "expansion_count[" << n << "][0] wrong";
				EXPECT_EQ( (2==n?4:1), expansion_indices.at(n).size() ) << "expansion_indices[" << n << "].size() wrong";
				for(uintMatrixIndex i=0; i < expansion_indices.at(n).size(); ++i ){
					EXPECT_EQ( 0, expansion_indices.at(n).at(i) ) << "expansion_count[" << n << "][" << i << "] wrong";
				}
			}
		}

		EXPECT_EQ(dim_indices_reduced.at(2).at(0), dim_indices_reduced.at(2).at(1)) << "Max difference for the separation doesn't seem to work in case " << test_case;
		EXPECT_EQ(2, dim_indices_count.at(2).at( dim_indices_reduced.at(2).at(1) )) << "Max difference for the separation doesn't seem to work in case " << test_case;

		bool count_two_left(true);
		for( auto count : dim_indices_count.at(2)){
			EXPECT_TRUE(1 == count || (count_two_left && 2 == count)) << "dim_indices_count[2] has wrong counts in case " << test_case;
			if(2 == count){
				count_two_left = false;
			}
		}
		std::vector<uintMatrixIndex> tmp_counts = dim_indices_count.at(2);
		for( auto ind : dim_indices_reduced.at(2)){
			EXPECT_TRUE( tmp_counts.at(ind)-- )  << "dim_indices_count[2] does not fit with dim_indices_reduced in case " << test_case;
		}
	}
}

void ProbabilityEstimatesTest::TestLogArrayCalcExpand(){
	LogArrayCalc<3> test;

	// Set dimensions
	for( uintMarginId n=3; n--;){
		test.dim_size_.at(n) = 3;
	}

	// Randomly fill test
	mt19937_64 rgen;
	rgen.seed(20180320);
	uniform_real_distribution<double> dist(0.5, 1.5);

	for( uintMarginId n=3;n--;){
		test.dim2_.at(n).resize(9);
		for( uintMatrixIndex p1=3; p1--; ){
			for( uintMatrixIndex p2=3; p2--; ){
				if( n == p1 && n == p2 ){
					test.dim2_.at(n).at(4) = 0.0; // Introduce zeros to see if the extension can cope with that
				}
				else{
					test.dim2_.at(n).at(p1*3+p2) = dist(rgen);
				}
			}
		}
	}

	// Calculate marginal sums for later comparison
	std::array<std::vector<double>, 3> marginal_sums;
	for( uintMarginId n=3; n--;){
		marginal_sums.at(n).resize(9, 0.0);
	}

	double value;
	for( uintMatrixIndex p1=test.dim_size_.at(0); p1--; ){
		for( uintMatrixIndex p2=test.dim_size_.at(1); p2--; ){
			for( uintMatrixIndex p3=test.dim_size_.at(2); p3--; ){
				// Double precision cannot be expected after the logs used
				value = test.dim2_.at(2).at(p3*3+p2)*test.dim2_.at(1).at(p3*3+p1)*test.dim2_.at(0).at(p2*3+p1);

				marginal_sums.at(0).at(p2*test.dim_size_.at(0)+p1) += value;
				marginal_sums.at(1).at(p3*test.dim_size_.at(0)+p1) += value;
				marginal_sums.at(2).at(p3*test.dim_size_.at(1)+p2) += value;
			}
		}
	}

	// Expand LogArray
	std::array<std::vector<uintMatrixIndex>, 3> dim_indices_reduced, dim_indices_count;
	dim_indices_reduced.at(0) = {0,1,2};
	dim_indices_reduced.at(1) = {0,1,1,2};
	dim_indices_reduced.at(2) = {0,1,1,1,2};
	dim_indices_count.at(0) = {1,1,1};
	dim_indices_count.at(1) = {1,2,1};
	dim_indices_count.at(2) = {1,3,1};

	test.Expand(dim_indices_reduced, dim_indices_count);

	// Calculate expanded marginal sums
	std::array<std::vector<double>, 3> expanded_marginal_sums;
	expanded_marginal_sums.at(0).resize(test.dim_size_.at(1)*test.dim_size_.at(0));
	expanded_marginal_sums.at(1).resize(test.dim_size_.at(2)*test.dim_size_.at(0));
	expanded_marginal_sums.at(2).resize(test.dim_size_.at(2)*test.dim_size_.at(1));

	for( uintMatrixIndex p1=test.dim_size_.at(0); p1--; ){
		for( uintMatrixIndex p2=test.dim_size_.at(1); p2--; ){
			for( uintMatrixIndex p3=test.dim_size_.at(2); p3--; ){
				// Double precision cannot be expected after the logs used
				value = test.dim2_.at(2).at(p3*4+p2)*test.dim2_.at(1).at(p3*3+p1)*test.dim2_.at(0).at(p2*3+p1);

				expanded_marginal_sums.at(0).at(p2*test.dim_size_.at(0)+p1) += value;
				expanded_marginal_sums.at(1).at(p3*test.dim_size_.at(0)+p1) += value;
				expanded_marginal_sums.at(2).at(p3*test.dim_size_.at(1)+p2) += value;
			}
		}
	}

	// Compare marginal sums
	uintMarginId dim_a(3), dim_b(2);
	for(uintMarginId n=3; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = 2;
		}

		for(auto i=dim_indices_reduced.at(dim_a).size(); i--; ){
			for(auto j=dim_indices_reduced.at(dim_b).size(); j--; ){
				EXPECT_DOUBLE_EQ( marginal_sums.at(n).at( dim_indices_reduced.at(dim_a).at(i)*dim_indices_count.at(dim_b).size() + dim_indices_reduced.at(dim_b).at(j) ),
						expanded_marginal_sums.at(n).at( i*dim_indices_reduced.at(dim_b).size() + j) * dim_indices_count.at(dim_a).at( dim_indices_reduced.at(dim_a).at(i) ) * dim_indices_count.at(dim_b).at( dim_indices_reduced.at(dim_b).at(j) ) )
						<< '[' << n << "][" << i << "][" << j << ']';
			}
		}
	}
}

void ProbabilityEstimatesTest::TestCombineDimIndices(){
	array< vector<uintMatrixIndex>, 2 > dim_indeces1, dim_indeces2;

	for(uintMarginId cur_dim=2; cur_dim--; ){
		dim_indeces1.at(cur_dim).push_back(0);
		dim_indeces1.at(cur_dim).push_back(1);
		dim_indeces1.at(cur_dim).push_back(1);
		dim_indeces1.at(cur_dim).push_back(2);
		dim_indeces1.at(cur_dim).push_back(2);

		dim_indeces2.at(cur_dim).push_back(0);
		dim_indeces2.at(cur_dim).push_back(1);
		dim_indeces2.at(cur_dim).push_back(0);
	}

	ProbabilityEstimatesSubClasses::CombineDimIndices(dim_indeces1, dim_indeces2);

	for(uintMarginId cur_dim=2; cur_dim--; ){
		EXPECT_EQ(5, dim_indeces1.at(cur_dim).size());
		EXPECT_EQ(0, dim_indeces1.at(cur_dim).at(0));
		EXPECT_EQ(1, dim_indeces1.at(cur_dim).at(1));
		EXPECT_EQ(1, dim_indeces1.at(cur_dim).at(2));
		EXPECT_EQ(0, dim_indeces1.at(cur_dim).at(3));
		EXPECT_EQ(0, dim_indeces1.at(cur_dim).at(4));
	}
}

void ProbabilityEstimatesTest::ResetProbabilityEstimates(){
	test_.SetVariables(1);
}

void ProbabilityEstimatesTest::RandomFillingTestCounts(
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts,
		const uintMatrixIndex dim1,
		const uintMatrixIndex dim2,
		const uintMatrixIndex dim3,
		const uintMatrixIndex dim4,
		const uintMatrixIndex dim5,
		const uintSeed seed){
	mt19937_64 rgen;
	rgen.seed(seed);
	uniform_int_distribution<uintMatrixCount> dist(0, 10);
	counts.Clear();

	for( auto p1=dim1; p1--; ){
		for( auto p2=dim2; p2--; ){
			for( auto p3=dim3; p3--; ){
				for( auto p4=dim4; p4--; ){
					for( auto p5=dim5; p5--; ){
						// Chose random values for every position
						counts[p1][p2][p3][p4][p5] = dist(rgen)*100; // Multiply everything by 100 so we are above minimum counts and bins are not combined
					}
				}
			}
		}
	}

	// In case we have zeros at the borders
	ShrinkVect(counts);
}

void ProbabilityEstimatesTest::SetOffset(
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts,
		const uintMatrixIndex dim1,
		const uintMatrixIndex dim2,
		const uintMatrixIndex dim3,
		const uintMatrixIndex dim4,
		const uintMatrixIndex dim5){
	for( auto &v1 : counts ){
		for( auto &v2 : v1 ){
			for( auto &v3 : v2 ){
				for( auto &v4 : v3 ){
					v4.Shift(dim5);
				}
				v3.Shift(dim4);
			}
			v2.Shift(dim3);
		}
		v1.Shift(dim2);
	}
	counts.Shift(dim1);
}

void ProbabilityEstimatesTest::CheckIPFResult(
		uintNumFits iterations,
		double precision,
		const vector<Vect<Vect<uintMatrixCount>>> &margins,
		const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position,
		const vector< pair<bool, bool> > &margin_def,
		const Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts,
		string context ){
	EXPECT_LE( iterations, kMaxIterations ) << "Maximum iterations exceeded\n";
	EXPECT_LT( precision, kPrecisionAim ) << "Precision aim not achieved\n";

	// Get margins
	vector<Vect<Vect<double>>> estimated_margins;
	Vect<SeqQualityStats<double>> estimated_margins_quality_position;
	CalculateMargins(margin_def, estimated_margins, estimated_margins_quality_position, estimated_counts);

	// Check all defined margins
	for(uintMarginId i=0; i < margins.size(); ++i){
		if(margin_def.at(i).first){
			for( auto pos1=margin_quality_position.from(); pos1 < margin_quality_position.to(); ++pos1 ){
				for( auto pos2=margin_quality_position.at(pos1).from(); pos2 < margin_quality_position.at(pos1).to(); ++pos2 ){
					EXPECT_NEAR( margin_quality_position.at(pos1).at(pos2), estimated_margins_quality_position.at(pos1).at(pos2), CalculateAbsolutPrecision(margin_quality_position.at(pos1).at(pos2), estimated_margins_quality_position.at(pos1).at(pos2), precision) ) << "(quality_position) margins[" << i << "][" << pos1 << "][" << pos2 << "] does not match within precision for " << context << '\n';
				}
			}
		}
		else{
			for( auto pos1=margins.at(i).from(); pos1 < margins.at(i).to(); ++pos1 ){
				for( auto pos2=margins.at(i).at(pos1).from(); pos2 < margins.at(i).at(pos1).to(); ++pos2 ){
					EXPECT_NEAR( margins.at(i).at(pos1).at(pos2), estimated_margins.at(i).at(pos1).at(pos2), CalculateAbsolutPrecision(margins.at(i).at(pos1).at(pos2), estimated_margins.at(i).at(pos1).at(pos2), precision) ) << "margins[" << i << "][" << pos1 << "][" << pos2 << "] does not match within precision for " << context << '\n';
				}
			}
		}
	}
}

void ProbabilityEstimatesTest::SetMarginDefQual(vector< pair<bool, bool> > &margin_def){
	// {use margin_quality_position, flip dimensions compared to dimension order}
	margin_def.push_back({false, true}); // er, qual
	margin_def.push_back({false, true}); // pq, qual
	margin_def.push_back({false, false}); // er, pq
	margin_def.push_back({true, true}); // pos, qual
	margin_def.push_back({false, true}); // pos, er
	margin_def.push_back({false, true}); // pos, pq
	margin_def.push_back({false, true}); // sq, qual
	margin_def.push_back({false, false}); // er, sq
	margin_def.push_back({false, true}); // sq, pq
	margin_def.push_back({false, false}); // pos, sq
}

void ProbabilityEstimatesTest::SetUpDataStorageErrorRate(ProbabilityEstimatesSubClasses::DataStorage<4> &data, array< vector<uintMatrixIndex>, 4 > &dim_indices, array< vector<uintMatrixIndex>, 4 > &initial_dim_indices, const Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts){
	vector<Vect<Vect<uintMatrixCount>>> margins;
	Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
	vector< pair<bool, bool> > margin_def;
	SetMarginDefErrorRate(margin_def);
	CalculateMargins(margin_def, margins, margin_quality_position, counts);

	DataStats stats(NULL);
	SetUpDataErrorRate(stats, margins);

	array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > margin_def2;
	test_.DefineMarginsErrorRate(stats, margin_def2, 0, 0);

	mutex print_mutex;
	data.SetUp(margin_def2, NULL, dim_indices, initial_dim_indices, print_mutex);
}

void ProbabilityEstimatesTest::TestDataStorageSetUpErrorRate(ProbabilityEstimatesSubClasses::DataStorage<4> &data, array< vector<uintMatrixIndex>, 4 > &dim_indices, array< vector<uintMatrixIndex>, 4 > &initial_dim_indices, bool middle_bin_inserted){

	EXPECT_EQ( 2 , dim_indices.at(0).size() );
	EXPECT_EQ( 0 , dim_indices.at(0).at(0) );
	EXPECT_EQ( 2 , dim_indices.at(0).at(1) );

	if(middle_bin_inserted){
		EXPECT_EQ( 3 , dim_indices.at(1).size() );
		EXPECT_EQ( 1 , dim_indices.at(1).at(0) );
		EXPECT_EQ( 2 , dim_indices.at(1).at(1) );
		EXPECT_EQ( 3 , dim_indices.at(1).at(2) );
	}
	else{
		EXPECT_EQ( 2 , dim_indices.at(1).size() );
		EXPECT_EQ( 1 , dim_indices.at(1).at(0) );
		EXPECT_EQ( 2 , dim_indices.at(1).at(1) );
	}

	EXPECT_EQ( 2 , dim_indices.at(2).size() );
	EXPECT_EQ( 0 , dim_indices.at(2).at(0) );
	EXPECT_EQ( 1 , dim_indices.at(2).at(1) );

	EXPECT_EQ( 2 , data.dim_size_.at(0) );
	EXPECT_EQ( 2 , data.dim_size_.at(1) );
	EXPECT_EQ( 2 , data.dim_size_.at(2) );

	EXPECT_EQ( 4 , data.data_.at(0).size() );
	EXPECT_DOUBLE_EQ( 3.0/36.0 , data.data_.at(0).at(0) );
	EXPECT_DOUBLE_EQ( 11.0/36.0 , data.data_.at(0).at(1) );
	EXPECT_DOUBLE_EQ( 7.0/36.0 , data.data_.at(0).at(2) );
	EXPECT_DOUBLE_EQ( 15.0/36.0 , data.data_.at(0).at(3) );

	EXPECT_EQ( 4 , data.data_.at(1).size() );
	EXPECT_DOUBLE_EQ( 4.0/36.0 , data.data_.at(1).at(0) );
	EXPECT_DOUBLE_EQ( 12.0/36.0 , data.data_.at(1).at(1) );
	EXPECT_DOUBLE_EQ( 6.0/36.0 , data.data_.at(1).at(2) );
	EXPECT_DOUBLE_EQ( 14.0/36.0 , data.data_.at(1).at(3) );

	EXPECT_EQ( 4 , data.data_.at(3).size() );
	EXPECT_DOUBLE_EQ( 6.0/36.0 , data.data_.at(3).at(0) );
	EXPECT_DOUBLE_EQ( 10.0/36.0 , data.data_.at(3).at(1) );
	EXPECT_DOUBLE_EQ( 8.0/36.0 , data.data_.at(3).at(2) );
	EXPECT_DOUBLE_EQ( 12.0/36.0 , data.data_.at(3).at(3) );

	if(middle_bin_inserted){
		EXPECT_EQ( 0 , initial_dim_indices.at(1).at(0) ) << "initial_dim_indices[" << 0 << "] is wrong";
		EXPECT_EQ( 0 , initial_dim_indices.at(1).at(1) ) << "initial_dim_indices[" << 0 << "] is wrong";
		EXPECT_EQ( 1 , initial_dim_indices.at(1).at(2) ) << "initial_dim_indices[" << 0 << "] is wrong";

		for(uintMarginId cur_dim=0; cur_dim < 3; cur_dim+=2){
			EXPECT_EQ( 0 , initial_dim_indices.at(cur_dim).at(0) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
			EXPECT_EQ( 1 , initial_dim_indices.at(cur_dim).at(1) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
		}
	}
	else{
		for(uintMarginId cur_dim=0; cur_dim < 3; ++cur_dim){
			EXPECT_EQ( 0 , initial_dim_indices.at(cur_dim).at(0) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
			EXPECT_EQ( 1 , initial_dim_indices.at(cur_dim).at(1) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
		}
	}
}

void ProbabilityEstimatesTest::IterativeProportionalFittingQual(uintBaseCall base, const vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataQual(stats, base, margins, margin_quality_position);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFQuality, kTemplateSegment, 0, base, 0, 0, kMaxIterations, kPrecisionAim);
}

void ProbabilityEstimatesTest::GetIPFResultQual( const ProbabilityEstimates &estimate, uintBaseCall base, Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts, const vector<Vect<Vect<uintMatrixCount>>> &margins ){
	estimated_counts.Clear();

	auto tot_counts = SumVect(margins.at(0)); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uintMatrixIndex>, 5 > &dim_indices(estimate.quality_.at(kTemplateSegment).at(0).at(base).dim_indices_);

	// Get result matrix
	for( uintMatrixIndex q=0; q < dim_indices.at(0).size(); ++q ){
		for( uintMatrixIndex er=0; er < dim_indices.at(4).size(); ++er ){
			for( uintMatrixIndex pq=0; pq < dim_indices.at(2).size(); ++pq ){
				for( uintMatrixIndex pos=0; pos < dim_indices.at(3).size(); ++pos ){
					for( uintMatrixIndex sq=0; sq < dim_indices.at(1).size(); ++sq ){
						estimated_counts[dim_indices.at(0).at(q)][dim_indices.at(4).at(er)][dim_indices.at(2).at(pq)][dim_indices.at(3).at(pos)][dim_indices.at(1).at(sq)] = estimate.quality_.at(kTemplateSegment).at(0).at(base).estimates_.GetMatrixElement({q, sq, pq, pos, er}) * tot_counts;
					}
				}
			}
		}
	}

	ShrinkVect(estimated_counts);
}

void ProbabilityEstimatesTest::IPFStepWiseQual(uintBaseCall base, const vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position){
	DataStats stats(NULL);
	SetUpDataQual(stats, base, margins, margin_quality_position);

	std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 10 > margin_defs;
	auto alternative_margin = test_.DefineMarginsQuality(stats, margin_defs, kTemplateSegment, 0, base);

	array< vector<uintMatrixIndex>, 5 > dim_indices;
	array< vector<uintMatrixIndex>, 5 > initial_dim_indices;
	DataStorage<5> data;
	mutex print_mutex;
	data.SetUp(margin_defs, alternative_margin, dim_indices, initial_dim_indices, print_mutex);

	test_.quality_.at(kTemplateSegment).at(0).at(base).estimates_.SetUp(dim_indices);

	IPFStepQual<1,0,2,3,4>(base, data); // (seq_qual, qual) margin
	IPFStepQual<2,0,1,3,4>(base, data); // (prev_qual, qual) margin
	IPFStepQual<3,0,1,2,4>(base, data); // (pos, qual) margin
	IPFStepQual<4,0,1,2,3>(base, data); // (error_rate, qual) margin
	IPFStepQual<2,1,0,3,4>(base, data); // (prev_qual, seq_qual) margin
	IPFStepQual<3,1,0,2,4>(base, data); // (pos, seq_qual) margin
	IPFStepQual<4,1,0,2,3>(base, data); // (error_rate, seq_qual) margin
	IPFStepQual<2,3,0,1,4>(base, data); // (pos, prev_qual) margin: Flipped indeces for U1, U2 (Test here that it doesn't matter)
	IPFStepQual<4,2,0,1,3>(base, data); // (error_rate, prev_qual) margin
	IPFStepQual<4,3,0,1,2>(base, data); // (pos, error_rate) margin
}

void ProbabilityEstimatesTest::IterativeProportionalFittingBaseCall( const vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position ){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataBaseCall(stats, margins, margin_quality_position);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFBaseCall, kTemplateSegment, 0, 0, 0, 0, kMaxIterations, kPrecisionAim);
}

void ProbabilityEstimatesTest::GetIPFResultBaseCall( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts, vector<Vect<Vect<uintMatrixCount>>> &margins ){
	estimated_counts.Clear();

	auto tot_counts = SumVect(margins.at(0)); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uintMatrixIndex>, 5 > &dim_indices(estimate.base_call_.at(kTemplateSegment).at(0).at(0).at(0).dim_indices_);

	// Get result matrix
	for( uintMatrixIndex bc=0; bc < dim_indices.at(0).size(); ++bc ){
		for( uintMatrixIndex q=0; q < dim_indices.at(1).size(); ++q ){
			for( uintMatrixIndex pos=0; pos < dim_indices.at(2).size(); ++pos ){
				for( uintMatrixIndex num=0; num < dim_indices.at(3).size(); ++num ){
					for( uintMatrixIndex er=0; er < dim_indices.at(4).size(); ++er ){
						estimated_counts[dim_indices.at(0).at(bc)][dim_indices.at(1).at(q)][dim_indices.at(2).at(pos)][dim_indices.at(3).at(num)][dim_indices.at(4).at(er)] = estimate.base_call_.at(kTemplateSegment).at(0).at(0).at(0).estimates_.GetMatrixElement({bc, q, pos, num, er}) * tot_counts;
					}
				}
			}
		}
	}

	ShrinkVect(estimated_counts);
}

void ProbabilityEstimatesTest::IterativeProportionalFittingDomError(const vector<Vect<Vect<uintMatrixCount>>> &margins){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataDomError(stats, margins);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFDominantError, kTemplateSegment, 0, 0, 0, 0, kMaxIterations, kPrecisionAim);
}

void ProbabilityEstimatesTest::GetIPFResultDomError( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts, vector<Vect<Vect<uintMatrixCount>>> &margins ){
	estimated_counts.Clear();

	auto tot_counts = SumVect(margins.at(0)); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uintMatrixIndex>, 4 > &dim_indices(estimate.dom_error_.at(0).at(0).at(0).dim_indices_);

	// Get result matrix
	for( uintMatrixIndex dist=0; dist < dim_indices.at(1).size(); ++dist ){
		for( uintMatrixIndex dom_err=0; dom_err < dim_indices.at(0).size(); ++dom_err ){
			for( uintMatrixIndex gc=0; gc < dim_indices.at(2).size(); ++gc ){
				for( uintMatrixIndex start_rate=0; start_rate < dim_indices.at(3).size(); ++start_rate ){
					estimated_counts[dim_indices.at(0).at(dom_err)][dim_indices.at(1).at(dist)][dim_indices.at(2).at(gc)][dim_indices.at(3).at(start_rate)][0] = estimate.dom_error_.at(0).at(0).at(0).estimates_.GetMatrixElement({dom_err, dist, gc, start_rate}) * tot_counts;
				}
			}
		}
	}

	ShrinkVect(estimated_counts);
}

void ProbabilityEstimatesTest::SetMarginDefErrorRate(vector< pair<bool, bool> > &margin_def){
	// {use margin_quality_position, flip dimensions compared to dimension order}
	margin_def.push_back({false, true}); // dist, er
	margin_def.push_back({false, true}); // gc, er
	margin_def.push_back({false, false}); // dist, gc
	margin_def.push_back({false, true}); // sr, er
	margin_def.push_back({false, false}); // dist, sr
	margin_def.push_back({false, false}); // gc, sr
}

void ProbabilityEstimatesTest::IterativeProportionalFittingErrorRate(const vector<Vect<Vect<uintMatrixCount>>> &margins){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataErrorRate(stats, margins);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFErrorRate, kTemplateSegment, 0, 0, 0, 0, kMaxIterations, kPrecisionAim);
}

void ProbabilityEstimatesTest::GetIPFResultErrorRate( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts, vector<Vect<Vect<uintMatrixCount>>> &margins ){
	estimated_counts.Clear();

	auto tot_counts = SumVect(margins.at(0)); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uintMatrixIndex>, 4 > &dim_indices(estimate.error_rate_.at(0).at(0).dim_indices_);

	// Get result matrix
	for( uintMatrixIndex dist=0; dist < dim_indices.at(1).size(); ++dist ){
		for( uintMatrixIndex er=0; er < dim_indices.at(0).size(); ++er ){
			for( uintMatrixIndex gc=0; gc < dim_indices.at(2).size(); ++gc ){
				for( uintMatrixIndex start_rate=0; start_rate < dim_indices.at(3).size(); ++start_rate ){
					estimated_counts[dim_indices.at(0).at(er)][dim_indices.at(1).at(dist)][dim_indices.at(2).at(gc)][dim_indices.at(3).at(start_rate)][0] = estimate.error_rate_.at(0).at(0).estimates_.GetMatrixElement({er, dist, gc, start_rate}) * tot_counts;
				}
			}
		}
	}

	ShrinkVect(estimated_counts);
}

namespace reseq{
	TEST_F(ProbabilityEstimatesTest, DataStorageSetUp){
		ProbabilityEstimatesSubClasses::DataStorage<4> data;

		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		uintMatrixCount n(0);
		for(auto i=0; i<3; ++i){
			if(1 != i){
				for(auto j=1; j<3; ++j){
					for(auto k=0; k<2; ++k){
						counts[i][j][k][0][0] = ++n*data.kMinCountsPer1dBin; // Multiply everything by kMinCountsPer1dBin so we are above minimum counts and bins are not combined
					}
				}
			}
		}

		array< vector<uintMatrixIndex>, 4 > dim_indices;
		array< vector<uintMatrixIndex>, 4 > initial_dim_indices;
		SetUpDataStorageErrorRate(data, dim_indices, initial_dim_indices, counts);

		TestDataStorageSetUpErrorRate(data, dim_indices, initial_dim_indices);

		// Move some stuff into a higher/lower bin, so that it will be shifted back because it does not reach minimum counts and values should remain the same
		for(auto i=0; i<3; ++i){
			if(1 != i){
				for(auto k=0; k<2; ++k){
					counts.at(i)[0][k][0][0] = 24;
					counts.at(i).at(1).at(k).at(0).at(0) -= 24;

					counts.at(i)[3][k][0][0] = 24;
					counts.at(i).at(2).at(k).at(0).at(0) -= 24;
				}
			}
		}

		SetUpDataStorageErrorRate(data, dim_indices, initial_dim_indices, counts);

		TestDataStorageSetUpErrorRate(data, dim_indices, initial_dim_indices);

		// Move some stuff in between bin 0 and 2 for index i, so that it will not be combined by the end combining, but by the second combining step, this slightly changes dim_indices and initial_dim_indices
		for(auto i=0; i<3; ++i){
			if(1 != i){
				for(auto k=0; k<2; ++k){
					counts.at(i)[3][k][0][0] = counts.at(i).at(2).at(k).at(0).at(0)+24;
					counts.at(i).at(2).at(k).at(0).at(0) = 24;
					counts.at(i).at(1).at(k).at(0).at(0) -= 24;
				}
			}
		}

		SetUpDataStorageErrorRate(data, dim_indices, initial_dim_indices, counts);

		TestDataStorageSetUpErrorRate(data, dim_indices, initial_dim_indices, true);
	}

	TEST_F(ProbabilityEstimatesTest, SupportFunctions){
		TestLogArrayCalcNormalize();
		TestDataStorageReduceAndExpand();
		TestDataStorageReduceAndExpand2();
		TestLogArrayCalcExpand();
		TestCombineDimIndices();
	}

	TEST_F(ProbabilityEstimatesTest, IterativeProportionalFittingStepWise){
		uintBaseCall base(1);
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,7635921665);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IPFStepWiseQual(base, margins, margin_quality_position);
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQuality){
		uintBaseCall base(0);
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,103741084);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		string save_file = test_dir+"saveTest.reseq.ipf";
		ASSERT_TRUE( test_.Save(save_file.c_str()) );

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultQual(test_, base, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_file.c_str()) );
		iterations = GetIterationsQual(test2, base);
		precision = GetPrecisionQual(test2, base);
		GetIPFResultQual(test2, base, estimated_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "loading" );

		EXPECT_EQ( 0, remove(save_file.c_str()) ) << "Error deleting file: " << save_file << '\n';

		std::vector<double> prob;
		double prob_sum;
		test_.Quality(0, 0, base).Draw(prob,prob_sum,{200, 200, 200, 200},0.5); // Check for a crash if random values are queried (tile_id=0 must exist, otherwise crash is expected)
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithOffset){
		uintBaseCall base(2);
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,213576235);
		SetOffset(counts, 37, 0, 37, 0, 0);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultQual(test_, base, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithZeroMarginValues){
		uintBaseCall base(3);
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,3290428530);
		SetOffset(counts, 37, 0, 37, 0, 0);
		counts.at(38).at(2).Clear(); // Set quality 38, error rate 2 margin to zero

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultQual(test_, base, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithZeroOverAllMargins){
		uintBaseCall base(0);
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,9527958903);
		SetOffset(counts, 37, 0, 37, 0, 0);
		counts.at(38).Clear(); // Set quality 38 to zero for all margins

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultQual(test_, base, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforBaseCalls){
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,5,3,5,3,5,768293592);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		// {use margin_quality_position, flip dimensions compared to dimension order}
		vector< pair<bool, bool> > margin_def;
		margin_def.push_back({false, false}); // bc, qual
		margin_def.push_back({false, false}); // bc, pos
		margin_def.push_back({true, true}); // pos, qual
		margin_def.push_back({false, false}); // bc, num
		margin_def.push_back({false, true}); // num, qual
		margin_def.push_back({false, true}); // num, pos
		margin_def.push_back({false, false}); // bc, er
		margin_def.push_back({false, true}); // er, qual
		margin_def.push_back({false, false}); // pos, er
		margin_def.push_back({false, false}); // num, er
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingBaseCall(margins, margin_quality_position);
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		string save_file = test_dir+"saveTest.reseq.ipf";
		ASSERT_TRUE( test_.Save(save_file.c_str()) );

		auto iterations = GetIterationsBaseCall(test_);
		auto precision = GetPrecisionBaseCall(test_);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultBaseCall(test_, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_file.c_str()) );
		iterations = GetIterationsBaseCall(test2);
		precision = GetPrecisionBaseCall(test2);
		GetIPFResultBaseCall(test2, estimated_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "loading" );

		EXPECT_EQ( 0, remove(save_file.c_str()) ) << "Error deleting file: " << save_file << '\n';
	}

	TEST_F(ProbabilityEstimatesTest, IPFforDominantError){
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,5,3,3,5,1,793256239);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		// {use margin_quality_position, flip dimensions compared to dimension order}
		vector< pair<bool, bool> > margin_def;
		margin_def.push_back({false, true}); // dist, domErr
		margin_def.push_back({false, true}); // gc, domErr
		margin_def.push_back({false, false}); // dist, gc
		margin_def.push_back({false, true}); // sr, domErr
		margin_def.push_back({false, false}); // dist, sr
		margin_def.push_back({false, false}); // gc, sr
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingDomError(margins);
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		string save_file = test_dir+"saveTest.reseq.ipf";
		ASSERT_TRUE( test_.Save(save_file.c_str()) );

		auto iterations = GetIterationsDomError(test_);
		auto precision = GetPrecisionDomError(test_);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultDomError(test_, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_file.c_str()) );
		iterations = GetIterationsDomError(test2);
		precision = GetPrecisionDomError(test2);
		GetIPFResultDomError(test2, estimated_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "loading" );

		EXPECT_EQ( 0, remove(save_file.c_str()) ) << "Error deleting file: " << save_file << '\n';
	}

	TEST_F(ProbabilityEstimatesTest, IPFforErrorRate){
		Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > counts;
		RandomFillingTestCounts(counts,5,3,5,5,1,3465435782);

		vector<Vect<Vect<uintMatrixCount>>> margins;
		Vect<SeqQualityStats<uintMatrixCount>> margin_quality_position;
		// {use margin_quality_position, flip dimensions compared to dimension order}
		vector< pair<bool, bool> > margin_def;
		SetMarginDefErrorRate(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingErrorRate(margins);
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		string save_file = test_dir+"saveTest.reseq.ipf";
		ASSERT_TRUE( test_.Save(save_file.c_str()) );

		auto iterations = GetIterationsErrorRate(test_);
		auto precision = GetPrecisionErrorRate(test_);

		Vect< Vect< Vect< Vect< Vect<double> > > > > estimated_counts;
		GetIPFResultErrorRate(test_, estimated_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_file.c_str()) );
		iterations = GetIterationsErrorRate(test2);
		precision = GetPrecisionErrorRate(test2);
		GetIPFResultErrorRate(test2, estimated_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, estimated_counts, "loading" );

		EXPECT_EQ( 0, remove(save_file.c_str()) ) << "Error deleting file: " << save_file << '\n';
	}
}
