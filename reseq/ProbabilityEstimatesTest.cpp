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

#include "utilities.h"
using reseq::utilities::Divide;

inline double ProbabilityEstimatesTest::CalculateAbsolutPrecision(uint64_t correct_margin, double marginal_sum, double precision){
	// Precision is defined to be always >1 so this has to be reverted first
	double sum_factor = correct_margin > marginal_sum ? precision : 1.0/precision;
	// We are interested in the absolut difference and not absolut value
	sum_factor = fabs( sum_factor - 1);
	// Precision is not evaluated for correct_margings equal to 0, because it would be infinity, as we still check for this cases everything considered a perfect fit (<epsilon) is accepted here
	return max(epsilon_, marginal_sum*sum_factor);
}

void ProbabilityEstimatesTest::SetUpDataQual(DataStats &stats, uint16_t base, const vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position){
	stats.qualities_.base_quality_for_error_rate_per_tile_reference_[template_segment_][base][0] = margins[0];
	stats.qualities_.base_quality_for_preceding_quality_per_tile_reference_[template_segment_][base][0] = margins[1];
	stats.qualities_.preceding_quality_for_error_rate_per_tile_reference_[template_segment_][base][0] = margins[2];
	stats.qualities_.base_quality_stats_per_tile_reference_[template_segment_][base][0] = margin_quality_position;
	stats.qualities_.error_rate_for_position_per_tile_reference_[template_segment_][base][0] = margins[4];
	stats.qualities_.preceding_quality_for_position_per_tile_reference_[template_segment_][base][0] = margins[5];
	stats.qualities_.base_quality_for_sequence_quality_per_tile_reference_[template_segment_][base][0] = margins[6];
	stats.qualities_.sequence_quality_for_error_rate_per_tile_reference_[template_segment_][base][0] = margins[7];
	stats.qualities_.preceding_quality_for_sequence_quality_per_tile_reference_[template_segment_][base][0] = margins[8];
	stats.qualities_.sequence_quality_for_position_per_tile_reference_[template_segment_][base][0] = margins[9];
}

void ProbabilityEstimatesTest::SetUpDataBaseCall(DataStats &stats, const vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position){
	stats.errors_.called_bases_by_base_quality_per_tile_[template_segment_][0][0][0] = margins[0];
	stats.errors_.called_bases_by_position_per_tile_[template_segment_][0][0][0] = margins[1];
	stats.qualities_.base_quality_stats_per_tile_per_error_reference_[template_segment_][0][0][0] = margin_quality_position;
	stats.errors_.called_bases_by_error_num_per_tile_[template_segment_][0][0][0] = margins[3];
	stats.errors_.error_num_by_quality_per_tile_[template_segment_][0][0][0] = margins[4];
	stats.errors_.error_num_by_position_per_tile_[template_segment_][0][0][0] = margins[5];
	stats.errors_.called_bases_by_error_rate_per_tile_[template_segment_][0][0][0] = margins[6];
	stats.qualities_.base_quality_for_error_rate_per_tile_per_error_reference_[template_segment_][0][0][0] = margins[7];
	stats.qualities_.error_rate_for_position_per_tile_per_error_reference_[template_segment_][0][0][0] = margins[8];
	stats.errors_.error_num_by_error_rate_per_tile_[template_segment_][0][0][0] = margins[9];
}

void ProbabilityEstimatesTest::SetUpDataDomError(DataStats &stats, const vector<Vect<Vect<uint64_t>>> &margins){
	stats.coverage_.dominant_errors_by_distance_[0][0][0] = margins[0];
	stats.coverage_.dominant_errors_by_gc_[0][0][0] = margins[1];
	stats.coverage_.gc_by_distance_de_[0][0][0] = margins[2];
}

void ProbabilityEstimatesTest::SetUpDataErrorRate(DataStats &stats, const vector<Vect<Vect<uint64_t>>> &margins){
	stats.coverage_.error_rates_by_distance_[0][0] = margins[0];
	stats.coverage_.error_rates_by_gc_[0][0] = margins[1];
	stats.coverage_.gc_by_distance_er_[0][0] = margins[2];
}

void ProbabilityEstimatesTest::SetUp(){
	BasicTestClass::SetUp();

	ResetProbabilityEstimates();
}

void ProbabilityEstimatesTest::TestLogArrayCalcNormalize(){
	LogArrayCalc<3> test, test_normalized;

	// Set dimensions
	for( uint16_t n=3; n--;){
		test.dim_size_[n] = 3;
		test_normalized.dim_size_[n] = 3;
	}

	// Randomly fill test
	mt19937_64 rgen;
	rgen.seed(20180315);
	uniform_real_distribution<double> dist(0.5, 1.5);

	for( uint16_t n=3;n--;){
		test.dim2_[n].resize(9);
		for( uint16_t p1=3; p1--; ){
			for( uint16_t p2=3; p2--; ){
				if( 1 == p1 && 1 == p2 ){
					test.dim2_[n][4] = 0.0; // Introduce zeros to see if the normalization can cope with that
				}
				else{
					test.dim2_[n][p1*3+p2] = dist(rgen);
				}

				// Introduce factors that cancel each other out that are close to the numeric limit to see whether that is a problem for the normalization
				if(2==n){
					test.dim2_[n][p1*3+p2] *= numeric_limits<double>::max()/2.0;
				}
				else if(1==n){
					test.dim2_[n][p1*3+p2] *= 2.0/numeric_limits<double>::max();
				}
			}
		}
	}

	// Copy values and renormalize new LogArray
	for( uint16_t n=3;n--;){
		test_normalized.dim2_[n].resize(9);
		for( uint16_t p1=3; p1--; ){
			for( uint16_t p2=3; p2--; ){
				test_normalized.dim2_[n][p1*3+p2] = test.dim2_[n][p1*3+p2];
			}
		}
	}
	test_normalized.Normalize();

	// Now check if values got as small as they are supposed to be (Nothing can be higher than the max of the real distribution 1.5)
	for( uint16_t n=3;n--;){
		for( uint16_t p1=3; p1--; ){
			for( uint16_t p2=3; p2--; ){
				EXPECT_TRUE( 1.6 > test_normalized.dim2_[n][p1*3+p2] ); // Chose 1.6, that is good enough and very generous towards precision issues
			}
		}
	}

	// Now check if the results are still the same
	for( uint16_t p1=3; p1--; ){
		for( uint16_t p2=3; p2--; ){
			for( uint16_t p3=3; p3--; ){
				// Double precision cannot be expected after the logs used
				EXPECT_FLOAT_EQ( test.dim2_[2][p3*3+p2]*test.dim2_[1][p3*3+p1]*test.dim2_[0][p2*3+p1],
						test_normalized.dim2_[2][p3*3+p2]*test_normalized.dim2_[1][p3*3+p1]*test_normalized.dim2_[0][p2*3+p1] );
			}
		}
	}
}

void ProbabilityEstimatesTest::TestDataStorageReduceAndExpand(){
	// Prepare test data
	DataStorage<3> test;
	for(uint16_t n=3; n--; ){
		test.dim_size_[n] = 3+2*n;
	}
	test.data_[0].resize(5*3, 0.0); // 1,0
	test.data_[1].resize(7*3, 0.0); // 2,0
	test.data_[2].resize(7*5, 0.0); // 2,1

	double value;
	uint16_t red_i, red_j, red_k;
	for(auto i=test.dim_size_[0]; i--; ){
		for(auto j=test.dim_size_[1]; j--; ){
			for(auto k=test.dim_size_[2]; k--; ){
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

				if( !(0==red_i && 0==red_j || 1==red_i && 1==red_k || 2==red_j && 2==red_k) ){
					value = (red_i+1)*(red_j+4)*(red_k+7);

					test.data_[0][j*test.dim_size_[0]+i] += value;
					test.data_[1][k*test.dim_size_[0]+i] += value;
					test.data_[2][k*test.dim_size_[1]+j] += value;
				}
			}
		}
	}

	// Test Reduce function
	DataStorage<3> reduced_data;
	std::array<std::vector<uint16_t>, 3> dim_indices_reduced, dim_indices_count;

	test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 3);

	// Check for correctness
	for(uint16_t n=3; n--; ){
		EXPECT_EQ(3+2*n, dim_indices_reduced[n].size()) << "dim_indices_reduced[" << n << "].size() wrong";
		for(auto i=dim_indices_reduced[n].size(); i--; ){
			// Index range is from 0 to 2
			EXPECT_TRUE(dim_indices_reduced[n][i] < 3) << "dim_indices_reduced[" << n << "][" << i << "] invalid index";
		}
		for(auto i=dim_indices_reduced[n].size(); --i; ){
			// First index has to be different then all other indeces
			EXPECT_TRUE(dim_indices_reduced[n][i] != dim_indices_reduced[n][0]) << "dim_indices_reduced[" << n << "][" << i << "] is equal to first index";
		}
		for(auto i=dim_indices_reduced[n].size()-1; i--; ){
			if(i < 4){
				// Last index has to be different then all other indeces
				EXPECT_TRUE(dim_indices_reduced[n][i] != dim_indices_reduced[n][dim_indices_reduced[n].size()-1]) << "dim_indices_reduced[" << n << "][" << i << "] is equal to last index";
			}
			else{
				// Special case for dimension 3, where the last index exists three times
				EXPECT_TRUE(dim_indices_reduced[n][i] == dim_indices_reduced[n][dim_indices_reduced[n].size()-1]) << "dim_indices_reduced[" << n << "][" << i << "] is not equal to last index";
			}
		}
		for(auto i=dim_indices_reduced[n].size()-(2==n?4:2); --i; ){
			// Middle indeces have to be equal
			EXPECT_TRUE(dim_indices_reduced[n][i] == dim_indices_reduced[n][dim_indices_reduced[n].size()-(2==n?4:2)]) << "dim_indices_reduced[" << n << "][" << i << "] is not equal to middle indeces";
		}

		EXPECT_EQ(3, dim_indices_count[n].size()) << "dim_indices_count[" << n << "].size() wrong";
		EXPECT_EQ(1, dim_indices_count[n][ dim_indices_reduced[n][0] ]) << "dim_indices_count[" << n << "][" << dim_indices_reduced[n][0] << "] wrong";
		EXPECT_EQ((0==n?1:3), dim_indices_count[n][ dim_indices_reduced[n][1] ]) << "dim_indices_count[" << n << "][" << dim_indices_reduced[n][1] << "] wrong";
		EXPECT_EQ((2==n?3:1), dim_indices_count[n][ dim_indices_reduced[n][dim_indices_reduced[n].size()-1] ]) << "dim_indices_count[" << n << "][" << dim_indices_reduced[n][dim_indices_reduced[n].size()-1] << "] wrong";
	}

	uint16_t dim_a(3), dim_b(2);
	for(uint16_t n=3; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = 2;
		}

		for(auto i=dim_indices_reduced[dim_a].size(); i--; ){
			for(auto j=dim_indices_reduced[dim_b].size(); j--; ){
				EXPECT_DOUBLE_EQ( reduced_data.data_[n][ dim_indices_reduced[dim_a][i]*dim_indices_count[dim_b].size() + dim_indices_reduced[dim_b][j] ],
						test.data_[n][ i*dim_indices_reduced[dim_b].size() + j] * dim_indices_count[dim_a][ dim_indices_reduced[dim_a][i] ] * dim_indices_count[dim_b][ dim_indices_reduced[dim_b][j] ] );
			}
		}
	}

	// Test Expand function
	std::array<std::vector<uint16_t>, 3> expansion_indices, expansion_count;
	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uint16_t n=3; n--; ){
		EXPECT_EQ((2==n?6:3), dim_indices_count[n].size()) << "dim_indices_count[" << n << "].size() wrong";
		if(0==n){
			for( auto count : dim_indices_count[n]){
				EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after first expansion";
			}
			for( auto ind=dim_indices_reduced[n].size(); ind--; ){
				EXPECT_EQ(ind, dim_indices_reduced[n][ind]) << "dim_indices_reduced[" << n << "] did not go back to original";
			}
		}
		else if(2==n){
			bool count_two_left(true);
			for( auto count : dim_indices_count[n]){
				EXPECT_TRUE(1 == count || count_two_left && 2 == count) << "dim_indices_count[" << n << "] has wrong counts after first expansion";
				if(2 == count){
					count_two_left = false;
				}
			}
			std::vector<uint16_t> tmp_counts = dim_indices_count[n];
			for( auto ind : dim_indices_reduced[n]){
				EXPECT_TRUE( tmp_counts[ind]-- )  << "dim_indices_count[" << n << "] does not fit with dim_indices_reduced after first expansion";
			}
		}

		EXPECT_EQ( 1, expansion_count[n][expansion_indices[n][dim_indices_reduced[n][0]]] ) << "expansion_count[" << n << "] wrong";
		EXPECT_EQ( (0==n?1:3), expansion_count[n][expansion_indices[n][dim_indices_reduced[n][1]]] ) << "expansion_count[" << n << "] wrong";
		EXPECT_EQ( (2==n?3:1), expansion_count[n][expansion_indices[n][dim_indices_reduced[n][dim_indices_reduced[n].size()-1]]] ) << "expansion_count[" << n << "] wrong";
	}

	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uint16_t n=3; n--; ){
		EXPECT_EQ((2==n?6:3+2*n), dim_indices_count[n].size()) << "dim_indices_count[" << n << "].size() wrong";
		if(2!=n){
			for( auto count : dim_indices_count[n]){
				EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after second expansion";
			}
			for( auto ind=dim_indices_reduced[n].size(); ind--; ){
				EXPECT_EQ(ind, dim_indices_reduced[n][ind]) << "dim_indices_reduced[" << n << "] did not go back to original";
			}
		}
	}

	EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
	for(uint16_t n=3; n--; ){
		EXPECT_EQ(3+2*n, dim_indices_count[n].size()) << "dim_indices_count[" << n << "].size() wrong";
		for( auto count : dim_indices_count[n]){
			EXPECT_EQ(1, count) << "dim_indices_count[" << n << "] has wrong counts after third expansion";
		}
		for( auto ind=dim_indices_reduced[n].size(); ind--; ){
			EXPECT_EQ(ind, dim_indices_reduced[n][ind]) << "dim_indices_reduced[" << n << "] did not go back to original";
		}
	}

	EXPECT_FALSE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";
}

void ProbabilityEstimatesTest::TestDataStorageReduceAndExpand2(){
	// Prepare test data
	DataStorage<3> test;
	test.dim_size_[0] = 1;
	test.dim_size_[1] = 1;
	test.dim_size_[2] = 5;
	test.data_[0].resize(1, 0.0); // 1,0
	test.data_[1].resize(5, 0.0); // 2,0
	test.data_[2].resize(5, 0.0); // 2,1

	double value;

	for(auto k=test.dim_size_[2]; k--; ){
		if(k < 2){
			value=200+k;
		}
		else{
			value=96+2*k;
		}

		test.data_[0][0] += value;
		test.data_[1][k] += value;
		test.data_[2][k] += value;
	}

	// Test Reduce function
	std::array<std::vector<uint16_t>, 3> dim_indices_reduced, dim_indices_count;

	for( uint16_t test_case=0; test_case<2; ++test_case){
		for( uint16_t n = 3; n--; ){
			dim_indices_count[n].clear();
			dim_indices_reduced[n].clear();
		}

		if(0 == test_case){
			DataStorage<3> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 2);
		}
		else{
			// Test Expand function
			DataStorage<3> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 1);

			std::array<std::vector<uint16_t>, 3> expansion_indices, expansion_count;

			EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 1) ) << "Expansion returns wrong value";

			for( uint16_t n = 3; n--; ){
				EXPECT_EQ( 1, expansion_count[n].size() ) << "expansion_count[" << n << "].size() wrong";
				EXPECT_EQ( (2==n?5:1), expansion_count[n][0] ) << "expansion_count[" << n << "][0] wrong";
				EXPECT_EQ( (2==n?2:1), expansion_indices[n].size() ) << "expansion_indices[" << n << "].size() wrong";
				for(uint16_t i=0; i < expansion_indices[n].size(); ++i ){
					EXPECT_EQ( 0, expansion_indices[n][i] ) << "expansion_count[" << n << "][" << i << "] wrong";
				}
			}
		}

		EXPECT_EQ(1, dim_indices_count[0][0]) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(1, dim_indices_count[1][0]) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(0, dim_indices_reduced[0][0]) << "Error in simple dimension in case " << test_case;
		EXPECT_EQ(0, dim_indices_reduced[1][0]) << "Error in simple dimension in case " << test_case;

		if(2 == dim_indices_count[2][0]){
			EXPECT_EQ(3, dim_indices_count[2][1]) << "dim_indices_count[2] wrong: " << dim_indices_count[2][0] << ' ' << dim_indices_count[2][1] << " in case " << test_case;
			for(uint16_t i=0; i<2; ++i){
				EXPECT_EQ(0, dim_indices_reduced[2][i]) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
			for(uint16_t i=2; i<5; ++i){
				EXPECT_EQ(1, dim_indices_reduced[2][i]) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
		}
		else{
			EXPECT_EQ(3, dim_indices_count[2][0]) << "dim_indices_count[2] wrong: " << dim_indices_count[2][0] << ' ' << dim_indices_count[2][1] << " in case " << test_case;
			EXPECT_EQ(2, dim_indices_count[2][1]) << "dim_indices_count[2] wrong: " << dim_indices_count[2][0] << ' ' << dim_indices_count[2][1] << " in case " << test_case;
			for(uint16_t i=0; i<2; ++i){
				EXPECT_EQ(1, dim_indices_reduced[2][i]) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
			for(uint16_t i=2; i<5; ++i){
				EXPECT_EQ(2, dim_indices_reduced[2][i]) << "dim_indices_reduced[2] wrong in case " << test_case;
			}
		}
	}

	for( uint16_t test_case=0; test_case<2; ++test_case){
		for( uint16_t n = 3; n--; ){
			dim_indices_count[n].clear();
			dim_indices_reduced[n].clear();
		}

		if(0 == test_case){
			DataStorage<3> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 4);
		}
		else{
			// Test Expand function
			DataStorage<3> reduced_data;

			test.Reduce(reduced_data, dim_indices_reduced, dim_indices_count, 1);

			std::array<std::vector<uint16_t>, 3> expansion_indices, expansion_count;

			EXPECT_TRUE( test.Expand(reduced_data, dim_indices_reduced, dim_indices_count, expansion_indices, expansion_count, 2) ) << "Expansion returns wrong value";

			for( uint16_t n = 3; n--; ){
				EXPECT_EQ( 1, expansion_count[n].size() ) << "expansion_count[" << n << "].size() wrong";
				EXPECT_EQ( (2==n?5:1), expansion_count[n][0] ) << "expansion_count[" << n << "][0] wrong";
				EXPECT_EQ( (2==n?4:1), expansion_indices[n].size() ) << "expansion_indices[" << n << "].size() wrong";
				for(uint16_t i=0; i < expansion_indices[n].size(); ++i ){
					EXPECT_EQ( 0, expansion_indices[n][i] ) << "expansion_count[" << n << "][" << i << "] wrong";
				}
			}
		}

		EXPECT_EQ(dim_indices_reduced[2][0], dim_indices_reduced[2][1]) << "Max difference for the separation doesn't seem to work in case " << test_case;
		EXPECT_EQ(2, dim_indices_count[2][ dim_indices_reduced[2][1] ]) << "Max difference for the separation doesn't seem to work in case " << test_case;

		bool count_two_left(true);
		for( auto count : dim_indices_count[2]){
			EXPECT_TRUE(1 == count || count_two_left && 2 == count) << "dim_indices_count[2] has wrong counts in case " << test_case;
			if(2 == count){
				count_two_left = false;
			}
		}
		std::vector<uint16_t> tmp_counts = dim_indices_count[2];
		for( auto ind : dim_indices_reduced[2]){
			EXPECT_TRUE( tmp_counts[ind]-- )  << "dim_indices_count[2] does not fit with dim_indices_reduced in case " << test_case;
		}
	}
}

void ProbabilityEstimatesTest::TestLogArrayCalcExpand(){
	LogArrayCalc<3> test;

	// Set dimensions
	for( uint16_t n=3; n--;){
		test.dim_size_[n] = 3;

	}

	// Randomly fill test
	mt19937_64 rgen;
	rgen.seed(20180320);
	uniform_real_distribution<double> dist(0.5, 1.5);

	for( uint16_t n=3;n--;){
		test.dim2_[n].resize(9);
		for( uint16_t p1=3; p1--; ){
			for( uint16_t p2=3; p2--; ){
				if( n == p1 && n == p2 ){
					test.dim2_[n][4] = 0.0; // Introduce zeros to see if the extension can cope with that
				}
				else{
					test.dim2_[n][p1*3+p2] = dist(rgen);
				}
			}
		}
	}

	// Calculate marginal sums for later comparison
	std::array<std::vector<double>, 3> marginal_sums;
	for( uint16_t n=3; n--;){
		marginal_sums[n].resize(9, 0.0);
	}

	double value;
	for( uint16_t p1=test.dim_size_[0]; p1--; ){
		for( uint16_t p2=test.dim_size_[1]; p2--; ){
			for( uint16_t p3=test.dim_size_[2]; p3--; ){
				// Double precision cannot be expected after the logs used
				value = test.dim2_[2][p3*3+p2]*test.dim2_[1][p3*3+p1]*test.dim2_[0][p2*3+p1];

				marginal_sums[0][p2*test.dim_size_[0]+p1] += value;
				marginal_sums[1][p3*test.dim_size_[0]+p1] += value;
				marginal_sums[2][p3*test.dim_size_[1]+p2] += value;
			}
		}
	}

	// Expand LogArray
	std::array<std::vector<uint16_t>, 3> dim_indices_reduced, dim_indices_count;
	dim_indices_reduced[0] = {0,1,2};
	dim_indices_reduced[1] = {0,1,1,2};
	dim_indices_reduced[2] = {0,1,1,1,2};
	dim_indices_count[0] = {1,1,1};
	dim_indices_count[1] = {1,2,1};
	dim_indices_count[2] = {1,3,1};

	test.Expand(dim_indices_reduced, dim_indices_count);

	// Calculate expanded marginal sums
	std::array<std::vector<double>, 3> expanded_marginal_sums;
	expanded_marginal_sums[0].resize(test.dim_size_[1]*test.dim_size_[0]);
	expanded_marginal_sums[1].resize(test.dim_size_[2]*test.dim_size_[0]);
	expanded_marginal_sums[2].resize(test.dim_size_[2]*test.dim_size_[1]);

	for( uint16_t p1=test.dim_size_[0]; p1--; ){
		for( uint16_t p2=test.dim_size_[1]; p2--; ){
			for( uint16_t p3=test.dim_size_[2]; p3--; ){
				// Double precision cannot be expected after the logs used
				value = test.dim2_[2][p3*4+p2]*test.dim2_[1][p3*3+p1]*test.dim2_[0][p2*3+p1];

				expanded_marginal_sums[0][p2*test.dim_size_[0]+p1] += value;
				expanded_marginal_sums[1][p3*test.dim_size_[0]+p1] += value;
				expanded_marginal_sums[2][p3*test.dim_size_[1]+p2] += value;
			}
		}
	}

	// Compare marginal sums
	uint16_t dim_a(3), dim_b(2);
	for(uint16_t n=3; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = 2;
		}

		for(auto i=dim_indices_reduced[dim_a].size(); i--; ){
			for(auto j=dim_indices_reduced[dim_b].size(); j--; ){
				EXPECT_DOUBLE_EQ( marginal_sums[n][ dim_indices_reduced[dim_a][i]*dim_indices_count[dim_b].size() + dim_indices_reduced[dim_b][j] ],
						expanded_marginal_sums[n][ i*dim_indices_reduced[dim_b].size() + j] * dim_indices_count[dim_a][ dim_indices_reduced[dim_a][i] ] * dim_indices_count[dim_b][ dim_indices_reduced[dim_b][j] ] )
						<< '[' << n << "][" << i << "][" << j << ']';
			}
		}
	}
}

void ProbabilityEstimatesTest::TestCombineDimIndices(){
	array< vector<uint16_t>, 2 > dim_indeces1, dim_indeces2;

	for(uint16_t cur_dim=2; cur_dim--; ){
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

	for(uint16_t cur_dim=2; cur_dim--; ){
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
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts,
		const uint16_t dim1,
		const uint16_t dim2,
		const uint16_t dim3,
		const uint16_t dim4,
		const uint16_t dim5,
		const uint64_t seed){
	mt19937_64 rgen;
	rgen.seed(seed);
	uniform_int_distribution<uint16_t> dist(0, 10);
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
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts,
		const uint16_t dim1,
		const uint16_t dim2,
		const uint16_t dim3,
		const uint16_t dim4,
		const uint16_t dim5){
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
		uint32_t iterations,
		double precision,
		const vector<Vect<Vect<uint64_t>>> &margins,
		const Vect<SeqQualityStats<uint64_t>> &margin_quality_position,
		const vector< pair<bool, bool> > &margin_def,
		const Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts,
		string context ){
	EXPECT_LE( iterations, max_iterations_ ) << "Maximum iterations exceeded\n";
	EXPECT_LT( precision, precision_aim_ ) << "Precision aim not achieved\n";

	// Get margins
	vector<Vect<Vect<double>>> comp_margins;
	Vect<SeqQualityStats<double>> comp_margin_quality_position;
	CalculateMargins(margin_def, comp_margins, comp_margin_quality_position, comp_counts);

	// Check all defined margins
	double marginal_sum;
	uint16_t dim1 = 0;
	uint16_t dim2 = 1;
	for(uint16_t i=0; i < margins.size(); ++i){
		for( auto pos1=margins[i].from(); pos1 < margins[i].to(); ++pos1 ){
			for( auto pos2=margins[i][pos1].from(); pos2 < margins[i][pos1].to(); ++pos2 ){
				EXPECT_NEAR( margins[i][pos1][pos2], comp_margins[i][pos1][pos2], CalculateAbsolutPrecision(margins[i][pos1][pos2], comp_margins[i][pos1][pos2], precision) ) << "margins[" << i << "][" << pos1 << "][" << pos2 << "] does not match within precision for " << context << '\n';
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

void ProbabilityEstimatesTest::SetUpDataStorageErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, array< vector<uint64_t>, 3 > &dim_indices, array< vector<uint16_t>, 3 > &initial_dim_indices, const Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts){
	vector<Vect<Vect<uint64_t>>> margins;
	Vect<SeqQualityStats<uint64_t>> margin_quality_position;
	vector< pair<bool, bool> > margin_def;
	SetMarginDefErrorRate(margin_def);
	CalculateMargins(margin_def, margins, margin_quality_position, counts);

	DataStats stats(NULL);
	SetUpDataErrorRate(stats, margins);

	array< pair<const Vect<Vect<uint64_t>> *, bool>, 3 > margin_def2;
	test_.DefineMarginsErrorRate(stats, margin_def2, 0, 0);

	mutex print_mutex;
	data.SetUp(margin_def2, NULL, dim_indices, initial_dim_indices, print_mutex);
}

void ProbabilityEstimatesTest::TestDataStorageSetUpErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, array< vector<uint64_t>, 3 > &dim_indices, array< vector<uint16_t>, 3 > &initial_dim_indices, bool middle_bin_inserted){

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
	EXPECT_DOUBLE_EQ( 15.0/36.0 , data.data_.at(0)[3] );

	EXPECT_EQ( 4 , data.data_.at(1).size() );
	EXPECT_DOUBLE_EQ( 4.0/36.0 , data.data_.at(1).at(0) );
	EXPECT_DOUBLE_EQ( 12.0/36.0 , data.data_.at(1).at(1) );
	EXPECT_DOUBLE_EQ( 6.0/36.0 , data.data_.at(1).at(2) );
	EXPECT_DOUBLE_EQ( 14.0/36.0 , data.data_.at(1)[3] );

	EXPECT_EQ( 4 , data.data_.at(2).size() );
	EXPECT_DOUBLE_EQ( 6.0/36.0 , data.data_.at(2).at(0) );
	EXPECT_DOUBLE_EQ( 10.0/36.0 , data.data_.at(2).at(1) );
	EXPECT_DOUBLE_EQ( 8.0/36.0 , data.data_.at(2).at(2) );
	EXPECT_DOUBLE_EQ( 12.0/36.0 , data.data_.at(2)[3] );

	if(middle_bin_inserted){
		EXPECT_EQ( 0 , initial_dim_indices.at(1).at(0) ) << "initial_dim_indices[" << 0 << "] is wrong";
		EXPECT_EQ( 0 , initial_dim_indices.at(1).at(1) ) << "initial_dim_indices[" << 0 << "] is wrong";
		EXPECT_EQ( 1 , initial_dim_indices.at(1).at(2) ) << "initial_dim_indices[" << 0 << "] is wrong";

		for(uint16_t cur_dim=0; cur_dim < initial_dim_indices.size(); cur_dim+=2){
			EXPECT_EQ( 0 , initial_dim_indices.at(cur_dim).at(0) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
			EXPECT_EQ( 1 , initial_dim_indices.at(cur_dim).at(1) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
		}
	}
	else{
		for(uint16_t cur_dim=0; cur_dim < initial_dim_indices.size(); ++cur_dim){
			EXPECT_EQ( 0 , initial_dim_indices.at(cur_dim).at(0) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
			EXPECT_EQ( 1 , initial_dim_indices.at(cur_dim).at(1) ) << "initial_dim_indices[" << cur_dim << "] is wrong";
		}
	}
}

void ProbabilityEstimatesTest::IterativeProportionalFittingQual(uint16_t base, const vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataQual(stats, base, margins, margin_quality_position);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFQuality, template_segment_, 0, base, 0, 0, max_iterations_, precision_aim_);
}

void ProbabilityEstimatesTest::GetIPFResultQual( const ProbabilityEstimates &estimate, uint16_t base, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, vector<Vect<Vect<uint64_t>>> &margins ){
	comp_counts.Clear();

	auto tot_counts = SumVect(margins[0]); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uint64_t>, 5 > &dim_indices(estimate.quality_[template_segment_][0][base].dim_indices_);

	// Get result matrix
	for( uint32_t q=0; q < dim_indices[0].size(); ++q ){
		for( uint32_t er=0; er < dim_indices[4].size(); ++er ){
			for( uint32_t pq=0; pq < dim_indices[2].size(); ++pq ){
				for( uint32_t pos=0; pos < dim_indices[3].size(); ++pos ){
					for( uint32_t sq=0; sq < dim_indices[1].size(); ++sq ){
						comp_counts[dim_indices[0][q]][dim_indices[4][er]][dim_indices[2][pq]][dim_indices[3][pos]][dim_indices[1][sq]] = estimate.quality_[template_segment_][0][base].estimates_.GetMatrixElement({q, sq, pq, pos, er}) * tot_counts;
					}
				}
			}
		}
	}

	ShrinkVect(comp_counts);
}

void ProbabilityEstimatesTest::IPFStepWiseQual(uint16_t base, const vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position){
	DataStats stats(NULL);
	SetUpDataQual(stats, base, margins, margin_quality_position);

	std::array< std::pair<const Vect<Vect<uint64_t>> *, bool>, 10 > margin_defs;
	auto alternative_margin = test_.DefineMarginsQuality(stats, margin_defs, template_segment_, 0, base);

	array< vector<uint64_t>, 5 > dim_indices;
	array< vector<uint16_t>, 5 > initial_dim_indices;
	DataStorage<5> data;
	mutex print_mutex;
	data.SetUp(margin_defs, alternative_margin, dim_indices, initial_dim_indices, print_mutex);

	test_.quality_[template_segment_][0][base].estimates_.SetUp(dim_indices);

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

void ProbabilityEstimatesTest::IterativeProportionalFittingBaseCall( const vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position ){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataBaseCall(stats, margins, margin_quality_position);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFBaseCall, template_segment_, 0, 0, 0, 0, max_iterations_, precision_aim_);
}

void ProbabilityEstimatesTest::GetIPFResultBaseCall( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, vector<Vect<Vect<uint64_t>>> &margins ){
	comp_counts.Clear();

	auto tot_counts = SumVect(margins[0]); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uint64_t>, 5 > &dim_indices(estimate.base_call_[template_segment_][0][0][0].dim_indices_);

	// Get result matrix
	for( uint32_t bc=0; bc < dim_indices[0].size(); ++bc ){
		for( uint32_t q=0; q < dim_indices[1].size(); ++q ){
			for( uint32_t pos=0; pos < dim_indices[2].size(); ++pos ){
				for( uint32_t num=0; num < dim_indices[3].size(); ++num ){
					for( uint32_t er=0; er < dim_indices[4].size(); ++er ){
						comp_counts[dim_indices[0][bc]][dim_indices[1][q]][dim_indices[2][pos]][dim_indices[3][num]][dim_indices[4][er]] = estimate.base_call_[template_segment_][0][0][0].estimates_.GetMatrixElement({bc, q, pos, num, er}) * tot_counts;
					}
				}
			}
		}
	}

	ShrinkVect(comp_counts);
}

void ProbabilityEstimatesTest::IterativeProportionalFittingDomError(const vector<Vect<Vect<uint64_t>>> &margins){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataDomError(stats, margins);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFDominantError, template_segment_, 0, 0, 0, 0, max_iterations_, precision_aim_);
}

void ProbabilityEstimatesTest::GetIPFResultDomError( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, vector<Vect<Vect<uint64_t>>> &margins ){
	comp_counts.Clear();

	auto tot_counts = SumVect(margins[0]); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uint64_t>, 3 > &dim_indices(estimate.dom_error_[0][0][0].dim_indices_);

	// Get result matrix
	for( uint32_t dist=0; dist < dim_indices[1].size(); ++dist ){
		for( uint32_t dom_err=0; dom_err < dim_indices[0].size(); ++dom_err ){
			for( uint32_t gc=0; gc < dim_indices[2].size(); ++gc ){
				comp_counts[dim_indices[0][dom_err]][dim_indices[1][dist]][dim_indices[2][gc]][0][0] = estimate.dom_error_[0][0][0].estimates_.GetMatrixElement({dom_err, dist, gc}) * tot_counts;
			}
		}
	}

	ShrinkVect(comp_counts);
}

void ProbabilityEstimatesTest::SetMarginDefErrorRate(vector< pair<bool, bool> > &margin_def){
	// {use margin_quality_position, flip dimensions compared to dimension order}
	margin_def.push_back({false, true}); // dist, er
	margin_def.push_back({false, true}); // gc, er
	margin_def.push_back({false, false}); // dist, gc
}

void ProbabilityEstimatesTest::IterativeProportionalFittingErrorRate(const vector<Vect<Vect<uint64_t>>> &margins){
	// Run iterative proportional fitting
	DataStats stats(NULL);
	SetUpDataErrorRate(stats, margins);
	test_.IterativeProportionalFitting(stats, ProbabilityEstimates::kIPFErrorRate, template_segment_, 0, 0, 0, 0, max_iterations_, precision_aim_);
}

void ProbabilityEstimatesTest::GetIPFResultErrorRate( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, vector<Vect<Vect<uint64_t>>> &margins ){
	comp_counts.Clear();

	auto tot_counts = SumVect(margins[0]); // As the data conversion before the IPF normalizes the sum of the matrix to 1 we have to revert it here

	const std::array< std::vector<uint64_t>, 3 > &dim_indices(estimate.error_rate_[0][0].dim_indices_);

	// Get result matrix
	for( uint32_t dist=0; dist < dim_indices[1].size(); ++dist ){
		for( uint32_t er=0; er < dim_indices[0].size(); ++er ){
			for( uint32_t gc=0; gc < dim_indices[2].size(); ++gc ){
				comp_counts[dim_indices[0][er]][dim_indices[1][dist]][dim_indices[2][gc]][0][0] = estimate.error_rate_[0][0].estimates_.GetMatrixElement({er, dist, gc}) * tot_counts;
			}
		}
	}

	ShrinkVect(comp_counts);
}

namespace reseq{
	TEST_F(ProbabilityEstimatesTest, DataStorageSetUp){
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		uint64_t n(0);
		for(auto i=0; i<3; ++i){
			if(1 != i){
				for(auto j=1; j<3; ++j){
					for(auto k=0; k<3; ++k){
						if(2 != k){
							counts[i][j][k][0][0] = ++n*100; // Multiply everything by 100 so we are above minimum counts and bins are not combined
						}
						else{
							counts[i][j][k][0][0] = 0;
						}
					}
				}
			}
		}

		array< vector<uint64_t>, 3 > dim_indices;
		array< vector<uint16_t>, 3 > initial_dim_indices;
		ProbabilityEstimatesSubClasses::DataStorage<3> data;
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
					counts.at(i).at(3).at(k).at(0).at(0) = counts.at(i).at(2).at(k).at(0).at(0)+24;
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
		uint16_t base(1);
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,7635921665);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IPFStepWiseQual(base, margins, margin_quality_position);
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQuality){
		uint16_t base(0);
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,103741084);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);
		ASSERT_TRUE( test_.Save(save_test_file_.c_str()) );

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultQual(test_, base, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_test_file_.c_str()) );
		iterations = GetIterationsQual(test2, base);
		precision = GetPrecisionQual(test2, base);
		GetIPFResultQual(test2, base, comp_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "loading" );

		EXPECT_EQ( 0, remove(save_test_file_.c_str()) ) << "Error deleting file: " << save_test_file_ << '\n';

		std::vector<double> prob;
		double prob_sum;
		test_.Quality(0, 0, base).Draw(prob,prob_sum,{200, 200, 200, 200},0.5); // Check for a crash if random values are queried (tile_id=0 must exist, otherwise crash is expected)
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithOffset){
		uint16_t base(2);
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,213576235);
		SetOffset(counts, 37, 0, 37, 0, 0);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultQual(test_, base, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithZeroMarginValues){
		uint16_t base(3);
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,3290428530);
		SetOffset(counts, 37, 0, 37, 0, 0);
		counts[38][2].Clear(); // Set quality 38, error rate 2 margin to zero

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultQual(test_, base, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforQualityWithZeroOverAllMargins){
		uint16_t base(0);
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,3,5,3,5,1,9527958903);
		SetOffset(counts, 37, 0, 37, 0, 0);
		counts[38].Clear(); // Set quality 38 to zero for all margins

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		vector< pair<bool, bool> > margin_def;
		SetMarginDefQual(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingQual(base, margins, margin_quality_position);

		auto iterations = GetIterationsQual(test_, base);
		auto precision = GetPrecisionQual(test_, base);

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultQual(test_, base, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "estimation" );
	}

	TEST_F(ProbabilityEstimatesTest, IPFforBaseCalls){
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,5,3,5,3,5,768293592);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
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
		ASSERT_TRUE( test_.Save(save_test_file_.c_str()) );

		auto iterations = GetIterationsBaseCall(test_);
		auto precision = GetPrecisionBaseCall(test_);

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultBaseCall(test_, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_test_file_.c_str()) );
		iterations = GetIterationsBaseCall(test2);
		precision = GetPrecisionBaseCall(test2);
		GetIPFResultBaseCall(test2, comp_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( iterations, precision, margins, margin_quality_position, margin_def, comp_counts, "loading" );

		EXPECT_EQ( 0, remove(save_test_file_.c_str()) ) << "Error deleting file: " << save_test_file_ << '\n';
	}

	TEST_F(ProbabilityEstimatesTest, IPFforDominantError){
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,5,3,3,1,1,793256239);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		// {use margin_quality_position, flip dimensions compared to dimension order}
		vector< pair<bool, bool> > margin_def;
		margin_def.push_back({false, true}); // dist, domErr
		margin_def.push_back({false, true}); // gc, domErr
		margin_def.push_back({false, false}); // dist, gc
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingDomError(margins);
		ASSERT_TRUE( test_.Save(save_test_file_.c_str()) );

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultDomError(test_, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( GetIterationsDomError(test_), GetPrecisionDomError(test_), margins, margin_quality_position, margin_def, comp_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_test_file_.c_str()) );
		GetIPFResultDomError(test2, comp_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( GetIterationsDomError(test2), GetPrecisionDomError(test2), margins, margin_quality_position, margin_def, comp_counts, "loading" );

		EXPECT_EQ( 0, remove(save_test_file_.c_str()) ) << "Error deleting file: " << save_test_file_ << '\n';
	}

	TEST_F(ProbabilityEstimatesTest, IPFforErrorRate){
		Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > counts;
		RandomFillingTestCounts(counts,5,3,5,1,1,3465435782);

		vector<Vect<Vect<uint64_t>>> margins;
		Vect<SeqQualityStats<uint64_t>> margin_quality_position;
		// {use margin_quality_position, flip dimensions compared to dimension order}
		vector< pair<bool, bool> > margin_def;
		SetMarginDefErrorRate(margin_def);
		CalculateMargins(margin_def, margins, margin_quality_position, counts);

		IterativeProportionalFittingErrorRate(margins);
		ASSERT_TRUE( test_.Save(save_test_file_.c_str()) );

		Vect< Vect< Vect< Vect< Vect<double> > > > > comp_counts;
		GetIPFResultErrorRate(test_, comp_counts, margins);

		test_.PrepareResult();
		CheckIPFResult( GetIterationsErrorRate(test_), GetPrecisionErrorRate(test_), margins, margin_quality_position, margin_def, comp_counts, "estimation" );

		ProbabilityEstimates test2;
		ASSERT_TRUE( test2.Load(save_test_file_.c_str()) );
		GetIPFResultErrorRate(test2, comp_counts, margins);
		test2.PrepareResult();
		CheckIPFResult( GetIterationsErrorRate(test2), GetPrecisionErrorRate(test2), margins, margin_quality_position, margin_def, comp_counts, "loading" );

		EXPECT_EQ( 0, remove(save_test_file_.c_str()) ) << "Error deleting file: " << save_test_file_ << '\n';
	}
}
