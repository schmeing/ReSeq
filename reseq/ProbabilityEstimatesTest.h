#ifndef PROBABILITYESTIMATESTEST_H
#define PROBABILITYESTIMATESTEST_H
#include "ProbabilityEstimates.h"

#include <array>
#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"
#include "CMakeConfig.h"
#include "Vect.hpp"

namespace reseq{
	class ProbabilityEstimatesTest : public BasicTestClass{
	public: 
		static void Register();

	private:
		const double epsilon_ = 1.0e-13;

		const uint16_t template_segment_ = 0;
		const uint16_t max_iterations_ = 150;
		const double precision_aim_ = 1.05;

		inline double CalculateAbsolutPrecision(uint64_t correct_margin, double marginal_sum, double precision);
		void SetUpDataQual(DataStats &stats, uint16_t base, const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position);
		void SetUpDataBaseCall(DataStats &stats, const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position);
		void SetUpDataDomError(DataStats &stats, const std::vector<Vect<Vect<uint64_t>>> &margins);
		void SetUpDataErrorRate(DataStats &stats, const std::vector<Vect<Vect<uint64_t>>> &margins);

	protected:
		const std::string save_test_file_ = std::string(PROJECT_SOURCE_DIR)+"/test/saveTest.reseq.ipf";

		// Class instance to test
		ProbabilityEstimates test_;

		// Google test
		virtual void SetUp();

		// Own functions
		void TestLogArrayCalcNormalize();
		void TestDataStorageReduceAndExpand();
		void TestDataStorageReduceAndExpand2();
		void TestLogArrayCalcExpand();
		void TestCombineDimIndices();

		void ResetProbabilityEstimates();

		void RandomFillingTestCounts( Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts, const uint16_t dim1, const uint16_t dim2, const uint16_t dim3, const uint16_t dim4, const uint16_t dim5, const uint64_t seed);
		void SetOffset( Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts, const uint16_t dim1, const uint16_t dim2, const uint16_t dim3, const uint16_t dim4, const uint16_t dim5);
		template<typename T> void CalculateMargins( const std::vector< std::pair<bool, bool> > &margin_def, std::vector<Vect<Vect<T>>> &margins, Vect<SeqQualityStats<T>> &margin_quality_position, const Vect< Vect< Vect< Vect< Vect<T> > > > > &counts ){
			margin_quality_position.Clear();
			margins.clear();
			margins.resize(margin_def.size());
			std::array<uint64_t, 5> pos;

			for( pos[0]=counts.to(); counts.from() < pos[0]--; ){ // counting backwards to get the least amount of storage allocations
				for( pos[1]=counts[pos[0]].to(); counts[pos[0]].from() < pos[1]--; ){
					for( pos[2]=counts[pos[0]][pos[1]].to(); counts[pos[0]][pos[1]].from() < pos[2]--; ){
						for( pos[3]=counts[pos[0]][pos[1]][pos[2]].to(); counts[pos[0]][pos[1]][pos[2]].from() < pos[3]--; ){
							for( pos[4]=counts[pos[0]][pos[1]][pos[2]][pos[3]].to(); counts[pos[0]][pos[1]][pos[2]][pos[3]].from() < pos[4]--; ){
								// Fill all defined margins
								uint16_t dim1 = 0;
								uint16_t dim2 = 1;
								for(uint16_t i=0; i < margin_def.size(); ++i){
									if(margin_def[i].second){
										// Switch dimension1 and 2
										auto tmp = dim1;
										dim1 = dim2;
										dim2 = tmp;
									}

									if(margin_def[i].first){
										margin_quality_position[pos[dim1]][pos[dim2]] += counts[pos[0]][pos[1]][pos[2]][pos[3]][pos[4]];
									}
									else{
										margins[i][pos[dim1]][pos[dim2]] += counts[pos[0]][pos[1]][pos[2]][pos[3]][pos[4]];

									}

									if(margin_def[i].second){
										// Switch dimension1 and 2 back
										auto tmp = dim1;
										dim1 = dim2;
										dim2 = tmp;
									}

									if( ++dim1 == dim2 ){
										dim1 = 0;
										++dim2;
									}
								}
							}
						}
					}
				}
			}

			for(auto &vect : margins){
				ShrinkVect(vect);
			}
			ShrinkVect(margin_quality_position);
		}
		void CheckIPFResult( uint32_t iterations, double precision, const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position, const std::vector< std::pair<bool, bool> > &margin_def, const Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::string context );

		void SetMarginDefQual(std::vector< std::pair<bool, bool> > &margin_def);
		void SetUpDataStorageErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, std::array< std::vector<uint64_t>, 3 > &dim_indices, std::array< std::vector<uint16_t>, 3 > &initial_dim_indices, const Vect< Vect< Vect< Vect< Vect<uint64_t> > > > > &counts);
		void TestDataStorageSetUpErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, std::array< std::vector<uint64_t>, 3 > &dim_indices, std::array< std::vector<uint16_t>, 3 > &initial_dim_indices, bool middle_bin_inserted=false);
		void IterativeProportionalFittingQual(uint16_t base, const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position);
		void GetIPFResultQual( const ProbabilityEstimates &estimate, uint16_t base, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uint64_t>>> &margins);
		template<uint16_t U1, uint16_t U2, uint16_t L1, uint16_t L2, uint16_t L3> void IPFStepQual( uint16_t base, const ProbabilityEstimatesSubClasses::DataStorage<5> &data ){
			test_.quality_[template_segment_][0][base].IPFStep<U1,U2,L1,L2,L3>(data);

			double marginal_sum, marginal;
			for( uint32_t dim1 = 0; dim1 < data.size<U1>(); ++dim1 ){
				for( uint32_t dim2 = 0; dim2 < data.size<U2>(); dim2++ ){
					marginal_sum = ProbabilityEstimatesSubClasses::SumLogTerms<U1,U2,L1,L2,L3>(data, test_.quality_[template_segment_][0][base].estimates_, dim1, dim2);
					marginal = data.get<U1,U2>(dim1, dim2);

					EXPECT_NEAR( marginal, marginal_sum, epsilon_ ) << "margin<" << U1 << ',' << U2 << "> not correct directly after its fitting step\n";
				}
			}
		}
		void IPFStepWiseQual(uint16_t base, const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position);
		uint32_t GetIterationsQual( const ProbabilityEstimates &estimate, uint16_t base ){
			return estimate.quality_[template_segment_][0][base].steps_;
		}
		double GetPrecisionQual( const ProbabilityEstimates &estimate, uint16_t base ){
			return estimate.quality_[template_segment_][0][base].precision_;
		}

		void IterativeProportionalFittingBaseCall( const std::vector<Vect<Vect<uint64_t>>> &margins, const Vect<SeqQualityStats<uint64_t>> &margin_quality_position );
		void GetIPFResultBaseCall( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uint64_t>>> &margins );
		uint32_t GetIterationsBaseCall( const ProbabilityEstimates &estimate ){
			return estimate.base_call_[template_segment_][0][0][0].steps_;
		}
		double GetPrecisionBaseCall( const ProbabilityEstimates &estimate ){
			return estimate.base_call_[template_segment_][0][0][0].precision_;
		}

		void IterativeProportionalFittingDomError( const std::vector<Vect<Vect<uint64_t>>> &margins );
		void GetIPFResultDomError( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uint64_t>>> &margins );
		uint32_t GetIterationsDomError( const ProbabilityEstimates &estimate ){
			return estimate.dom_error_[0][0][0].steps_;
		}
		double GetPrecisionDomError( const ProbabilityEstimates &estimate ){
			return estimate.dom_error_[0][0][0].precision_;
		}

		void SetMarginDefErrorRate(std::vector< std::pair<bool, bool> > &margin_def);
		void IterativeProportionalFittingErrorRate( const std::vector<Vect<Vect<uint64_t>>> &margins );
		void GetIPFResultErrorRate( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uint64_t>>> &margins );
		uint32_t GetIterationsErrorRate( const ProbabilityEstimates &estimate ){
			return estimate.error_rate_[0][0].steps_;
		}
		double GetPrecisionErrorRate( const ProbabilityEstimates &estimate ){
			return estimate.error_rate_[0][0].precision_;
		}
	};
}

#endif // PROBABILITYESTIMATESTEST_H
