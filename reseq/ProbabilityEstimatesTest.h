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
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	class ProbabilityEstimatesTest : public BasicTestClass{
	public: 
		static void Register();

	private:
		const double kEpsilon = 1.0e-13;

		const uintTempSeq kTemplateSegment = 0;
		const uintNumFits kMaxIterations = 150;
		const double kPrecisionAim = 1.05;

		inline double CalculateAbsolutPrecision(uintMatrixCount correct_margin, double marginal_sum, double precision);
		void SetUpDataQual(DataStats &stats, uintBaseCall base, const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position);
		void SetUpDataBaseCall(DataStats &stats, const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position);
		void SetUpDataDomError(DataStats &stats, const std::vector<Vect<Vect<uintMatrixCount>>> &margins);
		void SetUpDataErrorRate(DataStats &stats, const std::vector<Vect<Vect<uintMatrixCount>>> &margins);

	protected:
		const std::string kSaveTestFile = std::string(PROJECT_SOURCE_DIR)+"/test/saveTest.reseq.ipf";

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

		void RandomFillingTestCounts( Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts, const uintMatrixIndex dim1, const uintMatrixIndex dim2, const uintMatrixIndex dim3, const uintMatrixIndex dim4, const uintMatrixIndex dim5, const uintSeed seed);
		void SetOffset( Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts, const uintMatrixIndex dim1, const uintMatrixIndex dim2, const uintMatrixIndex dim3, const uintMatrixIndex dim4, const uintMatrixIndex dim5);
		template<typename T> void CalculateMargins( const std::vector< std::pair<bool, bool> > &margin_def, std::vector<Vect<Vect<T>>> &margins, Vect<SeqQualityStats<T>> &margin_quality_position, const Vect< Vect< Vect< Vect< Vect<T> > > > > &counts ){
			margin_quality_position.Clear();
			margins.clear();
			margins.resize(margin_def.size());
			std::array<uintMatrixIndex, 5> pos;

			for( pos[0]=counts.to(); counts.from() < pos[0]--; ){ // counting backwards to get the least amount of storage allocations
				for( pos[1]=counts[pos[0]].to(); counts[pos[0]].from() < pos[1]--; ){
					for( pos[2]=counts[pos[0]][pos[1]].to(); counts[pos[0]][pos[1]].from() < pos[2]--; ){
						for( pos[3]=counts[pos[0]][pos[1]][pos[2]].to(); counts[pos[0]][pos[1]][pos[2]].from() < pos[3]--; ){
							for( pos[4]=counts[pos[0]][pos[1]][pos[2]][pos[3]].to(); counts[pos[0]][pos[1]][pos[2]][pos[3]].from() < pos[4]--; ){
								// Fill all defined margins
								uintMarginId dim1 = 0;
								uintMarginId dim2 = 1;
								for(uintMarginId i=0; i < margin_def.size(); ++i){
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
		void CheckIPFResult( uintNumFits iterations, double precision, const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position, const std::vector< std::pair<bool, bool> > &margin_def, const Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::string context );

		void SetMarginDefQual(std::vector< std::pair<bool, bool> > &margin_def);
		void SetUpDataStorageErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, std::array< std::vector<uintMatrixIndex>, 3 > &dim_indices, std::array< std::vector<uintMatrixIndex>, 3 > &initial_dim_indices, const Vect< Vect< Vect< Vect< Vect<uintMatrixCount> > > > > &counts);
		void TestDataStorageSetUpErrorRate(ProbabilityEstimatesSubClasses::DataStorage<3> &data, std::array< std::vector<uintMatrixIndex>, 3 > &dim_indices, std::array< std::vector<uintMatrixIndex>, 3 > &initial_dim_indices, bool middle_bin_inserted=false);
		void IterativeProportionalFittingQual(uintBaseCall base, const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position);
		void GetIPFResultQual( const ProbabilityEstimates &estimate, uintBaseCall base, Vect< Vect< Vect< Vect< Vect<double> > > > > &estimated_counts, const std::vector<Vect<Vect<uintMatrixCount>>> &margins);
		template<uintMarginId U1, uintMarginId U2, uintMarginId L1, uintMarginId L2, uintMarginId L3> void IPFStepQual( uintBaseCall base, const ProbabilityEstimatesSubClasses::DataStorage<5> &data ){
			test_.quality_[kTemplateSegment][0][base].IPFStep<U1,U2,L1,L2,L3>(data);

			double marginal_sum, marginal;
			for( uintMatrixIndex dim1 = 0; dim1 < data.size<U1>(); ++dim1 ){
				for( uintMatrixIndex dim2 = 0; dim2 < data.size<U2>(); dim2++ ){
					marginal_sum = ProbabilityEstimatesSubClasses::SumLogTerms<U1,U2,L1,L2,L3>(data, test_.quality_[kTemplateSegment][0][base].estimates_, dim1, dim2);
					marginal = data.get<U1,U2>(dim1, dim2);

					EXPECT_NEAR( marginal, marginal_sum, kEpsilon ) << "margin<" << U1 << ',' << U2 << "> not correct directly after its fitting step\n";
				}
			}
		}
		void IPFStepWiseQual(uintBaseCall base, const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position);
		uintNumFits GetIterationsQual( const ProbabilityEstimates &estimate, uintBaseCall base ){
			return estimate.quality_[kTemplateSegment][0][base].steps_;
		}
		double GetPrecisionQual( const ProbabilityEstimates &estimate, uintBaseCall base ){
			return estimate.quality_[kTemplateSegment][0][base].precision_;
		}

		void IterativeProportionalFittingBaseCall( const std::vector<Vect<Vect<uintMatrixCount>>> &margins, const Vect<SeqQualityStats<uintMatrixCount>> &margin_quality_position );
		void GetIPFResultBaseCall( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uintMatrixCount>>> &margins );
		uintNumFits GetIterationsBaseCall( const ProbabilityEstimates &estimate ){
			return estimate.base_call_[kTemplateSegment][0][0][0].steps_;
		}
		double GetPrecisionBaseCall( const ProbabilityEstimates &estimate ){
			return estimate.base_call_[kTemplateSegment][0][0][0].precision_;
		}

		void IterativeProportionalFittingDomError( const std::vector<Vect<Vect<uintMatrixCount>>> &margins );
		void GetIPFResultDomError( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uintMatrixCount>>> &margins );
		uintNumFits GetIterationsDomError( const ProbabilityEstimates &estimate ){
			return estimate.dom_error_[0][0][0].steps_;
		}
		double GetPrecisionDomError( const ProbabilityEstimates &estimate ){
			return estimate.dom_error_[0][0][0].precision_;
		}

		void SetMarginDefErrorRate(std::vector< std::pair<bool, bool> > &margin_def);
		void IterativeProportionalFittingErrorRate( const std::vector<Vect<Vect<uintMatrixCount>>> &margins );
		void GetIPFResultErrorRate( const ProbabilityEstimates &estimate, Vect< Vect< Vect< Vect< Vect<double> > > > > &comp_counts, std::vector<Vect<Vect<uintMatrixCount>>> &margins );
		uintNumFits GetIterationsErrorRate( const ProbabilityEstimates &estimate ){
			return estimate.error_rate_[0][0].steps_;
		}
		double GetPrecisionErrorRate( const ProbabilityEstimates &estimate ){
			return estimate.error_rate_[0][0].precision_;
		}
	};
}

#endif // PROBABILITYESTIMATESTEST_H
