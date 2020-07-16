#ifndef PROBABILITYESTIMATES_H
#define PROBABILITYESTIMATES_H

#include <algorithm>
#include <array>
#include <atomic>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "SeqQualityStats.hpp"
#include "utilities.hpp"
#include "Vect.hpp"
#include "DataStats.h"

namespace reseq{
	class ProbabilityEstimatesTest;

	namespace ProbabilityEstimatesSubClasses{
		template<uintMarginId A, uintMarginId B, uintMarginId N> inline uintMarginId MapDim2To1(){
			return A+B*(2*N-B-3)/2-1;
		}

		inline void MaxPrecision(double &maximum_precision, double desired_value, double given_value){
			if( desired_value ){
				double precision( std::abs(given_value-desired_value)/desired_value );
				if( precision > maximum_precision){
					maximum_precision = precision;
				}
			}
		}

		template<size_t N> inline void CombineDimIndices(std::array< std::vector<uintMatrixIndex>, N > &dim_indeces1, const std::array< std::vector<uintMatrixIndex>, N > &dim_indeces2){
			for( uintMarginId cur_dim = 0; cur_dim < N; ++cur_dim ){
				for(uintMatrixIndex bin = 0; bin < dim_indeces1.at(cur_dim).size(); ++bin){
					dim_indeces1.at(cur_dim).at(bin) = dim_indeces2.at(cur_dim).at( dim_indeces1.at(cur_dim).at(bin) );
				}
			}
		}

		template<uintMarginId N> class DataStorage{
		public:
			// Definitions
			const uintMatrixCount kMinCountsPer1dBin = 100;
			const uintMatrixIndex kMaxBinsPerDimension = 200;

		private:
			static const uintMarginId kNumMargins = N*(N-1)/2;

			// Storage
			std::array<std::vector<double>, kNumMargins> data_;
			std::array<uintMatrixIndex, N> dim_size_;

			void CombineEndBinsToMinCount(std::array< std::vector<uintMatrixIndex>, N > &dim_indices, std::array<Vect<uintMatrixCount>, N> &dim_control);
			void CombineBins(std::array<std::vector<uintMatrixIndex>, N> &initial_dim_indices_reduced, std::array<Vect<uintMatrixCount>, N> &dim_control);

			void FillDiffMatrix(std::vector<double> &diff_matrix, uintMarginId dimension) const;
			double MaxDiff(const std::vector<double> &diff_matrix, uintMatrixIndex combined_index, std::vector<uintMatrixIndex> &same_combined_index, const std::vector<uintMatrixIndex> &dim_indices_reduced, const std::vector<uintMatrixIndex> &dim_indices_count) const;

			// Google test
			friend class reseq::ProbabilityEstimatesTest;
		public:
			DataStorage(){
				dim_size_.fill(0);
			}

			bool SetUp(const std::array< std::pair< const Vect<Vect<uintMatrixCount>> *, bool>, N*(N-1)/2 > &margins, const Vect<SeqQualityStats<uintMatrixCount>> *alternative_margin, std::array< std::vector<uintMatrixIndex>, N > &dim_indices, std::array<std::vector<uintMatrixIndex>, N> &initial_dim_indices_reduced, std::mutex &print_mutex);

			template<uintMarginId A, uintMarginId B> inline double get(uintMatrixIndex index_a, uintMatrixIndex index_b) const{
				if( A > B ){
					return data_.at(MapDim2To1<A,B,N>()).at(index_a*dim_size_.at(B)+index_b);
				}
				else{
					return data_.at(MapDim2To1<B,A,N>()).at(index_b*dim_size_.at(A)+index_a);
				}
			}

			template<uintMarginId A> inline uintMatrixIndex size() const{
				return dim_size_.at(A);
			}

			void Reduce(DataStorage<N> &reduced_data, const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_reduced, const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count) const;
			void Reduce(DataStorage<N> &reduced_data, std::array<std::vector<uintMatrixIndex>, N> &dim_indices_reduced, std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count, uintMatrixIndex max_size_dim) const;
			bool Expand(DataStorage<N> &reduced_data, std::array<std::vector<uintMatrixIndex>, N> &dim_indices_reduced, std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count, std::array<std::vector<uintMatrixIndex>, N> &expansion_indices, std::array<std::vector<uintMatrixIndex>, N> &expansion_count, uintMatrixIndex max_duplications) const;
		};

		template<uintMarginId N> class LogArrayCalc{
		private:
			template<uintMarginId M> friend class LogArrayResult;

			static const uintMarginId kNumMargins = N*(N-1)/2;

			std::array<std::vector<double>, kNumMargins> dim2_;
			std::array<uintMatrixIndex, N> dim_size_;

			inline void ResizeMarginalSums(std::array< std::vector<double>, kNumMargins> &marginal_sums) const {
				for(auto n=kNumMargins;n--;){
					marginal_sums.at(n).resize(dim2_.at(n).size(),0.0);
				}
			}

			// boost serialization
			friend class boost::serialization::access;
			template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
				ar & dim2_;
				ar & dim_size_;
			}

			// Google test
			friend class reseq::ProbabilityEstimatesTest;

		public:
			LogArrayCalc(){
				dim_size_.fill(0);
			}

			template<typename T> void SetUp(const std::array< std::vector<T>, N > &dim_indices){
				// Fill in dim2_
				uintMarginId dim_a(N), dim_b(N-1);
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					dim2_.at(n).clear();
					dim2_.at(n).resize(dim_indices.at(dim_a).size()*dim_indices.at(dim_b).size(), 1.0);
				}

				// Set dim_size_
				for(uintMarginId n=N; n--; ){
					dim_size_.at(n) = dim_indices.at(n).size();
				}
			}

			void Normalize(){
				double log_dim0(0.0);
				std::array<std::vector<double>, N> log_dim1;
				for( auto n=N; n--; ){
					log_dim1.at(n).resize(dim_size_.at(n), 0.0);
				}

				// Normalize dim2 so that the product of all values in a row/column are equal to 1
				// Shift the normalization to dim1 and from there to dim0 and then redistribute equally in the other direction
				// First change values so that the values exp(log_dim1[0]+log_dim1[1])*dim2[01] don't change
				uintMarginId dim_a(N), dim_b(N-1);
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					// Normalize in direction dim_b
					for(uintMatrixIndex i = dim_size_.at(dim_a); i--; ){
						// Find the value to normalize
						double prod(0.0);
						uintMatrixCount count(0);
						for(uintMatrixIndex j = dim_size_.at(dim_b); j--; ){
							if(0.0 != dim2_.at(n).at(i*dim_size_.at(dim_b)+j)){
								prod += std::log(dim2_.at(n).at(i*dim_size_.at(dim_b)+j)); // Use the logarithm here so that we don't run into a problem with overflow or underflow
								++count;
							}
						}
						if(count){
							prod = prod/count;
						}

						// Normalize so the product = 1.0
						log_dim1.at(dim_a).at(i) += prod;
						prod = 1/std::exp(prod);
						for(uintMatrixIndex j = dim_size_.at(dim_b); j--; ){
							dim2_.at(n).at(i*dim_size_.at(dim_b)+j) *= prod;
						}
					}

					// Normalize in direction dim_a
					for(uintMatrixIndex j = dim_size_.at(dim_b); j--; ){
						// Find the value to normalize
						double prod(0.0);
						uintMatrixCount count(0);
						for(uintMatrixIndex i = dim_size_.at(dim_a); i--; ){
							if(0.0 != dim2_.at(n).at(i*dim_size_.at(dim_b)+j)){
								prod += std::log(dim2_.at(n).at(i*dim_size_.at(dim_b)+j));
								++count;
							}
						}
						if(count){
							prod = prod/count;
						}

						// Normalize so the product = 1.0
						log_dim1.at(dim_b).at(j) += prod;
						prod = 1/std::exp(prod);
						for(uintMatrixIndex i = dim_size_.at(dim_a); i--; ){
							dim2_.at(n).at(i*dim_size_.at(dim_b)+j) *= prod;
						}
					}
				}

				// Normalize dim1 so that the values log_dim0+log_dim1 don't change
				for(uintMarginId n=N; n--; ){
					// Find the value to normalize
					double sum(0.0);
					for(uintMatrixIndex i = dim_size_.at(n); i--; ){
						sum += log_dim1.at(n).at(i);
					}
					sum = sum/dim_size_.at(n);

					// Normalize so the product = 1.0
					log_dim0 += sum;
					for(uintMatrixIndex i = dim_size_.at(n); i--; ){
						log_dim1.at(n).at(i) -= sum;
					}
				}

				// Now redistribute log_dim0 and log_dim1 equally to dim2_ again, so that exp(log_dim0+log_dim1)*dim2_ doesn't change, but log_dim0=log_dim1=0 afterwards
				// First store the factor to apply in log_dim1
				log_dim0 /= N; // dim0 is separated onto N dim1 variables
				for(uintMarginId n=N; n--; ){
					for(uintMatrixIndex i = dim_size_.at(n); i--; ){
						// dim1 variables are each separated into N-1 dim2 variables
						log_dim1.at(n).at(i) = std::exp( (log_dim1.at(n).at(i)+log_dim0)/(N-1) );
					}
				}

				// Now apply those factors
				dim_a = N;
				dim_b = N-1;
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					// Normalize in direction dim_b
					for(uintMatrixIndex i = dim_size_.at(dim_a); i--; ){
						for(uintMatrixIndex j = dim_size_.at(dim_b); j--; ){
							dim2_.at(n).at(i*dim_size_.at(dim_b)+j) *= log_dim1.at(dim_a).at(i) * log_dim1.at(dim_b).at(j);
						}
					}
				}
			}

			void Expand(const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_reduced, const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count, const std::array<std::vector<uintMatrixIndex>, N> *dim_indices_new_count=NULL){
				// Each value dim2[0] * dim2[1] * ... is divided by the number of new values it will become and then filled into those values
				// This happens by equally splitting the division over all N-1 margins containing the variable

				// Prepare the margin/((count/new_count)^(1/(N-1))) already in multipliers, so that margin*multipliers
				std::array<std::vector<double>, N> multipliers;
				for( auto n=N; n--; ){
					multipliers.at(n).resize(dim_indices_reduced.at(n).size());

					for(uintMatrixIndex i = dim_indices_reduced.at(n).size(); i--; ){
						multipliers.at(n).at(i) = std::pow( (dim_indices_new_count?static_cast<double>((*dim_indices_new_count).at(n).at(i)):1.0)/dim_indices_count.at(n).at(dim_indices_reduced.at(n).at(i)), 1.0/(N-1));
					}
				}

				// Now copy the value to all new values and apply multipliers from the two relevant dimensions
				std::vector<double> old_values;
				uintMarginId dim_a(N), dim_b(N-1);
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					old_values = std::move(dim2_.at(n));
					dim2_.at(n).resize(dim_indices_reduced.at(dim_a).size()*dim_indices_reduced.at(dim_b).size());
					for(uintMatrixIndex i = dim_indices_reduced.at(dim_a).size(); i--; ){
						for(uintMatrixIndex j = dim_indices_reduced.at(dim_b).size(); j--; ){
							dim2_.at(n).at(i*dim_indices_reduced.at(dim_b).size()+j) = old_values.at( dim_indices_reduced.at(dim_a).at(i)*dim_indices_count.at(dim_b).size()+dim_indices_reduced.at(dim_b).at(j) ) * multipliers.at(dim_a).at(i) * multipliers.at(dim_b).at(j);
						}
					}
				}

				// Update size array
				for( auto n=N; n--; ){
					dim_size_.at(n) = dim_indices_reduced.at(n).size();
				}
			}

			bool CheckConsistency(const std::string &descriptor, std::mutex &print_mutex){
				uintMarginId dim_a(N), dim_b(N-1);
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					for(auto i = dim2_.at(n).size(); i--; ){
						if( std::isnan(dim2_.at(n).at(i)) ){
							print_mutex.lock();
							printErr << "dim2_[" << n << "][" << i/dim_size_.at(dim_b) << "][" << i%dim_size_.at(dim_b) << "] is NaN for " << descriptor << std::endl;
							print_mutex.unlock();
							return false;
						}
					}
				}

				return true;
			}

			void Clear(){
				for(auto n=kNumMargins; n--;){
					dim2_.at(n).clear();
					dim2_.at(n).shrink_to_fit();
				}
				dim_size_.fill(0);
			}

			template<uintMarginId A, uintMarginId B> inline double &dim2( decltype(dim2_.at(0).size()) index1, decltype(dim2_.at(0).size()) index2 ){
				if( A > B ){
					return dim2_.at(MapDim2To1<A,B,N>()).at(index1*dim_size_.at(B)+index2);
				}
				else{
					return dim2_.at(MapDim2To1<B,A,N>()).at(index2*dim_size_.at(A)+index1);
				}
			}
			template<uintMarginId A, uintMarginId B> inline const double &dim2( decltype(dim2_.at(0).size()) index1, decltype(dim2_.at(0).size()) index2 ) const{
				if( A > B ){
					return dim2_.at(MapDim2To1<A,B,N>()).at(index1*dim_size_.at(B)+index2);
				}
				else{
					return dim2_.at(MapDim2To1<B,A,N>()).at(index2*dim_size_.at(A)+index1);
				}
			}
			inline double GetMatrixElement(const std::array<uintMatrixIndex, N> &index) const{
				double prod(1.0);

				uintMarginId dim_a(N), dim_b(N-1);
				for(uintMarginId n=kNumMargins; n--; ){
					if(--dim_a == dim_b){
						--dim_b;
						dim_a = N-1;
					}

					prod *= dim2_.at(n).at(index.at(dim_a)*dim_size_.at(dim_b)+index.at(dim_b));
				}

				return prod;
			}
		};

		template<uintMarginId N> class LogArrayResult{
		private:
			static const uintMarginId kNumMargins = (N-1); // Here we only take the margins with the first parameter in them, as the rest is irrelevant

			std::array<std::vector<double>, kNumMargins> dim2_;
			std::array<std::pair<uintMatrixIndex, uintMatrixIndex>, N-1> limits_;
			std::vector<uintMatrixIndex> par0_indeces_;

			inline double Likelihood(uintMatrixIndex ind0, const std::array<uintMatrixIndex, N-1> &index) const{
				double prod(dim2_.at(0).at(index.at(0)*par0_indeces_.size()+ind0));
				for(uintMarginId n=1; n<kNumMargins; ++n ){
					prod *= dim2_.at(n).at(index.at(n)*par0_indeces_.size()+ind0);
				}

				return prod;
			}

			inline void AdjustIndeces(std::array<uintMatrixIndex, N-1> &index) const{
				for(auto n=limits_.size(); n--; ){
					if(index.at(n) < limits_.at(n).first){
						index.at(n) = 0; // Set it equal to the first element
					}
					else if(index.at(n) >= limits_.at(n).second){
						index.at(n) = limits_.at(n).second-limits_.at(n).first-1; // Set it equal to the last element
					}
					else{
						index.at(n) -= limits_.at(n).first; // Only subtract offset
					}
				}
			}

		public:
			LogArrayResult(){
				limits_.fill({0,0});
			}
			void GetResults(const LogArrayCalc<N> &calc, const std::array< std::vector<uintMatrixIndex>, N > &dim_indices ){
				if(dim_indices.at(0).size()){
					// Order first parameter by increasing mean, so that when called from the back the values with the highest probability come first
					std::vector<std::pair<double, uintMatrixIndex>> dim_order;
					dim_order.resize(dim_indices.at(0).size());
					for(auto j=dim_order.size(); j--; ){
						dim_order.at(j) = {0.0, j};
					}

					uintMarginId dim_a;
					double sum;
					for(uintMarginId n=kNumMargins; n--; ){
						dim_a = n+1;

						for(auto j = dim_order.size(); j--; ){
							sum = 0.0;
							for(auto i = dim_indices.at(dim_a).size(); i--; ){
								sum += calc.dim2_.at(n).at(i*dim_order.size()+j);
							}
							dim_order.at(j).first += sum/dim_indices.at(dim_a).size();
						}
					}
					std::sort(dim_order.begin(), dim_order.end());

					par0_indeces_.resize(dim_order.size());
					for(auto k=dim_order.size(); k--; ){
						par0_indeces_.at( dim_order.at(k).second ) = k;
					}

					// Get limits {from,to} for other parameter
					for(uintMarginId n=kNumMargins; n--; ){
						limits_.at(n) = {std::numeric_limits<uintMatrixIndex>::max(), 0};

						for(auto ind : dim_indices.at(n+1) ){
							if( ind < limits_.at(n).first ){
								limits_.at(n).first = ind;
							}
							if( ind > limits_.at(n).second ){
								limits_.at(n).second = ind;
							}
						}
						++limits_.at(n).second; // To is one past the last element
					}

					// Copy values to correct position
					for(uintMarginId n=kNumMargins; n--; ){
						dim_a = n+1;

						dim2_.at(n).resize( (limits_.at(n).second - limits_.at(n).first) * par0_indeces_.size() );
						for(auto i = dim_indices.at(dim_a).size(); i--; ){
							for(auto j = dim_indices.at(0).size(); j--; ){
								dim2_.at(n).at((dim_indices.at(dim_a).at(i)-limits_.at(n).first)*par0_indeces_.size() + par0_indeces_.at(j)) = calc.dim2_.at(n).at(i*dim_indices.at(0).size()+j);
							}
						}
					}

					// Set par0_indeces to contain the link from result indeces to original indeces (ignoring the calculation indeces we needed before)
					for(auto k=dim_order.size(); k--; ){
						par0_indeces_.at( k ) = dim_indices.at(0).at( dim_order.at(k).second );
					}
				}
				else{
					// No data exists so set limits to 0
					for(uintMarginId n=kNumMargins; n--; ){
						limits_.at(n) = {0, 0};
					}
				}
			}

			void ImputeMissingValues(){
				uintMatrixIndex last_filled_index;
				for(uintMarginId n=kNumMargins; n--; ){
					last_filled_index = 0;
					for(uintMatrixIndex i = 1; i < limits_.at(n).second-limits_.at(n).first; ++i ){
						bool filled(false);
						for(auto j = 0; j < par0_indeces_.size(); ++j ){
							if(0.0 != dim2_.at(n).at(i*par0_indeces_.size()+j)){
								filled=true;
								break;
							}
						}

						if(filled){
							// This line is filled, so impute all lines between this one and the last filled line
							for(auto imp_ind = last_filled_index+1; imp_ind < i; ++imp_ind){
								for(auto j = 0; j < par0_indeces_.size(); ++j ){
									dim2_.at(n).at(imp_ind*par0_indeces_.size()+j) = dim2_.at(n).at(last_filled_index*par0_indeces_.size()+j)*(imp_ind-last_filled_index)/(i-last_filled_index) + dim2_.at(n).at(i*par0_indeces_.size()+j)*(i-imp_ind)/(i-last_filled_index); // Impute linearly
								}
							}
							last_filled_index = i;
						}
					}
				}
			}

			inline uintMatrixIndex Draw(std::vector<double> &prob, double &prob_sum, std::array<uintMatrixIndex, N-1> index, double random_number) const{
				prob_sum = 0.0;
				
				if(par0_indeces_.size()){
					AdjustIndeces(index);

					prob.clear();
					prob.resize(par0_indeces_.size());

					for( uintMatrixIndex ind0 = 0; ind0 < prob.size(); ++ind0 ){
						prob.at(ind0) = Likelihood(ind0, index);
						prob_sum += prob.at(ind0);
					}

					random_number *= prob_sum; // Multiply the random number [0,1[ by the sum, which is the same as dividing all likelihoods by the sum to make them probabilities

					double sum(0.0);
					uintMatrixIndex ind0 = prob.size();
					while(sum <= random_number && --ind0){
						sum += prob.at(ind0);
					}

					return par0_indeces_.at(ind0);
				}
				else{
					return 0;
				}
			}
			
			inline uintMatrixIndex MaxValue() const{
				if(par0_indeces_.size()){
					return *std::max_element(par0_indeces_.begin(), par0_indeces_.end());
				}
				else{
					return 0;
				}
			}
			
			inline uintMatrixIndex MostLikely() const{
				if(par0_indeces_.size()){
					return par0_indeces_.at(par0_indeces_.size()-1);
				}
				else{
					return 0;
				}
			}

			inline const std::vector<uintMatrixIndex> &PossibleValues() const{
				return par0_indeces_;
			}

			inline void ModifyPar0(uintMatrixIndex par0_index, double multiplier){
				// Find which column in the matrix corresponds to the given index
				uintMatrixIndex corrected_index = 0;
				while(corrected_index < par0_indeces_.size() && par0_indeces_.at(corrected_index) != par0_index){
					++corrected_index;
				}

				if(corrected_index < par0_indeces_.size() && par0_indeces_.at(corrected_index) == par0_index){
					// par0_index exists
					for(auto matrix_index = corrected_index; matrix_index < dim2_.at(0).size(); matrix_index += par0_indeces_.size()){
						dim2_.at(0).at(matrix_index) *= multiplier;
					}
				}
			}
		};

		template<uintMarginId G1, uintMarginId G2, uintMarginId L> double SumLogTerms(
				const DataStorage<3> &data,
				const LogArrayCalc<3> &mat,
				uintMatrixIndex given_index1,
				uintMatrixIndex given_index2
				){
			double sum = 0.0;
			for( uintMatrixIndex i = 0; i < data.size<L>(); ++i ){
				sum += mat.template dim2<G1,L>(given_index1,i) * mat.template dim2<G2,L>(given_index2,i);
			}
			return sum * mat.template dim2<G1,G2>(given_index1,given_index2);
		}

		template<uintMarginId G1, uintMarginId G2, uintMarginId L1, uintMarginId L2> double SumLogTerms(
				const DataStorage<4> &data,
				const LogArrayCalc<4> &mat,
				uintMatrixIndex given_index1,
				uintMatrixIndex given_index2
				){
			double tmp_sum;
			double sum = 0.0;
			for( uintMatrixIndex i = 0; i < data.size<L2>(); ++i ){
				tmp_sum = 0.0;
				for( uintMatrixIndex j = 0; j < data.size<L1>(); ++j ){
					tmp_sum += mat.template dim2<G1,L1>(given_index1,j) * mat.template dim2<G2,L1>(given_index2,j) * mat.template dim2<L2,L1>(i,j);
				}
				sum += tmp_sum * mat.template dim2<G1,L2>(given_index1,i) * mat.template dim2<G2,L2>(given_index2,i);
			}
			return sum * mat.template dim2<G1,G2>(given_index1,given_index2);
		}

		template<uintMarginId G1, uintMarginId G2, uintMarginId L1, uintMarginId L2, uintMarginId L3> double SumLogTerms(
				const DataStorage<5> &data,
				const LogArrayCalc<5> &mat,
				uintMatrixIndex given_index1,
				uintMatrixIndex given_index2
				){
			double tmp_sum, tmp_sum2;
			double sum = 0.0;
			for( uintMatrixIndex i = 0; i < data.size<L3>(); ++i ){
				tmp_sum = 0.0;
				for( uintMatrixIndex j = 0; j < data.size<L2>(); ++j ){
					tmp_sum2 = 0.0;
					for( uintMatrixIndex k = 0; k < data.size<L1>(); ++k ){
						tmp_sum2 += mat.template dim2<G1,L1>(given_index1,k) * mat.template dim2<G2,L1>(given_index2,k) * mat.template dim2<L3,L1>(i,k) * mat.template dim2<L2,L1>(j,k);
					}
					tmp_sum += tmp_sum2 * mat.template dim2<G1,L2>(given_index1,j) * mat.template dim2<G2,L2>(given_index2,j) * mat.template dim2<L3,L2>(i,j);
				}
				sum += tmp_sum * mat.template dim2<G1,L3>(given_index1,i) * mat.template dim2<G2,L3>(given_index2,i);
			}
			return sum * mat.template dim2<G1,G2>(given_index1,given_index2);
		}

		template<uintMarginId N> class LogIPF{
		private:
			static const uintMarginId kNumMargins = N*(N-1)/2;

		public:
			uintNumFits steps_;
			uintNumFits needed_updates_;
			double precision_;
			std::array< double, kNumMargins > margin_precision_;
			uintMarginId last_margin_;
			std::array< uintNumFits, kNumMargins > last_update_;
			std::array< uint16_t, kNumMargins > update_dist_;
			LogArrayCalc<N> estimates_;
			std::array< std::vector<uintMatrixIndex>, N > dim_indices_; // dim_indices_[Dimension][IndexInFullDataVector] = IndexInOriginalVectFromStats : Set in DataStorage.SetUp and later used for LogArrayResult.GetResult
			std::array< std::vector<uintMatrixIndex>, N > initial_dim_indices_reduced_; // initial_dim_indices_reduced_[Dimension][IndexInFullDataVector] = IndexInInitiallyReducedVector
			std::array< std::vector<uintMatrixIndex>, N > dim_indices_reduced_; // dim_indices_reduced_[Dimension][IndexInInitiallyReducedVector] = IndexInReducedVector

		private:
			// Temporary variables
			bool normalize_;
			double precision_aim_;
			std::vector<double> update_factors_;

			// Functions
			inline void UpdateEstimate(
					double &matrix_element,
					double factor
					){
				matrix_element *= factor;
				if(std::abs(matrix_element) > 1e200 ){
					normalize_ = true; // Values are getting to big and we risk overflow or underflow
				}
			}

			template<uintMarginId U1, uintMarginId U2, uintMarginId... L> inline void IPFCheckPrecision( const DataStorage<N> &data ){
				double sum;

				for( uintMatrixIndex i = 0; i < data.template size<U1>(); ++i ){
					for( uintMatrixIndex j = 0; j < data.template size<U2>(); ++j ){
						if( 0.0 != data.template get<U1,U2>(i,j) ){
							sum = SumLogTerms<U1,U2,L...>(data, estimates_, i, j);

							MaxPrecision(precision_, data.template get<U1,U2>(i,j), sum);
						}
					}
				}
			}

			template<uintMarginId U1, uintMarginId U2, uintMarginId... L> inline void IPFStep( const DataStorage<N> &data ){
				double sum;

				for( uintMatrixIndex i = 0; i < data.template size<U1>(); ++i ){
					for( uintMatrixIndex j = 0; j < data.template size<U2>(); ++j ){
						if( 0.0 == data.template get<U1,U2>(i,j) ){
							estimates_.template dim2<U1,U2>(i,j) = 0.0;
						}
						else{
							sum = SumLogTerms<U1,U2,L...>(data, estimates_, i, j);

							UpdateEstimate(estimates_.template dim2<U1,U2>(i,j), data.template get<U1,U2>(i,j)/sum);

							MaxPrecision(precision_, data.template get<U1,U2>(i,j), sum);
						}
					}
				}
			}

			template<uintMarginId U1, uintMarginId U2, uintMarginId... L> inline void IPFPotentialStep( const DataStorage<N> &data ){
				double sum;
				update_factors_.resize(data.template size<U1>()*data.template size<U2>());

				for( uintMatrixIndex i = 0; i < data.template size<U1>(); ++i ){
					for( uintMatrixIndex j = 0; j < data.template size<U2>(); ++j ){
						if( 0.0 == data.template get<U1,U2>(i,j) ){
							update_factors_.at(i*data.template size<U2>()+j) = 0.0;
						}
						else{
							sum = SumLogTerms<U1,U2,L...>(data, estimates_, i, j);

							// Store the factors in case we want to update
							update_factors_.at(i*data.template size<U2>()+j) = data.template get<U1,U2>(i,j)/sum;

							MaxPrecision(precision_, data.template get<U1,U2>(i,j), sum);
						}
					}
				}

				// Update in case the precision aim is violated
				if(precision_ > precision_aim_){
						for( uintMatrixIndex i = 0; i < data.template size<U1>(); ++i ){
							for( uintMatrixIndex j = 0; j < data.template size<U2>(); ++j ){
								UpdateEstimate(estimates_.template dim2<U1,U2>(i,j), update_factors_.at(i*data.template size<U2>()+j));
							}
						}
				}

				update_factors_.clear();
			}

			template<uint16_t UPDATE, uintMarginId U1, uintMarginId U2, uintMarginId... L> inline void IPFUpdateOrCheck( const DataStorage<N> &data ){
				if(1 == UPDATE){
					IPFStep<U1,U2,L...>(data);
				}
				else if(2 == UPDATE){
					IPFPotentialStep<U1,U2,L...>(data);
				}
				else{
					IPFCheckPrecision<U1,U2,L...>(data);
				}

				if(normalize_){
					estimates_.Normalize();
					normalize_ = false;
				}
			}

			template<uintMarginId S, uint16_t UPDATE=1> inline void IPFStepCaller( const DataStorage<N> &data ){
				precision_ = 0;

				IPFStepCallerTemp<S,UPDATE>( data );

				if(1 == UPDATE || (2 == UPDATE && precision_ > precision_aim_)){ // IPFStep or IPFPotentialStep with update
					last_margin_ = S;

					if(precision_ > margin_precision_.at(S)){
						// If the precision has become worse since the last update, update more often
						if(1 < update_dist_.at(S)){
							--update_dist_.at(S);
						}
					}
					else{
						if(steps_ >= last_update_.at(S)+update_dist_.at(S)*kNumMargins){
							// If precision did not get worse since the last update and we called it due to precision independent updates, update less often
							++update_dist_.at(S);
						}
					}
					last_update_.at(S) = steps_;

					margin_precision_.at(S) = precision_;

					++steps_;
				}
				else if(2 == UPDATE){ // IPFPotentialStep without update
					if(precision_ > margin_precision_.at(S)){
						margin_precision_.at(S) = precision_;
					}

					++needed_updates_;
				}
				else{ // IPFCheckPrecision
					margin_precision_.at(S) = precision_;

					++needed_updates_;
				}
			}

			template<uintMarginId S, uint16_t UPDATE=1> inline void IPFStepCallerTemp( const DataStorage<N> &data );

			template<uintMarginId S> inline void IPFStepList( const DataStorage<N> &data ){
				IPFStepCaller<S>(data);
			}

			template<uintMarginId S1, uintMarginId... S> inline typename std::enable_if<sizeof...(S) != 0, void>::type IPFStepList( const DataStorage<N> &data ){
				IPFStepList<S1>(data);
				IPFStepList<S...>(data);
			}

			template<uint16_t UPDATE=1> inline void IPFStep( const DataStorage<N> &data, uintMarginId step ){
				switch(step){
				case 0:
					IPFStepCaller<0,UPDATE>(data);
					break;
				case 1:
					IPFStepCaller<1,UPDATE>(data);
					break;
				case 2:
					IPFStepCaller<2,UPDATE>(data);
					break;
				case 3:
					IPFStepCaller<3,UPDATE>(data);
					break;
				case 4:
					IPFStepCaller<4,UPDATE>(data);
					break;
				case 5:
					IPFStepCaller<5,UPDATE>(data);
					break;
				case 6:
					IPFStepCaller<6,UPDATE>(data);
					break;
				case 7:
					IPFStepCaller<7,UPDATE>(data);
					break;
				case 8:
					IPFStepCaller<8,UPDATE>(data);
					break;
				case 9:
					IPFStepCaller<9,UPDATE>(data);
					break;
				}
			}

			inline void UpdatePrecision(){
				precision_ = *std::max_element(margin_precision_.begin(), margin_precision_.end());
			}

			inline void IPFPreparation( const DataStorage<N> &data ){
				// Fixed lists that were used before:
				switch(N){
				case 3:
					IPFStepList<0,2,1>(data);
					//IPFStepList<2,0,1>(data);
					break;
				case 4:
					IPFStepList<5,1,4,0,2,3>(data);
					break;
				case 5:
					IPFStepList<2,6,7,0,8,4,3,5,1,9>(data);
					break;
				}

				UpdatePrecision();
			}

			inline void IPFIteration( const DataStorage<N> &data ){
				uintMarginId update_id(0);
				double max_precision(0);

				for( uintMarginId ind=kNumMargins; ind--; ){
					if(last_update_.at(ind)+update_dist_.at(ind)*kNumMargins <= steps_){
						// In case one margin hasn't been run for quite a while run it now
						update_id = ind;
						break;
					}

					if(ind != last_margin_ && margin_precision_.at(ind) > max_precision){
						max_precision = margin_precision_.at(ind);
						update_id = ind;
					}
				}

				IPFStep(data, update_id);

				UpdatePrecision();
			}

			inline void IPFConfirmation( const DataStorage<N> &data ){
				std::array< std::pair<double,uintMarginId>, kNumMargins > margin_precision_sorted;
				for(auto n = kNumMargins; n--;){
					if( n == last_margin_ ){
						margin_precision_sorted.at(n) = {0.0,n}; // Make sure that the margin that was just run comes up last
					}
					else{
						margin_precision_sorted.at(n) = {margin_precision_.at(n),n};
					}
				}

				std::sort(margin_precision_sorted.begin(), margin_precision_sorted.end());
				double max_precision(0.0);
				for(auto n = kNumMargins; --n;){ // Loop backwards from highest to lowest and don't do the [0](last_margin_)
					IPFStep<2>(data, margin_precision_sorted.at(n).second);

					if(precision_ > max_precision){
						max_precision = precision_;
					}
					if(precision_ > precision_aim_){
						break;
					}
				}

				precision_ = max_precision;
			}

			inline void RunAllMargins( const DataStorage<N> &data ){
				std::array< std::pair<double,uintMarginId>, kNumMargins > last_update_sorted;
				for(auto n = kNumMargins; n--;){
					last_update_sorted.at(n) = {last_update_.at(n),n};
				}

				std::sort(last_update_sorted.begin(), last_update_sorted.end());
				for(uintMarginId n = 0; n < kNumMargins; ++n){
					IPFStep(data, last_update_sorted.at(n).second);
				}

				UpdatePrecision();
			}

			inline void UpdatePrecision( const DataStorage<N> &data ){
				for(auto n = kNumMargins; n--; ){
					if(n != last_margin_){
						IPFStep<0>(data, n);
					}
				}

				margin_precision_.at(last_margin_) = 0.0;

				UpdatePrecision();
			}

			inline uintMatrixIndex ReductionInfo(std::string &reduction_descriptor, const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count){
				// Prepare
				double step_reduction_factor(1.0);
				std::stringstream reduction_descriptor_stream;
				reduction_descriptor_stream << " (";

				for( uintMarginId n=0; n < N; ++n){
					step_reduction_factor *= static_cast<double>(dim_indices_reduced_.at(n).size()) / dim_indices_count.at(n).size();
					if(n){
						reduction_descriptor_stream << ',';
					}
					reduction_descriptor_stream << dim_indices_count.at(n).size();
				}

				reduction_descriptor_stream << ')';
				reduction_descriptor = reduction_descriptor_stream.str();

				return std::round(step_reduction_factor);
			}

			bool CheckConsistency(const std::string &descriptor, std::mutex &print_mutex){
				if( std::isnan(precision_) ){
					print_mutex.lock();
					printErr << "precision_ is NaN for " << descriptor << std::endl;
					print_mutex.unlock();
					return false;
				}
				else{
					return estimates_.CheckConsistency(descriptor, print_mutex);
				}
			}

			inline void ReconstructDimIndicesCount(std::array<std::vector<uintMatrixIndex>, N> &dim_indices_count, const std::array<std::vector<uintMatrixIndex>, N> &dim_indices_reduced){
				// Reconstruct dim_indices_count from dim_indices_reduced_
				for(auto n=N; n--; ){
					dim_indices_count.at(n).resize( *std::max_element(dim_indices_reduced.at(n).begin(),dim_indices_reduced.at(n).end()) + 1, 0 );
					for( auto ind : dim_indices_reduced.at(n) ){
						++dim_indices_count.at(n).at(ind);
					}
				}
			}

			// boost serialization
			friend class boost::serialization::access;
			template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
				ar & steps_;
				ar & needed_updates_;
				ar & precision_;
				ar & margin_precision_;
				ar & last_margin_;
				ar & last_update_;
				ar & update_dist_;
				ar & estimates_;
				ar & dim_indices_;
				ar & initial_dim_indices_reduced_;
				ar & dim_indices_reduced_;
			}

			// Google test
			friend class reseq::ProbabilityEstimatesTest;
		public:
			LogIPF():
				steps_(0),
				needed_updates_(0),
				precision_(std::numeric_limits<double>::max()),
				last_margin_(0),
				normalize_(false)
			{
				margin_precision_.fill(std::numeric_limits<double>::max());
				last_update_.fill(0);
				update_dist_.fill(2);
			}

			inline void Clear(){
				steps_ = 0;
				needed_updates_ = 0;
				precision_ = std::numeric_limits<double>::max();
				margin_precision_.fill(std::numeric_limits<double>::max());
				last_margin_ = 0;
				last_update_.fill(0);
				update_dist_.fill(2);

				estimates_.Clear();
				for(auto n=N; n--; ){
					dim_indices_.at(n).clear();
					dim_indices_.at(n).shrink_to_fit();
					initial_dim_indices_reduced_.at(n).clear();
					initial_dim_indices_reduced_.at(n).shrink_to_fit();
					dim_indices_reduced_.at(n).clear();
					dim_indices_reduced_.at(n).shrink_to_fit();
				}
			}

			inline void FullExpansion(){
				bool expansion_necessary(false);

				for(auto n=N; n-- && !expansion_necessary; ){
					for( auto ind = dim_indices_reduced_.at(n).size(); ind--;){
						if(ind != dim_indices_reduced_.at(n).at(ind)){
							expansion_necessary = true;
							break;
						}
					}
				}

				for(auto n=N; n-- && !expansion_necessary; ){
					for( auto ind = initial_dim_indices_reduced_.at(n).size(); ind--;){
						if(ind != initial_dim_indices_reduced_.at(n).at(ind)){
							expansion_necessary = true;
							break;
						}
					}
				}

				if(expansion_necessary){
					std::array<std::vector<uintMatrixIndex>, N> tmp_dim_indices;
					for( uintMarginId cur_dim = 0; cur_dim < N; ++cur_dim ){
						tmp_dim_indices.at(cur_dim) = initial_dim_indices_reduced_.at(cur_dim);
					}
					CombineDimIndices(tmp_dim_indices, dim_indices_reduced_);

					std::array<std::vector<uintMatrixIndex>, N> dim_indices_count;
					ReconstructDimIndicesCount(dim_indices_count, tmp_dim_indices);
					estimates_.Expand(tmp_dim_indices,dim_indices_count);
				}
			}

			void IterativeProportionalFitting(
					std::atomic<bool> &error_during_fitting,
					std::atomic<bool> &precision_improved,
					double precision_aim,
					uintNumFits max_iterations,
					const std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, kNumMargins > &margins,
					const Vect<SeqQualityStats<uintMatrixCount>> *alternative_margin,
					const std::string &descriptor,
					std::mutex &print_mutex
					){
				DataStorage<N> data;

				if(precision_ > precision_aim){
					// We have to do something, so data has to be set up
					if( !data.SetUp(margins, alternative_margin, dim_indices_, initial_dim_indices_reduced_, print_mutex) ){
						error_during_fitting = true;
						return;
					}

					precision_aim_ = precision_aim;

					// Reserve the needed space for temporary vectors
					uintMarginId max_dim1, max_dim2;
					if(dim_indices_.at(N-1).size() > dim_indices_.at(N-2).size()){
						max_dim1 = N-1;
						max_dim2 = N-2;
					}
					else{
						max_dim1 = N-2;
						max_dim2 = N-1;
					}
					for(auto n=N-2; n--; ){
						if(dim_indices_.at(n).size() > dim_indices_.at(max_dim2).size()){
							if(dim_indices_.at(n).size() > dim_indices_.at(max_dim1).size()){
								max_dim2 = max_dim1;
								max_dim1 = n;
							}
							else{
								max_dim2 = n;
							}
						}
					}

					update_factors_.reserve(dim_indices_.at(max_dim1).size()*dim_indices_.at(max_dim2).size());
				}

				DataStorage<N> reduced_data;
				std::array<std::vector<uintMatrixIndex>, N> dim_indices_count, expansion_indices, expansion_count;
				auto old_steps(steps_);
				if(steps_){
					print_mutex.lock();
					printInfo << "Loaded " << descriptor << " at step " << steps_ << " with precision " << precision_*100 << "%" << std::endl;
					print_mutex.unlock();
					steps_ = 0;

					if(precision_ > precision_aim){
						ReconstructDimIndicesCount(dim_indices_count, dim_indices_reduced_);
						data.Reduce(reduced_data, dim_indices_reduced_, dim_indices_count);
					}
				}
				else{
					data.Reduce(reduced_data, dim_indices_reduced_, dim_indices_count, 10);

					estimates_.SetUp(dim_indices_count);
					IPFPreparation(reduced_data);
				}

				if(precision_ > precision_aim){
					std::string reduction_descriptor;
					auto step_multiplier = ReductionInfo(reduction_descriptor, dim_indices_count);

					while(true){ // Break conditions at the end
						// Do Iterative Proportional Fitting
						while(steps_ <= kNumMargins*max_iterations*step_multiplier && precision_ > precision_aim && !error_during_fitting){
							// Update estimates
							IPFIteration(reduced_data);

							// Print information to screen
							if(!error_during_fitting && 0 == steps_%(50*step_multiplier)){
								print_mutex.lock();
								printInfo << descriptor << reduction_descriptor << " step " << steps_+old_steps << ": Current precision " << precision_*100 << "%" << std::endl;
								print_mutex.unlock();
							}
						}

						do{
							// Update estimates
							IPFConfirmation(reduced_data);

							// Print information to screen
							if(!error_during_fitting && 0 == steps_%(50*step_multiplier)){
								print_mutex.lock();
								printInfo << "Confirming " << descriptor << reduction_descriptor <<" step " << steps_+old_steps << ": Current precision " << precision_*100 << "%" << std::endl;
								print_mutex.unlock();
							}
						}while(steps_ <= kNumMargins*max_iterations*step_multiplier && precision_ > precision_aim && !error_during_fitting);

						// Finishing steps (for this expansion level) and preparation of next round
						if(steps_ > kNumMargins*max_iterations*step_multiplier){
							break;
						}

						if( data.Expand(reduced_data, dim_indices_reduced_, dim_indices_count, expansion_indices, expansion_count, N) ){
							estimates_.Expand(expansion_indices,expansion_count,&dim_indices_count);

							// Prepare next round
							auto step_multiplier_old = step_multiplier;
							step_multiplier = ReductionInfo(reduction_descriptor, dim_indices_count);

							// Adjust steps
							steps_ = steps_*step_multiplier/step_multiplier_old;
							old_steps = old_steps*step_multiplier/step_multiplier_old;
							needed_updates_ = needed_updates_*step_multiplier/step_multiplier_old;

							// Run every margin once to get the precision against the expanded data
							margin_precision_.fill(std::numeric_limits<double>::max());
							RunAllMargins( reduced_data );

							// Print information to screen
							if(!error_during_fitting && steps_%(50*step_multiplier) < (steps_-kNumMargins)%(50*step_multiplier)){
								print_mutex.lock();
								printInfo << descriptor << reduction_descriptor <<" step " << steps_+old_steps << ": Current precision " << precision_*100 << "%" << std::endl;
								print_mutex.unlock();
							}
						}
						else{
							// We already reached completely expanded data
							break;
						}
					}
				}

				if(steps_){
					estimates_.Normalize();

					// Print status
					if(!error_during_fitting){
						// Check if some result got NaN
						if( !CheckConsistency(descriptor, print_mutex) ){
							error_during_fitting = true;
						}
						else{
							precision_improved = true;

							steps_ += old_steps;
							print_mutex.lock();
							if(precision_ > precision_aim){
								printWarn << descriptor << " did not reach precision aim: " << precision_*100 << "%" << std::endl;
							}
							else{
								printInfo << "Finished iterative proportional fitting for " << descriptor << " after step " << steps_ << " with precision " << precision_*100 << "%" << std::endl;
							}
							print_mutex.unlock();
						}
					}

					update_factors_.shrink_to_fit();
				}
				else{
					steps_ = old_steps;
				}
			}
		};

		template<> template<uintMarginId S, uint16_t UPDATE> inline void LogIPF<3>::IPFStepCallerTemp( const DataStorage<3> &data ){
			switch(S){
			case 0:
				IPFUpdateOrCheck<UPDATE,0,1,2>(data);
				break;
			case 1:
				IPFUpdateOrCheck<UPDATE,0,2,1>(data);
				break;
			case 2:
				IPFUpdateOrCheck<UPDATE,1,2,0>(data);
				break;
			default:
				throw std::invalid_argument( "Received impossible step in IPFStepCaller" );
			}
		}

		template<> template<uintMarginId S, uint16_t UPDATE> inline void LogIPF<4>::IPFStepCallerTemp( const DataStorage<4> &data ){
			switch(S){
			case 0:
				IPFUpdateOrCheck<UPDATE,0,1,2,3>(data);
				break;
			case 1:
				IPFUpdateOrCheck<UPDATE,0,2,1,3>(data);
				break;
			case 2:
				IPFUpdateOrCheck<UPDATE,0,3,1,2>(data);
				break;
			case 3:
				IPFUpdateOrCheck<UPDATE,1,2,0,3>(data);
				break;
			case 4:
				IPFUpdateOrCheck<UPDATE,1,3,0,2>(data);
				break;
			case 5:
				IPFUpdateOrCheck<UPDATE,2,3,0,1>(data);
				break;
			default:
				throw std::invalid_argument( "Received impossible step in IPFStepCaller" );
			}
		}

		template<> template<uintMarginId S, uint16_t UPDATE> inline void LogIPF<5>::IPFStepCallerTemp( const DataStorage<5> &data ){
			switch(S){
			case 0:
				IPFUpdateOrCheck<UPDATE,0,1,2,3,4>(data);
				break;
			case 1:
				IPFUpdateOrCheck<UPDATE,0,2,1,3,4>(data);
				break;
			case 2:
				IPFUpdateOrCheck<UPDATE,0,3,1,2,4>(data);
				break;
			case 3:
				IPFUpdateOrCheck<UPDATE,0,4,1,2,3>(data);
				break;
			case 4:
				IPFUpdateOrCheck<UPDATE,1,2,0,3,4>(data);
				break;
			case 5:
				IPFUpdateOrCheck<UPDATE,1,3,0,2,4>(data);
				break;
			case 6:
				IPFUpdateOrCheck<UPDATE,1,4,0,2,3>(data);
				break;
			case 7:
				IPFUpdateOrCheck<UPDATE,2,3,0,1,4>(data);
				break;
			case 8:
				IPFUpdateOrCheck<UPDATE,2,4,0,1,3>(data);
				break;
			case 9:
				IPFUpdateOrCheck<UPDATE,3,4,0,1,2>(data);
				break;
			default:
				throw std::invalid_argument( "Received impossible step in IPFStepCaller" );
			}
		}
	}

	class ProbabilityEstimates{
	private:
		enum IPFDataSelector{
			kIPFQuality,
			kIPFSequenceQuality,
			kIPFBaseCall,
			kIPFDominantError,
			kIPFErrorRate,
			kIPFInDels
		};

		struct IPFThreadParams{
			IPFDataSelector selected_data;
			uintTempSeq template_segment;
			uintTileId tile_id;
			uintBaseCall ref_base;
			uintBaseCall dom_error;
			uintBaseCall last_ref_base;
		};

		// Mutex
		std::mutex print_mutex_;

		// bookkeeping
		std::atomic<uintNumFits> current_param_;
		std::atomic<bool> error_during_fitting_;
		std::atomic<bool> precision_improved_;

		uint64_t stats_creation_time_; // Store and redo fits if stats file has been updated

		// ipf calc
		std::array<std::vector<std::array<ProbabilityEstimatesSubClasses::LogIPF<5>,4>>, 2> quality_; // quality_[template_segment][tile_id][ref_base] : quality, sequence quality, previous quality, position, error rate
		std::array<std::vector<ProbabilityEstimatesSubClasses::LogIPF<4>>, 2> sequence_quality_; // sequence_quality_[template_segment][tile_id]: sequence quality, gc, mean error rate, fragment length
		std::array<std::vector<std::array<std::array<ProbabilityEstimatesSubClasses::LogIPF<5>,5>,4>>, 2> base_call_; // base_call_[template_segment][tile_id][ref_base][domError] : called base, quality, position, error number, error rate
		std::array<std::array<std::array<ProbabilityEstimatesSubClasses::LogIPF<4>,5>,5>,4> dom_error_; // dom_error_[ref_base][last_ref_base][dom_last5]: dominant error, distance, gc, start rate
		std::array<std::array<ProbabilityEstimatesSubClasses::LogIPF<4>,5>,4> error_rate_; // error_rate_[ref_base][dom_error]: error rate, distance, gc, start rate
		std::array<std::array<ProbabilityEstimatesSubClasses::LogIPF<4>,6>, 2> indels_; // indels_[Insertion/DeletionBefore][PreviousRegularCall]: indel, indel pos/length, position, gc

		// ipf results
		std::array<std::vector<std::array<ProbabilityEstimatesSubClasses::LogArrayResult<5>,4>>, 2> quality_result_; // quality_result_[template_segment][tile_id][ref_base] : quality, sequence quality, previous quality, position, error rate
		std::array<std::vector<ProbabilityEstimatesSubClasses::LogArrayResult<4>>, 2> sequence_quality_result_; // sequence_quality_result_[template_segment][tile_id]: sequence quality, gc, mean error rate, fragment length
		std::array<std::vector<std::array<std::array<ProbabilityEstimatesSubClasses::LogArrayResult<5>,5>,4>>, 2> base_call_result_; // base_call_result_[template_segment][tile_id][ref_base][domError] : called base, quality, position, error number, error rate
		std::array<std::array<std::array<ProbabilityEstimatesSubClasses::LogArrayResult<4>,5>,5>,4> dom_error_result_; // dom_error_result_[ref_base][last_ref_base][dom_last5]: dominant error, distance, gc, start rate
		std::array<std::array<ProbabilityEstimatesSubClasses::LogArrayResult<4>,5>,4> error_rate_result_; // error_rate_result_[ref_base][dom_error]: error rate, distance, gc, start rate
		std::array<std::array<ProbabilityEstimatesSubClasses::LogArrayResult<4>,6>, 2> indels_result_; // indels_result_[Insertion/DeletionBefore][PreviousRegularCall]: indel, indel pos/length, position, gc

		void SetVariables(uintTileId num_tiles);

		inline const Vect<SeqQualityStats<uintMatrixCount>> *DefineMarginsQuality(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 10 > &margins, uintTempSeq template_segment, uintTileId tile_id, uintBaseCall ref_base){
			// Dimension order: quality, sequence quality, previous quality, position, error rate
			margins.at(0) = { &stats.Qualities().BaseQualityForSequenceQualityReference(template_segment, tile_id, ref_base), true };
			margins.at(1) = { &stats.Qualities().BaseQualityForPrecedingQualityReference(template_segment, tile_id, ref_base), true };
			margins.at(3) = { &stats.Qualities().BaseQualityForErrorRateReference(template_segment, tile_id, ref_base), true };
			margins.at(4) = { &stats.Qualities().PrecedingQualityForSequenceQualityReference(template_segment, tile_id, ref_base), false };
			margins.at(5) = { &stats.Qualities().SequenceQualityForPositionReference(template_segment, tile_id, ref_base), true };
			margins.at(6) = { &stats.Qualities().SequenceQualityForErrorRateReference(template_segment, tile_id, ref_base), true };
			margins.at(7) = { &stats.Qualities().PrecedingQualityForPositionReference(template_segment, tile_id, ref_base), true };
			margins.at(8) = { &stats.Qualities().PrecedingQualityForErrorRateReference(template_segment, tile_id, ref_base), true };
			margins.at(9) = { &stats.Qualities().ErrorRateForPositionReference(template_segment, tile_id, ref_base), false };

			margins.at(2) = { NULL, true };
			return &stats.Qualities().BaseQualityStatsReference(template_segment, tile_id, ref_base); // alternative_margin
		}
		inline const Vect<SeqQualityStats<uintMatrixCount>> *DefineMarginsSequenceQuality(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > &margins, uintTempSeq template_segment, uintTileId tile_id){
			// Dimension order: sequence quality, gc, mean error rate, fragment length
			margins.at(1) = { &stats.Qualities().SequenceQualityMeanForMeanErrorRatePerTileReference(template_segment, tile_id), true };
			margins.at(2) = { &stats.Qualities().SequenceQualityMeanForFragmentLengthPerTileReference(template_segment, tile_id), true };
			margins.at(3) = { &stats.Qualities().MeanErrorRateForGCPerTileReference(template_segment, tile_id), false };
			margins.at(4) = { &stats.Qualities().GCForFragmentLengthPerTileReference(template_segment, tile_id), true };
			margins.at(5) = { &stats.Qualities().MeanErrorRateForFragmentLengthPerTileReference(template_segment, tile_id), true };

			margins.at(0) = { NULL, true };
			return &stats.Qualities().SequenceQualityMeanForGCPerTileReference(template_segment, tile_id);
		}
		inline const Vect<SeqQualityStats<uintMatrixCount>> *DefineMarginsBaseCall(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 10 > &margins, uintTempSeq template_segment, uintTileId tile_id, uintBaseCall ref_base, uintBaseCall dom_error){
			// Dimension order: called base, quality, position, error number, error rate
			margins.at(0) = { &stats.Errors().CalledBasesByBaseQuality(template_segment, tile_id, ref_base, dom_error), false};
			margins.at(1) = { &stats.Errors().CalledBasesByPosition(template_segment, tile_id, ref_base, dom_error), false};
			margins.at(2) = { &stats.Errors().CalledBasesByErrorNum(template_segment, tile_id, ref_base, dom_error), false};
			margins.at(3) = { &stats.Errors().CalledBasesByErrorRate(template_segment, tile_id, ref_base, dom_error), false};
			margins.at(5) = { &stats.Errors().ErrorNumByBaseQuality(template_segment, tile_id, ref_base, dom_error), true};
			margins.at(6) = { &stats.Qualities().BaseQualityForErrorRateReference(template_segment, tile_id, ref_base, dom_error), true };
			margins.at(7) = { &stats.Errors().ErrorNumByPosition(template_segment, tile_id, ref_base, dom_error), true};
			margins.at(8) = { &stats.Qualities().ErrorRateForPositionReference(template_segment, tile_id, ref_base, dom_error), false };
			margins.at(9) = { &stats.Errors().ErrorNumByErrorRate(template_segment, tile_id, ref_base, dom_error), false};

			margins.at(4) = { NULL, true };
			return &stats.Qualities().BaseQualityStatsReference(template_segment, tile_id, ref_base, dom_error);
		}
		inline void DefineMarginsDominantError(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > &margins, uintTempSeq ref_base, uintBaseCall last_ref_base, uintBaseCall dom_last5){
			// Dimension order: dominant error, distance, gc, start rate
			margins.at(0) = { &stats.Coverage().DominantErrorsByDistance(ref_base, last_ref_base, dom_last5), true };
			margins.at(1) = { &stats.Coverage().DominantErrorsByGC(ref_base, last_ref_base, dom_last5), true };
			margins.at(2) = { &stats.Coverage().DominantErrorsByStartRates(ref_base, last_ref_base, dom_last5), true };
			margins.at(3) = { &stats.Coverage().GCByDistance(ref_base, last_ref_base, dom_last5), false };
			margins.at(4) = { &stats.Coverage().StartRatesByDistance(ref_base, last_ref_base, dom_last5), false };
			margins.at(5) = { &stats.Coverage().StartRatesByGC(ref_base, last_ref_base, dom_last5), false };
		}

		inline void DefineMarginsErrorRate(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > &margins, uintBaseCall ref_base, uintBaseCall dom_error){
			margins.at(0) = { &stats.Coverage().ErrorRatesByDistance(ref_base, dom_error), true };
			margins.at(1) = { &stats.Coverage().ErrorRatesByGC(ref_base, dom_error), true };
			margins.at(2) = { &stats.Coverage().ErrorRatesByStartRates(ref_base, dom_error), true };
			margins.at(3) = { &stats.Coverage().GCByDistance(ref_base, dom_error), false };
			margins.at(4) = { &stats.Coverage().StartRatesByDistance(ref_base, dom_error), false };
			margins.at(5) = { &stats.Coverage().StartRatesByGC(ref_base, dom_error), false };
		}
		inline void DefineMarginsIndels(const DataStats &stats, std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > &margins, uintInDelType indel_type, uintBaseCall last_call){
			margins.at(0) = { &stats.Errors().InDelByInDelPos(indel_type, last_call), true };
			margins.at(1) = { &stats.Errors().InDelByPosition(indel_type, last_call), true };
			margins.at(2) = { &stats.Errors().InDelByGC(indel_type, last_call), true };
			margins.at(3) = { &stats.Errors().InDelPosByPosition(indel_type, last_call), true };
			margins.at(4) = { &stats.Errors().InDelPosByGC(indel_type, last_call), true };
			margins.at(5) = { &stats.Errors().GCByPosition(indel_type, last_call), false };
		}

		template<size_t N> void UpdateMarginDefinitionToTotal(std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, N > &margins, const std::array< Vect<Vect<uintMatrixCount>>, N> &margin_sums){
			for(auto n=N; n--; ){
				margins.at(n).first = &(margin_sums.at(n));
				margins.at(n).second = true;
			}
		}

		template<size_t N> void AddMarginsFromDefinition(std::array< Vect<Vect<uintMatrixCount>>, N> &margin_sums, const std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, N > &margins, const Vect<SeqQualityStats<uintMatrixCount>> *alternative_margin){
			size_t N2(0);
			switch(N){
			case 3:
				N2=3;
				break;
			case 6:
				N2=4;
				break;
			case 10:
				N2=5;
				break;
			}

			uintMarginId dim_a(N2), dim_b(N2-1);
			for(uintMarginId n=N; n--; ){
				if(--dim_a == dim_b){
					--dim_b;
					dim_a = N2-1;
				}

				uintMatrixIndex index_first_from, index_first_to;
				if(margins.at(n).first){
					index_first_from = (*margins.at(n).first).from();
					index_first_to = (*margins.at(n).first).to();
				}
				else{
					index_first_from = (*alternative_margin).from();
					index_first_to = (*alternative_margin).to();
				}

				for(uintMatrixIndex index_first = index_first_from; index_first < index_first_to; ++index_first ){
					uintMatrixIndex index_second_from, index_second_to;
					if(margins.at(n).first){
						index_second_from = (*margins.at(n).first).at(index_first).from();
						index_second_to = (*margins.at(n).first).at(index_first).to();
					}
					else{
						index_second_from = (*alternative_margin).at(index_first).from();
						index_second_to = (*alternative_margin).at(index_first).to();
					}

					for(uintMatrixIndex index_second = index_second_from; index_second < index_second_to; ++index_second ){
						uintMatrixIndex ia, ib;

						if(margins.at(n).second){
							// Indices are flipped (like they are in data_)
							ia = index_first;
							ib = index_second;
						}
						else{
							ia = index_second;
							ib = index_first;
						}

						if(margins.at(n).first){
							margin_sums.at(n).at(ia).at(ib) += (*margins.at(n).first).at(index_first).at(index_second);
						}
						else{
							margin_sums.at(n).at(ia).at(ib) += (*alternative_margin).at(index_first).at(index_second);
						}
					}
				}
			}
		}

		void IterativeProportionalFitting(const DataStats &stats, IPFDataSelector selected_data, uintTempSeq template_segment, uintTileId tile_id, uintBaseCall ref_base, uintBaseCall dom_error, uintBaseCall last_ref_base, uintNumFits max_iterations, double precision_aim);
		static void IPFThread(ProbabilityEstimates &self, const DataStats &stats, const std::vector<IPFThreadParams> &params, uintNumFits max_iterations, double precision_aim);

		// boost serialization
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & stats_creation_time_;
			ar & quality_;
			ar & sequence_quality_;
			ar & base_call_;
			ar & dom_error_;
			ar & error_rate_;
			ar & indels_;
		}

		// Google test
		friend class ProbabilityEstimatesTest;

	public:
		ProbabilityEstimates():
			current_param_(0),
			error_during_fitting_(false),
			precision_improved_(false){
		}

		void PrepareResult();

		const ProbabilityEstimatesSubClasses::LogArrayResult<5> &Quality(uintTempSeq template_segment, uintTileId tile_id, uintBaseCall base) const{
			return quality_result_.at(template_segment).at(tile_id).at(base);
		}
		const ProbabilityEstimatesSubClasses::LogArrayResult<4> &SequenceQuality(uintTempSeq template_segment, uintTileId tile_id) const{
			return sequence_quality_result_.at(template_segment).at(tile_id);
		}
		const ProbabilityEstimatesSubClasses::LogArrayResult<5> &BaseCall(uintTempSeq template_segment, uintTileId tile_id, uintBaseCall ref_base, uintBaseCall dom_error) const{
			return base_call_result_.at(template_segment).at(tile_id).at(ref_base).at(dom_error);
		}
		const ProbabilityEstimatesSubClasses::LogArrayResult<4> &DominantError(uintBaseCall base, uintBaseCall prev_base, uintBaseCall dom_last5) const{
			return dom_error_result_.at(base).at(prev_base).at(dom_last5);
		}
		const ProbabilityEstimatesSubClasses::LogArrayResult<4> &ErrorRate(uintBaseCall base, uintBaseCall dom_last5) const{
			return error_rate_result_.at(base).at(dom_last5);
		}
		const ProbabilityEstimatesSubClasses::LogArrayResult<4> &InDels(uintInDelType indel_type, uintBaseCall last_call) const{
			return indels_result_.at(indel_type).at(last_call);
		}

		void ChangeErrorRate(double error_multiplier){
			printInfo << "Multiplying error rate with " << error_multiplier << std::endl;
			for(auto template_segment=base_call_result_.size(); template_segment--; ){
				for(auto tile_id=base_call_result_.at(template_segment).size(); tile_id--; ){
					for(auto ref_base=base_call_result_.at(template_segment).at(tile_id).size(); ref_base--; ){
						for(auto dom_error=base_call_result_.at(template_segment).at(tile_id).at(ref_base).size(); dom_error--; ){
							base_call_result_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).ModifyPar0(ref_base, 1.0/error_multiplier); // Didive the probability of the correct base by error_multiplier
						}
					}
				}
			}
		}

		bool Load( const char *archive_file );
		bool Save( const char *archive_file ) const;

		bool Estimate(const DataStats &stats, uintNumFits max_iterations, double precision_aim, uintNumThreads num_threads, const char *output, const char *input="");
	};

}

#endif // PROBABILITYESTIMATES_H
