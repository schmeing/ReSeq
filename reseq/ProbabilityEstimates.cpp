#include "ProbabilityEstimates.h"
using reseq::ProbabilityEstimatesSubClasses::DataStorage;
using reseq::ProbabilityEstimatesSubClasses::LogArrayCalc;
using reseq::ProbabilityEstimatesSubClasses::MaxPrecision;
using reseq::ProbabilityEstimates;
using reseq::SeqQualityStats;
using reseq::Vect;

//include <algorithm>
using std::max_element;
//include<array>
using std::array;
#include <cmath>
using std::log;
#include <exception>
using std::exception;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iterator>
using std::distance;
//include <limits>
using std::numeric_limits;
#include <random>
using std::mt19937_64;
using std::uniform_int_distribution;
#include <sstream>
using std::ostringstream;
using std::stringstream;
#include <stdio.h>
//include <string>
using std::string;
#include <thread>
using std::thread;
//include <utility>
using std::pair;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

#include <iomanip>
using std::setw;
//include <algorithm>
using std::min;
using std::max;

//include "utilities.hpp"
using reseq::uintMarginId;
using reseq::utilities::CreateDir;
using reseq::utilities::DeleteFile;
using reseq::utilities::FileExists;

template<uintMarginId N> void DataStorage<N>::CombineEndBinsToMinCount(array< vector<uintMatrixIndex>, N > &dim_indices, array<Vect<uintMatrixCount>, N> &dim_control){
	// Start with summing up from both ends until kMinCountsPer1dBin is reached for the outermost bins
	array<uintMatrixIndex, N> first_bin, last_bin;
	first_bin.at(0) = dim_control.at(0).from();
	last_bin.at(0) = dim_control.at(0).to()-1;
	for( uintMarginId cur_dim = 1; cur_dim < dim_control.size(); ++cur_dim){ // Do not do change first dimension!!!
		// Get the bin, where everything coming after will be added to (to reach kMinCountsPer1dBin) and update dim_control in the process
		last_bin.at(cur_dim) = dim_control.at(cur_dim).to()-1;
		while(dim_control.at(cur_dim).at(last_bin.at(cur_dim)) < kMinCountsPer1dBin && last_bin.at(cur_dim) > dim_control.at(cur_dim).from()){
			dim_control.at(cur_dim).at(last_bin.at(cur_dim)-1) += dim_control.at(cur_dim).at(last_bin.at(cur_dim));
			dim_control.at(cur_dim).at(last_bin.at(cur_dim)) = 0;
			--last_bin.at(cur_dim);
		}

		// Get the bin, where everything before will be added to (to reach kMinCountsPer1dBin)
		first_bin.at(cur_dim) = dim_control.at(cur_dim).from();
		while(dim_control.at(cur_dim).at(first_bin.at(cur_dim)) < kMinCountsPer1dBin && first_bin.at(cur_dim) < last_bin.at(cur_dim)){
			dim_control.at(cur_dim).at(first_bin.at(cur_dim)+1) += dim_control.at(cur_dim).at(first_bin.at(cur_dim));
			dim_control.at(cur_dim).at(first_bin.at(cur_dim)) = 0;
			++first_bin.at(cur_dim);
		}

		// Shrink and Update dim_control so it can be used later in CombineBins
		dim_control.at(cur_dim).SoftShrink(); // Space will be freed soon anyway as it is just a temporary variable, so don't do it here
		dim_control.at(cur_dim).Shift(-3000000); // Set offset to zero to really remove the elements in the beginning

		// Reduce dim_indices to the new limits
		dim_indices.at(cur_dim).resize(last_bin.at(cur_dim)+1);
		dim_indices.at(cur_dim).erase(dim_indices.at(cur_dim).begin(), dim_indices.at(cur_dim).begin()+first_bin.at(cur_dim));
	}

	// Reduce data_ to the new limits (This part does not have to be reverted later, as we anyways take the outermost bins for everything that is further out)
	uintMarginId dim_a(N), dim_b(N-1);
	for(uintMarginId n=kNumMargins; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = N-1;
		}

		// Add values from bins outside the new limits, to first or last bin (dim_a)
		for(uintMatrixIndex ia = 0; ia < first_bin.at(dim_a); ++ia ){
			for(uintMatrixIndex ib = 0; ib < dim_size_.at(dim_b); ++ib ){
				data_.at(n).at(first_bin.at(dim_a)*dim_size_.at(dim_b)+ib) += data_.at(n).at(ia*dim_size_.at(dim_b)+ib);
			}
		}
		for(uintMatrixIndex ia = dim_size_.at(dim_a); --ia > last_bin.at(dim_a); ){
			for(uintMatrixIndex ib = 0; ib < dim_size_.at(dim_b); ++ib ){
				data_.at(n).at(last_bin.at(dim_a)*dim_size_.at(dim_b)+ib) += data_.at(n).at(ia*dim_size_.at(dim_b)+ib);
			}
		}
		// Add values from bins outside the new limits, to first or last bin (dim_b)
		for(uintMatrixIndex ib = 0; ib < first_bin.at(dim_b); ++ib ){
			for(uintMatrixIndex ia = first_bin.at(dim_a); ia <= last_bin.at(dim_a); ++ia ){ // We can ignore the other indexes as we already added their values
				data_.at(n).at(ia*dim_size_.at(dim_b)+first_bin.at(dim_b)) += data_.at(n).at(ia*dim_size_.at(dim_b)+ib);
			}
		}
		for(uintMatrixIndex ib = dim_size_.at(dim_b); --ib > last_bin.at(dim_b); ){
			for(uintMatrixIndex ia = first_bin.at(dim_a); ia <= last_bin.at(dim_a); ++ia ){ // We can ignore the other indexes as we already added their values
				data_.at(n).at(ia*dim_size_.at(dim_b)+last_bin.at(dim_b)) += data_.at(n).at(ia*dim_size_.at(dim_b)+ib);
			}
		}

		// Remove bins outside the new limits
		if(dim_indices.at(dim_b).size() != dim_size_.at(dim_b)){
			// Inner dimension has changed, so instead of deleting small bits here and there, we start from the beginning and copy the new value there and finally cut the rest
			for(uintMatrixIndex ia = 0; ia < dim_indices.at(dim_a).size(); ++ia ){
				for(uintMatrixIndex ib = 0; ib < dim_indices.at(dim_b).size(); ++ib ){
					data_.at(n).at(ia*dim_indices.at(dim_b).size()+ib) = data_.at(n).at((ia+first_bin.at(dim_a))*dim_size_.at(dim_b)+ib+first_bin.at(dim_b));
				}
			}
			data_.at(n).resize(dim_indices.at(dim_a).size() * dim_indices.at(dim_b).size());
		}
		else if(dim_indices.at(dim_a).size() != dim_size_.at(dim_a)){
			// Only the outer dimension changed so we can just remove the bins at the beginning and end
			data_.at(n).resize( (first_bin.at(dim_a) + dim_indices.at(dim_a).size()) * dim_size_.at(dim_b));
			data_.at(n).erase( data_.at(n).begin(), data_.at(n).begin() + first_bin.at(dim_a)*dim_size_.at(dim_b) );
		}
	}

	// Update dim_size_ with the new limits (Still needed to reduce data, so only done afterwards)
	for(uintMarginId cur_dim=N; cur_dim--;){
		dim_size_.at(cur_dim) = dim_indices.at(cur_dim).size();
	}
}

template<uintMarginId N> void DataStorage<N>::CombineBins(array<vector<uintMatrixIndex>, N> &initial_dim_indices_reduced, array<Vect<uintMatrixCount>, N> &dim_control){
	array<vector<uintMatrixIndex>, N> dim_indices_count;

	for( uintMarginId cur_dim = 0; cur_dim < N; ++cur_dim ){
		initial_dim_indices_reduced.at(cur_dim).clear();
		initial_dim_indices_reduced.at(cur_dim).reserve(dim_control.at(cur_dim).size());
		dim_indices_count.at(cur_dim).reserve(dim_control.at(cur_dim).size());
	}

	// First dimension remains unchanged as we later read out this values based on all the other dimensions
	for(uintMatrixIndex bin = 0; bin < dim_control.at(0).size(); ++bin){
		initial_dim_indices_reduced.at(0).push_back(bin);
		dim_indices_count.at(0).push_back(1);
	}

	// Handle other dimensions
	for( uintMarginId cur_dim = 1; cur_dim < N; ++cur_dim ){ // Do not do anything for first dimension
		// First and Last bin are above kMinCountsPer1dBin already from call to CombineEndBinsToMinCount, except if only one element is left in the vector
		if( 1== dim_control.at(cur_dim).size() ){
			// Counts in the single bin could be lower than kMinCountsPer1dBin, but we cannot do anything about it
			initial_dim_indices_reduced.at(cur_dim).push_back(0);
			dim_indices_count.at(cur_dim).push_back(1);
		}
		else{
			// Going from left to right combine bins with a count lower than kMinCountsPer1dBin with the neighbouring bin that has fewer counts
			uintMatrixIndex new_bin(0), combined_bins(1);
			for(uintMatrixIndex old_bin = 0; old_bin < dim_control.at(cur_dim).size(); ++old_bin){
				if( kMinCountsPer1dBin > dim_control.at(cur_dim).at(old_bin) ){
					if( dim_control.at(cur_dim).at(old_bin+1) > dim_control.at(cur_dim).at(new_bin-1) ){ // old_bin+1 and new_bin-1 guaranteed to be valid as CombineEndBinsToMinCount called earlier makes sure first and last bin are above kMinCountsPer1dBin and don't end up here
						// Next bin bigger, so add to previous bin
						for(auto i=combined_bins; i--; ){
							initial_dim_indices_reduced.at(cur_dim).push_back(new_bin-1);
						}
						dim_indices_count.at(cur_dim).at(new_bin-1) += combined_bins;
						dim_control.at(cur_dim).at(new_bin-1) += dim_control.at(cur_dim).at(old_bin);

						// Reset combined_bins as we have integrated it into dim_indices_count
						combined_bins = 1;
					}
					else{
						// Previous bin bigger, so add to next bin
						dim_control.at(cur_dim).at(old_bin+1) += dim_control.at(cur_dim).at(old_bin);
						++combined_bins; // Next bin now contains one more bin
					}
				}
				else{
					// Insert bin
					for(auto i=combined_bins; i--; ){
						initial_dim_indices_reduced.at(cur_dim).push_back(new_bin);
					}
					dim_indices_count.at(cur_dim).push_back(combined_bins);
					dim_control.at(cur_dim).at(new_bin) = dim_control.at(cur_dim).at(old_bin); // Update control bin on the fly (Everything before new_bin is the new stuff, everything after old_bin is still untouched and everything in between is junk)

					// Prepare for next bin
					++new_bin;
					combined_bins = 1;
				}
			}
			dim_control.at(cur_dim).PopBack( dim_control.at(cur_dim).size() - dim_indices_count.at(cur_dim).size() );
		}
	}

	DataStorage<N> reduced_data;
	Reduce(reduced_data, initial_dim_indices_reduced, dim_indices_count);

	// Reduce further if necessary, so that we end up below kMaxBinsPerDimension
	array<vector<uintMatrixIndex>, N> dim_indices_reduced;
	reduced_data.Reduce(*this, dim_indices_reduced, dim_indices_count, kMaxBinsPerDimension);

	// Integrate dim_indices_reduced from second reduction into initial_dim_indices_reduced from first reduction
	CombineDimIndices(initial_dim_indices_reduced, dim_indices_reduced);
}

template<uintMarginId N> void DataStorage<N>::FillDiffMatrix(vector<double> &diff_matrix, uintMarginId dimension) const{
	// Fill diff_matrix
	diff_matrix.clear();
	diff_matrix.resize(dim_size_.at(dimension)*dim_size_.at(dimension), 0.0);

	uintMarginId dim_a(N), dim_b(N-1);
	for(uintMarginId n=kNumMargins; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = N-1;
		}

		if(dimension == dim_a || dimension == dim_b){
			double value1, value2, sum_diff;
			for(uintMatrixIndex ind1=0; ind1 < dim_size_.at(dimension); ind1++){
				for(uintMatrixIndex ind2=ind1+1; ind2 < dim_size_.at(dimension); ind2++){
					sum_diff = 0.0;
					if(dimension == dim_a){
						for(uintMatrixIndex j=0; j < dim_size_.at(dim_b); ++j){
							value1 = data_.at(n).at(ind1*dim_size_.at(dim_b)+j);
							value2 = data_.at(n).at(ind2*dim_size_.at(dim_b)+j);

							if(value1 || value2){
								sum_diff += std::abs(value1-value2)/(value1+value2);
							}
						}
					}
					else{
						for(uintMatrixIndex j=0; j < dim_size_.at(dim_a); ++j){
							value1 = data_.at(n).at(j*dim_size_.at(dim_b)+ind1);
							value2 = data_.at(n).at(j*dim_size_.at(dim_b)+ind2);

							if(value1 || value2){
								sum_diff += std::abs(value1-value2)/(value1+value2);
							}
						}
					}

					diff_matrix.at(ind1*dim_size_.at(dimension)+ind2) += sum_diff/dim_size_.at((dimension == dim_a?dim_b:dim_a));
				}
			}
		}
	}
}

template<uintMarginId N> double DataStorage<N>::MaxDiff(const vector<double> &diff_matrix, uintMatrixIndex combined_index, vector<uintMatrixIndex> &same_combined_index, const vector<uintMatrixIndex> &dim_indices_reduced, const vector<uintMatrixIndex> &dim_indices_count) const{
	double max_diff(0.0);

	if(1 < dim_indices_count.at(combined_index)){ // We need a minimum of 2 values to have a difference
		same_combined_index.clear();
		for(uintMatrixIndex current_index = dim_indices_reduced.size(); current_index--;){
			if(combined_index == dim_indices_reduced.at(current_index)){
				for(auto comparison_index : same_combined_index){
					// We go backwards so all comparison_index are higher than current_index which is important as diff_matrix is triangular
					if( diff_matrix.at(current_index*dim_indices_reduced.size()+comparison_index) > max_diff){
						max_diff = diff_matrix.at(current_index*dim_indices_reduced.size()+comparison_index);
					}
				}

				same_combined_index.push_back(current_index);
			}
		}
	}

	return max_diff;
}

template<uintMarginId N> bool DataStorage<N>::SetUp(const array< pair< const Vect<Vect<uintMatrixCount>> *, bool>, N*(N-1)/2 > &margins, const Vect<SeqQualityStats<uintMatrixCount>> *alternative_margin, array< vector<uintMatrixIndex>, N > &dim_indices, array<vector<uintMatrixIndex>, N> &initial_dim_indices_reduced, std::mutex &print_mutex){
	// Make sure dim_indices are empty
	for(uintMarginId n=N; n--;){
		dim_indices.at(n).clear();
	}

	// Fill the existing indices(qualities, positions, etc.) for each dimension into vector
	uintMarginId dim_a(N), dim_b(N-1);
	for(uintMarginId n=kNumMargins; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = N-1;
		}

		uintMarginId dim_first;
		if(margins.at(n).second){
			// Indices are flipped (like they are in data_)
			dim_first = dim_a;
		}
		else{
			dim_first = dim_b;
		}

		if( !dim_indices.at(dim_first).size() ){
			if(margins.at(n).first){
				GetIndicesFirstDimension(*margins.at(n).first, dim_indices.at(dim_first));
			}
			else{
				GetIndicesFirstDimension(*alternative_margin, dim_indices.at(dim_first));
			}
		}

		if( 0 == dim_b ){
			// It's time to also use the second dimension if a vector where this dimension is first, does not exist
			if( !dim_indices.at(dim_a).size() ){
				if(margins.at(n).first){
					GetIndicesSecondDimension(*margins.at(n).first, dim_indices.at(dim_a));
				}
				else{
					GetIndicesSecondDimension(*alternative_margin, dim_indices.at(dim_a));
				}
			}
			else if(0 == n && !dim_indices.at(0).size()){
				if(margins.at(0).first){
					GetIndicesSecondDimension(*margins.at(0).first, dim_indices.at(0));
				}
				else{
					GetIndicesSecondDimension(*alternative_margin, dim_indices.at(0));
				}
			}
		}
	}

	// Get the number of counts in the matrix, so we can later normalize to 1
	uintMatrixCount sum_matrix;
	if(margins.at(0).first){
		sum_matrix = SumVect(*margins.at(0).first);
	}
	else{
		sum_matrix = SumVect(*alternative_margin);
	}

	// Fill in data and check that counts are distributed the same way within a variable in all margins where it appears
	for(uintMarginId cur_dim=N; cur_dim--;){
		dim_size_.at(cur_dim) = dim_indices.at(cur_dim).size();
	}
	array<Vect<uintMatrixCount>, N> dim_control;
	Vect<uintMatrixCount> dim_control_a, dim_control_b;

	dim_a = N;
	dim_b = N-1;
	for(uintMarginId n=kNumMargins; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = N-1;
		}

		// Clear the controls from last iteration
		dim_control_a.Clear();
		dim_control_b.Clear();

		// Set the vectors to the proper size
		data_.at(n).resize(dim_size_.at(dim_a)*dim_size_.at(dim_b));

		for(uintMatrixIndex ia = 0; ia < dim_size_.at(dim_a); ++ia ){
			// Fill the vectors
			for(uintMatrixIndex ib = 0; ib < dim_size_.at(dim_b); ++ib ){
				uintMatrixIndex index_first, index_second;
				if(margins.at(n).second){
					// Indices are flipped (like they are in data_)
					index_first = dim_indices.at(dim_a).at(ia);
					index_second = dim_indices.at(dim_b).at(ib);
				}
				else{
					index_first = dim_indices.at(dim_b).at(ib);
					index_second = dim_indices.at(dim_a).at(ia);
				}

				if(margins.at(n).first){
					// Set data
					data_.at(n).at(ia*dim_size_.at(dim_b)+ib) = static_cast<double>((*margins.at(n).first).at(index_first)[index_second])/sum_matrix; // Second index might not exist for all first indeces
					// Set control
					dim_control_a[ia] += (*margins.at(n).first).at(index_first)[index_second];
					dim_control_b[ib] += (*margins.at(n).first).at(index_first)[index_second];
				}
				else{
					// Set data
					data_.at(n).at(ia*dim_size_.at(dim_b)+ib) = static_cast<double>((*alternative_margin).at(index_first)[index_second])/sum_matrix; // Second index might not exist for all first indeces
					// Set control
					dim_control_a[ia] += (*alternative_margin).at(index_first)[index_second];
					dim_control_b[ib] += (*alternative_margin).at(index_first)[index_second];
				}
			}
		}

		// Check control
		if(dim_control.at(dim_a).size()){
			if(! (dim_control_a == dim_control.at(dim_a)) ){
				print_mutex.lock();
				printErr << "The dimension " << dim_a << " is inconsistent in margin " << n << ":" << std::endl;
				for(auto element : dim_control_a){
					std::cout << element << ' ';
				}
				std::cout << std::endl;
				for(auto element : dim_control.at(dim_a)){
					std::cout << element << ' ';
				}
				std::cout << std::endl;
				print_mutex.unlock();
				return false;
			}
		}
		else{
			dim_control.at(dim_a) = dim_control_a;
		}

		if(dim_control.at(dim_b).size()){
			if(! (dim_control_b == dim_control.at(dim_b)) ){
				print_mutex.lock();
				printErr << "The dimension " << dim_b << " is inconsistent in margin " << n << ":" << std::endl;
				for(auto element : dim_control_b){
					std::cout << element << ' ';
				}
				std::cout << std::endl;
				for(auto element : dim_control.at(dim_b)){
					std::cout << element << ' ';
				}
				std::cout << std::endl;
				print_mutex.unlock();
				return false;
			}
		}
		else{
			dim_control.at(dim_b) = dim_control_b;
		}
	}

	// Initial reduction of data that is only reverted after all the fitting (Reversion in LogIPF.FullExpansion): Combining 1d-bins that have less than kMinCountsPer1dBin entries and reduce data so it has kMaxBinsPerDimension (First dimension must not be changed, as this are the values we want to read out)
	CombineEndBinsToMinCount(dim_indices, dim_control);
	CombineBins(initial_dim_indices_reduced, dim_control);

	return true;
}

template<uintMarginId N> void DataStorage<N>::Reduce(DataStorage<N> &reduced_data, const array<vector<uintMatrixIndex>, N> &dim_indices_reduced, const array<vector<uintMatrixIndex>, N> &dim_indices_count) const{
	// Set the new dimension sizes
	for(uintMarginId n=N; n--; ){
		reduced_data.dim_size_.at(n) = dim_indices_count.at(n).size();
	}

	// Sum the data based on dim_indices_reduced
	uintMarginId dim_a(N), dim_b(N-1);
	for(uintMarginId n=kNumMargins; n--; ){
		if(--dim_a == dim_b){
			--dim_b;
			dim_a = N-1;
		}

		reduced_data.data_.at(n).clear();
		reduced_data.data_.at(n).reserve( this->data_.at(n).size() ); // It will become this size after expansions
		reduced_data.data_.at(n).resize( reduced_data.dim_size_.at(dim_a)*reduced_data.dim_size_.at(dim_b), 0.0 );

		for(auto ind_a=dim_indices_reduced.at(dim_a).size(); ind_a--; ){
			for(auto ind_b=dim_indices_reduced.at(dim_b).size(); ind_b--; ){
				reduced_data.data_.at(n).at( dim_indices_reduced.at(dim_a).at(ind_a) * reduced_data.dim_size_.at(dim_b) + dim_indices_reduced.at(dim_b).at(ind_b) ) += data_.at(n).at(ind_a*dim_size_.at(dim_b)+ind_b);
			}
		}
	}
}

template<uintMarginId N> void DataStorage<N>::Reduce(DataStorage<N> &reduced_data, array<vector<uintMatrixIndex>, N> &dim_indices_reduced, array<vector<uintMatrixIndex>, N> &dim_indices_count, uintMatrixIndex max_size_dim) const{
	vector<double> diff_matrix;
	vector<uintMatrixIndex> seed_indeces;

	for(auto dimension=N; dimension--; ){
		if(dim_size_.at(dimension) > max_size_dim){
			FillDiffMatrix(diff_matrix, dimension);

			// Split the max_size_dim most separate values to generate a seed
			// Start by taking the two most different values
			double max_element(0.0);
			uintMatrixIndex max_ind1(0), max_ind2(0);
			for(uintMatrixIndex ind1=0; ind1 < dim_size_.at(dimension); ind1++){
				for(uintMatrixIndex ind2=ind1+1; ind2 < dim_size_.at(dimension); ind2++){
					if(diff_matrix.at(ind1*dim_size_.at(dimension)+ind2) > max_element){
						max_element = diff_matrix.at(ind1*dim_size_.at(dimension)+ind2);
						max_ind1 = ind1;
						max_ind2 = ind2;
					}
				}
			}

			seed_indeces.resize(2);
			seed_indeces.at(0) = max_ind1;
			seed_indeces.at(1) = max_ind2;

			dim_indices_reduced.at(dimension).clear();
			dim_indices_reduced.at(dimension).resize(dim_size_.at(dimension),max_size_dim); // max_size_dim cannot be a valid index, as it is the highest reduced index+1
			for( auto red_ind = 2; red_ind--; ){
				dim_indices_reduced.at(dimension).at(seed_indeces.at(red_ind)) = red_ind;
			}

			uintMatrixIndex max_ind;
			double max_value, min_value, current_value;
			while(seed_indeces.size() < max_size_dim){
				// Add the next seed index by choosing the one with the highest minimal differences to all already chosen seed values
				max_value = 0.0;

				for(auto ind=dim_size_.at(dimension); ind--; ){
					if(max_size_dim == dim_indices_reduced.at(dimension).at(ind)){ // No reduced index has been assigned yet, so it is not a seed index already
						// Get the minimal differences to all already chosen seed values
						min_value = std::numeric_limits<double>::max();

						for( auto red_ind = seed_indeces.size(); red_ind--; ){
							if(ind > seed_indeces.at(red_ind)){
								// Entry in triangular matrix is [seed][ind]
								current_value = diff_matrix.at(seed_indeces.at(red_ind)*dim_size_.at(dimension)+ind);
							}
							else{
								// Entry in triangular matrix is [ind][seed]
								current_value = diff_matrix.at(ind*dim_size_.at(dimension)+seed_indeces.at(red_ind));
							}

							if(current_value < min_value){
								min_value = current_value;
							}
						}

						// If it is higher than for another index, take this one
						if( min_value > max_value ){
							max_value = min_value;
							max_ind = ind;
						}
					}
				}

				dim_indices_reduced.at(dimension).at(max_ind) = seed_indeces.size();
				seed_indeces.push_back(max_ind);
			}

			dim_indices_count.at(dimension).clear();
			dim_indices_count.at(dimension).reserve(dim_size_.at(dimension)); // It will become this size after expansions
			dim_indices_count.at(dimension).resize(max_size_dim, 1);

			// Join everything with the closest seed value
			for(auto ind=dim_size_.at(dimension); ind--; ){
				if(max_size_dim == dim_indices_reduced.at(dimension).at(ind)){ // No reduced index has been assigned yet, so it is not one of the seed indeces
					uintMatrixIndex min_red_ind = 0;
					double current_value, min_value(std::numeric_limits<double>::max());

					// Find closest seed index
					for( auto red_ind = max_size_dim; red_ind--; ){
						if(ind > seed_indeces.at(red_ind)){
							// Entry in triangular matrix is [seed][ind]
							current_value = diff_matrix.at(seed_indeces.at(red_ind)*dim_size_.at(dimension)+ind);
						}
						else{
							// Entry in triangular matrix is [ind][seed]
							current_value = diff_matrix.at(ind*dim_size_.at(dimension)+seed_indeces.at(red_ind));
						}

						if(current_value < min_value){
							min_value = current_value;
							min_red_ind = red_ind;
						}
					}

					// Enter into returned vectors
					dim_indices_reduced.at(dimension).at(ind) = min_red_ind;
					++dim_indices_count.at(dimension).at(min_red_ind);
				}
			}
		}
		else{
			// We don't have more than max_size_dim entries in this dimension, so we just keep the indeces
			dim_indices_reduced.at(dimension).resize(dim_size_.at(dimension));
			for( auto i = dim_size_.at(dimension); i--; ){
				dim_indices_reduced.at(dimension).at(i) = i;
			}

			dim_indices_count.at(dimension).resize(dim_size_.at(dimension), 1);
		}
	}

	// Sum the data based on dim_indices_reduced
	Reduce(reduced_data, dim_indices_reduced, dim_indices_count);
}

template<uintMarginId N> bool DataStorage<N>::Expand(DataStorage<N> &reduced_data, array<vector<uintMatrixIndex>, N> &dim_indices_reduced, array<vector<uintMatrixIndex>, N> &dim_indices_count, array<vector<uintMatrixIndex>, N> &expansion_indices, array<vector<uintMatrixIndex>, N> &expansion_count, uintMatrixIndex max_duplications) const{
	// Check need of duplications
	uintMatrixIndex duplications_needed(0);
	array<uintMatrixIndex, N> dim_new_size;

	for(auto dimension=N; dimension--; ){
		for( auto tmp_value = dim_indices_count.at(dimension).size(); tmp_value < dim_size_.at(dimension); tmp_value *= 2 ){
			++duplications_needed;
		}
		dim_new_size.at(dimension) = dim_indices_count.at(dimension).size(); // Store current sizes
	}

	if( 0 == duplications_needed ){
		return false; // No expansion necessary anymore
	}

	// Distribute duplications
	if(duplications_needed > max_duplications){
		// Distribute max_duplications equally, in ties favoring the ones that need more
		uintMatrixIndex duplications(0), dupl_dimension;

		while(duplications < max_duplications){
			dupl_dimension = N;
			for(auto dimension=N; dimension--; ){
				if(dim_new_size.at(dimension) < dim_size_.at(dimension)){
					// This dimension still needs duplication
					if( N == dupl_dimension ){ // No dimension needed duplication yet
						dupl_dimension = dimension;
					}
					else{
						// Give duplication to the one that has less duplications, in a tie to the one that needs more
						if(dim_new_size.at(dimension) < dim_new_size.at(dupl_dimension) || (dim_new_size.at(dimension) == dim_new_size.at(dupl_dimension) && dim_size_.at(dimension) > dim_size_.at(dupl_dimension))){
							dupl_dimension = dimension;
						}
					}
				}
			}
			dim_new_size.at(dupl_dimension) *= 2; // We don't need to adjust it to the real size as the algorithm later doesn't care
			++duplications;
		}
	}
	else{
		// Give all the duplications needed as we need less than max_duplications
		for(auto dimension=N; dimension--; ){
			dim_new_size.at(dimension) = dim_size_.at(dimension);
		}
	}

	// Apply duplications
	vector<double> diff_matrix, max_diff;
	vector<uintMatrixIndex> same_combined_index;
	for(auto dimension=N; dimension--; ){
		if(dim_new_size.at(dimension) >= dim_size_.at(dimension)){
			// We can expand to the original data so we do
			expansion_indices.at(dimension) = std::move(dim_indices_reduced.at(dimension));
			dim_indices_reduced.at(dimension).resize(dim_size_.at(dimension));
			for(auto i=dim_size_.at(dimension); i--; ){
				dim_indices_reduced.at(dimension).at(i) = i;
			}

			expansion_count.at(dimension) = std::move(dim_indices_count.at(dimension));
			dim_indices_count.at(dimension).resize(dim_size_.at(dimension), 1);
		}
		else{
			// Prepare new sizes
			expansion_indices.at(dimension).reserve(dim_size_.at(dimension)); // It will become this size after final expansions
			expansion_indices.at(dimension).resize(dim_new_size.at(dimension));
			for( uintMatrixIndex ind=dim_indices_count.at(dimension).size(); ind--; ){
				expansion_indices.at(dimension).at(ind) = ind; // The indexes we already have are staying were they are
			}

			expansion_count.at(dimension).reserve(dim_size_.at(dimension)); // It will become this size after final expansions
			expansion_count.at(dimension) = dim_indices_count.at(dimension); // Store dim_indices_count, together with the modified dim_indices_count values can be adjusted proportionally

			dim_indices_count.at(dimension).resize(dim_new_size.at(dimension), 1);

			// Calculate difference matrix
			FillDiffMatrix(diff_matrix, dimension);

			// Fill the maximum difference between indeces that are currently combined in max_diff[combined_index]
			max_diff.clear();
			max_diff.reserve(dim_new_size.at(dimension));
			max_diff.resize(expansion_count.at(dimension).size(), 0.0);

			for(uintMatrixIndex combined_index = max_diff.size(); combined_index--; ){
				max_diff.at(combined_index) = MaxDiff(diff_matrix, combined_index, same_combined_index, dim_indices_reduced.at(dimension), dim_indices_count.at(dimension));
			}

			// Do the expansion by splitting the indexes with max difference
			while(max_diff.size() < dim_new_size.at(dimension)){
				// Pick combined index to split
				uintMatrixIndex split_index = 0;
				for( auto ind = max_diff.size(); --ind; ){
					// Pick the one with the highest difference and in a tie with the highest number of combined indices
					if( max_diff.at(ind) > max_diff.at(split_index) || (max_diff.at(ind) == max_diff.at(split_index) && dim_indices_count.at(dimension).at(ind) > dim_indices_count.at(dimension).at(split_index)) ){
						split_index = ind;
					}
				}

				// Find indexes to split in the combined index
				uintMatrixIndex split_ind1(0), split_ind2(0);
				same_combined_index.clear();
				for(uintMatrixIndex current_index = dim_size_.at(dimension); current_index-- && split_ind1 == split_ind2;){
					if(split_index == dim_indices_reduced.at(dimension).at(current_index)){
						for(auto comparison_index : same_combined_index){
							// We go backwards so all comparison_index are higher than current_index which is important as diff_matrix is triangular
							if( diff_matrix.at(current_index*dim_size_.at(dimension)+comparison_index) == max_diff.at(split_index)){
								split_ind1 = current_index;
								split_ind2 = comparison_index;
								break;
							}
						}

						same_combined_index.push_back(current_index);
					}
				}

				// Split indexes
				dim_indices_reduced.at(dimension).at(split_ind2) = max_diff.size();
				--dim_indices_count.at(dimension).at(split_index);
				expansion_indices.at(dimension).at(max_diff.size()) = expansion_indices.at(dimension).at(split_index); // Point to the original index (from before the expansion)

				// Distribute indexes between the two new combined indexes minimizing difference
				double diff1, diff2;
				for(uintMatrixIndex current_index = dim_size_.at(dimension); current_index--;){
					if(split_index == dim_indices_reduced.at(dimension).at(current_index)){
						if(current_index < split_ind1){
							diff1 = diff_matrix.at(current_index*dim_size_.at(dimension)+split_ind1);
						}
						else{
							diff1 = diff_matrix.at(split_ind1*dim_size_.at(dimension)+current_index);
						}

						if(current_index < split_ind2){
							diff2 = diff_matrix.at(current_index*dim_size_.at(dimension)+split_ind2);
						}
						else{
							diff2 = diff_matrix.at(split_ind2*dim_size_.at(dimension)+current_index);
						}

						if( diff2 < diff1 ){
							dim_indices_reduced.at(dimension).at(current_index) = max_diff.size();
							--dim_indices_count.at(dimension).at(split_index);
							++dim_indices_count.at(dimension).at(max_diff.size());
						}
					}
				}

				// Update max_diff
				max_diff.at(split_index) = MaxDiff(diff_matrix, split_index, same_combined_index, dim_indices_reduced.at(dimension), dim_indices_count.at(dimension));
				max_diff.push_back( MaxDiff(diff_matrix, max_diff.size(), same_combined_index, dim_indices_reduced.at(dimension), dim_indices_count.at(dimension)) );
			}
		}
	}

	// Sum the data based on new dim_indices_reduced
	Reduce(reduced_data, dim_indices_reduced, dim_indices_count);

	return true; // We did an expansion
}

void ProbabilityEstimates::SetVariables(uintTileId num_tiles){
	error_during_fitting_ = false;
	precision_improved_ = false;

	for(auto template_segment=0; template_segment<2; ++template_segment){
		quality_.at(template_segment).clear();
		quality_.at(template_segment).resize(num_tiles);
		quality_.at(template_segment).shrink_to_fit();

		sequence_quality_.at(template_segment).clear();
		sequence_quality_.at(template_segment).resize(num_tiles);
		sequence_quality_.at(template_segment).shrink_to_fit();

		base_call_.at(template_segment).clear();
		base_call_.at(template_segment).resize(num_tiles);
		base_call_.at(template_segment).shrink_to_fit();
	}

	for( auto ref_base=dom_error_.size(); ref_base--; ){
		for( auto last_ref_base=dom_error_.at(ref_base).size(); last_ref_base--; ){
			for( auto dom_last5=dom_error_.at(ref_base).at(last_ref_base).size(); dom_last5--; ){
				dom_error_.at(ref_base).at(last_ref_base).at(dom_last5).Clear();
			}
		}
		for( auto dom_error=error_rate_.at(ref_base).size(); dom_error--; ){
			error_rate_.at(ref_base).at(dom_error).Clear();
		}
	}

	for( auto type=indels_.size(); type--; ){
		for( auto last_call=indels_.at(type).size(); last_call--; ){
			indels_.at(type).at(last_call).Clear();
		}
	}
}

void ProbabilityEstimates::IterativeProportionalFitting(
		const DataStats &stats,
		IPFDataSelector selected_data,
		uintTempSeq template_segment,
		uintTileId tile_id,
		uintBaseCall ref_base,
		uintBaseCall dom_error,
		uintBaseCall last_ref_base,
		uintNumFits max_iterations,
		double precision_aim
		){
	switch(selected_data){
	case kIPFQuality:
		if( stats.Qualities().BaseQualityForErrorRateReference(template_segment, tile_id, ref_base).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: quality, sequence quality, previous quality, position, error rate)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 10 > margins;
			auto alternative_margin = DefineMarginsQuality(stats, margins, template_segment, tile_id, ref_base);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "Quality segment " << static_cast<uintTempSeqPrint>(template_segment) << ", tile_id " << tile_id << ", ref_base " << ref_base;

			// Run the iterative proportional fitting
			quality_.at(template_segment).at(tile_id).at(ref_base).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					alternative_margin,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	case kIPFSequenceQuality:
		if( stats.Qualities().SequenceQualityMeanForGCPerTileReference(template_segment, tile_id).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: quality, previous quality, position, error rate)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > margins;
			auto alternative_margin = DefineMarginsSequenceQuality(stats, margins, template_segment, tile_id);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "Sequence quality segment " << static_cast<uintTempSeqPrint>(template_segment) << ", tile_id " << tile_id;

			// Run the iterative proportional fitting
			sequence_quality_.at(template_segment).at(tile_id).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					alternative_margin,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	case kIPFBaseCall:
		if( stats.Qualities().BaseQualityForErrorRateReference(template_segment, tile_id, ref_base, dom_error).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: called base, quality, position, error number, error rate)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 10 > margins;
			auto alternative_margin = DefineMarginsBaseCall(stats, margins, template_segment, tile_id, ref_base, dom_error);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "Base-call segment " << static_cast<uintTempSeqPrint>(template_segment) << ", tile_id " << tile_id << ", ref_base " << ref_base << ", dom_error " << dom_error;

			// Run the iterative proportional fitting
			base_call_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					alternative_margin,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	case kIPFDominantError:
		if( stats.Coverage().GCByDistance(ref_base, last_ref_base, dom_error).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: dominant error, distance, gc)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > margins;
			DefineMarginsDominantError(stats, margins, ref_base, last_ref_base, dom_error);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "Dominant-error ref_base " << ref_base << ", last_ref_base " << last_ref_base << ", dom_last5 " << dom_error;

			// Run the iterative proportional fitting
			dom_error_.at(ref_base).at(last_ref_base).at(dom_error).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					NULL,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	case kIPFErrorRate:
		if( stats.Coverage().GCByDistance(ref_base, dom_error).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: error rate, distance, gc)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > margins;
			DefineMarginsErrorRate(stats, margins, ref_base, dom_error);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "Error-rate ref_base " << ref_base << ", dom_error " << dom_error;

			// Run the iterative proportional fitting
			error_rate_.at(ref_base).at(dom_error).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					NULL,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	case kIPFInDels:
		if( stats.Errors().InDelByInDelPos(template_segment, last_ref_base).size() ){
			// Something has to be done as data is not empty
			// Start with defining the margins (Dimension order: error rate, distance, gc)
			array< pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > margins;
			DefineMarginsIndels(stats, margins, template_segment, last_ref_base);

			// Write the description for printing
			stringstream descriptor;
			descriptor << "InDel type " << static_cast<uintInDelTypePrint>(template_segment) << ", last_call " << last_ref_base;

			// Run the iterative proportional fitting
			indels_.at(template_segment).at(last_ref_base).IterativeProportionalFitting(
					error_during_fitting_,
					precision_improved_,
					precision_aim,
					max_iterations,
					margins,
					NULL,
					descriptor.str(),
					print_mutex_
					);
			return;
		}
		break;
	}
}

void ProbabilityEstimates::IPFThread(ProbabilityEstimates &self, const DataStats &stats, const std::vector<IPFThreadParams> &params, uintNumFits max_iterations, double precision_aim){
	decltype(params.size()) cur_par(self.current_param_++);

	for(; cur_par < params.size() && !self.error_during_fitting_; cur_par = self.current_param_++){
		self.IterativeProportionalFitting(stats, params.at(cur_par).selected_data, params.at(cur_par).template_segment, params.at(cur_par).tile_id, params.at(cur_par).ref_base, params.at(cur_par).dom_error, params.at(cur_par).last_ref_base, max_iterations, precision_aim);
	}
}

void ProbabilityEstimates::PrepareResult(){
	for( auto template_segment=2; template_segment--; ){
		quality_result_.at(template_segment).resize(quality_.at(template_segment).size());
		for( auto tile_id=quality_.at(template_segment).size(); tile_id--; ){
			for( auto ref_base=4; ref_base--; ){
				quality_.at(template_segment).at(tile_id).at(ref_base).FullExpansion();
				quality_result_.at(template_segment).at(tile_id).at(ref_base).GetResults(quality_.at(template_segment).at(tile_id).at(ref_base).estimates_, quality_.at(template_segment).at(tile_id).at(ref_base).dim_indices_);
				quality_result_.at(template_segment).at(tile_id).at(ref_base).ImputeMissingValues();
			}
		}
		quality_.at(template_segment).clear(); // Free space by removing the calculation storage for the ipf
		quality_.at(template_segment).shrink_to_fit();

		sequence_quality_result_.at(template_segment).resize(sequence_quality_.at(template_segment).size());
		for( auto tile_id=sequence_quality_.at(template_segment).size(); tile_id--; ){
			sequence_quality_.at(template_segment).at(tile_id).FullExpansion();
			sequence_quality_result_.at(template_segment).at(tile_id).GetResults(sequence_quality_.at(template_segment).at(tile_id).estimates_, sequence_quality_.at(template_segment).at(tile_id).dim_indices_);
			sequence_quality_result_.at(template_segment).at(tile_id).ImputeMissingValues();
		}
		sequence_quality_.at(template_segment).clear(); // Free space by removing the calculation storage for the ipf
		sequence_quality_.at(template_segment).shrink_to_fit();

		base_call_result_.at(template_segment).resize(base_call_.at(template_segment).size());
		for( auto tile_id=base_call_.at(template_segment).size(); tile_id--; ){
			for( auto ref_base=4; ref_base--; ){
				for( auto dom_error=5; dom_error--; ){
					base_call_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).FullExpansion();
					base_call_result_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).GetResults(base_call_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).estimates_, base_call_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).dim_indices_);
					base_call_result_.at(template_segment).at(tile_id).at(ref_base).at(dom_error).ImputeMissingValues();
				}
			}
		}
		base_call_.at(template_segment).clear(); // Free space by removing the calculation storage for the ipf
		base_call_.at(template_segment).shrink_to_fit();
	}

	for( auto ref_base=4; ref_base--; ){
		for( auto last_base=5; last_base--; ){
			for( auto dom_last5=5; dom_last5--; ){
				dom_error_.at(ref_base).at(last_base).at(dom_last5).FullExpansion();
				dom_error_result_.at(ref_base).at(last_base).at(dom_last5).GetResults(dom_error_.at(ref_base).at(last_base).at(dom_last5).estimates_, dom_error_.at(ref_base).at(last_base).at(dom_last5).dim_indices_);
				dom_error_result_.at(ref_base).at(last_base).at(dom_last5).ImputeMissingValues();
				dom_error_.at(ref_base).at(last_base).at(dom_last5).Clear(); // Free space by removing the calculation storage for the ipf
			}
			error_rate_.at(ref_base).at(last_base).FullExpansion();
			error_rate_result_.at(ref_base).at(last_base).GetResults(error_rate_.at(ref_base).at(last_base).estimates_, error_rate_.at(ref_base).at(last_base).dim_indices_);
			error_rate_result_.at(ref_base).at(last_base).ImputeMissingValues();
			error_rate_.at(ref_base).at(last_base).Clear(); // Free space by removing the calculation storage for the ipf
		}
	}

	for( auto type=2; type--; ){
		for( auto last_call=6; last_call--; ){
			indels_.at(type).at(last_call).FullExpansion();
			indels_result_.at(type).at(last_call).GetResults(indels_.at(type).at(last_call).estimates_, indels_.at(type).at(last_call).dim_indices_);
			indels_result_.at(type).at(last_call).ImputeMissingValues();
			indels_.at(type).at(last_call).Clear(); // Free space by removing the calculation storage for the ipf
		}
	}
}

bool ProbabilityEstimates::Load( const char *archive_file ){
	if( !FileExists(archive_file) ){
		printErr << "File '" << archive_file << "' does not exists or no read permission given." << std::endl;
		return false;
	}

	try{
		// create and open an archive for input
		ifstream ifs(archive_file);
		boost::archive::text_iarchive ia(ifs);

		// read class state from archive
		ia >> *this;
	}
	catch(const exception& e){
		printErr << "Could not load probability estimates from '" << archive_file << "': " << e.what() << std::endl;
		return false;
	}

	error_during_fitting_ = false;
	precision_improved_ = false;

	return true;
}

bool ProbabilityEstimates::Save( const char *archive_file ) const{
	try{
		// create all missing directories
		CreateDir(archive_file);

		// create and open a character archive for output
		ofstream ofs(archive_file);
		boost::archive::text_oarchive oa(ofs);

		// save data to archive
		oa << *this;
	}
	catch(const exception& e){
		printErr<< "Could not save probability estimates to '" << archive_file << "': " << e.what() << std::endl;
		return false;
	}

	return true;
}

bool ProbabilityEstimates::Estimate(const DataStats &stats, uintNumFits max_iterations, double precision_aim, uintNumThreads num_threads, const char *output, const char *input){
	if( string("") == input){
		printInfo << "Starting new probability estimates" << std::endl;
		stats_creation_time_ = stats.CreationTime();
		SetVariables(stats.Tiles().NumTiles());
	}
	else{
		printInfo << "Loading probability estimates from '" << input << "'" << std::endl;
		if( !this->Load(input) ){
			return false;
		}

		if( stats_creation_time_ != stats.CreationTime() ){
			printInfo << "Statistics file has been updated, the probability estimates will be overwritten." << std::endl;
			stats_creation_time_ = stats.CreationTime();
			SetVariables(stats.Tiles().NumTiles());
		}
	}

	printInfo << "Aiming for a precision lower than " << precision_aim << "%" << std::endl;
	printInfo << "Using a maximum of " << max_iterations << " iterations" << std::endl;

	// Collect the different parameters determining the different matrices that have to be calculated
	vector<IPFThreadParams> params;
	params.reserve( (2*4 + 2*4*5 + 2)*stats.Tiles().NumTiles() + 4*4 + 4*4*5 + 2*6 );
	for(uintTempSeq template_segment=0; template_segment<2; ++template_segment){
		for(uintTileId tile_id=0; tile_id<stats.Tiles().NumTiles(); ++tile_id){
			for( uintBaseCall ref_base = base_call_.at(template_segment).at(tile_id).size(); ref_base--; ){
				params.push_back( {kIPFQuality,template_segment,tile_id,ref_base,0,0} );
			}
		}
	}
	for(uintTempSeq template_segment=0; template_segment<2; ++template_segment){
		for(uintTileId tile_id=0; tile_id<stats.Tiles().NumTiles(); ++tile_id){
			for( uintBaseCall ref_base = base_call_.at(template_segment).at(tile_id).size(); ref_base--; ){
				for( uintBaseCall dom_error = base_call_.at(template_segment).at(tile_id).at(ref_base).size(); dom_error--; ){
					params.push_back( {kIPFBaseCall,template_segment,tile_id,ref_base,dom_error,0} );
				}
			}
		}
	}

	for(uintTempSeq template_segment=0; template_segment<2; ++template_segment){
		for(uintTileId tile_id=0; tile_id<stats.Tiles().NumTiles(); ++tile_id){
			params.push_back( {kIPFSequenceQuality,template_segment,tile_id,0,0,0} );
		}
	}

	for( uintBaseCall ref_base = 4; ref_base--; ){
		for( uintBaseCall dom_error = 5; dom_error--; ){
			if(ref_base != dom_error){
				params.push_back( {kIPFErrorRate,0,0,ref_base,dom_error,0} );
			}
		}
		for( uintBaseCall dom_base = 4; dom_base--; ){
			for( uintBaseCall last_ref_base = 5; last_ref_base--; ){
				params.push_back( {kIPFDominantError,0,0,ref_base,dom_base,last_ref_base} );
			}
		}
	}

	std::array< std::pair<const Vect<Vect<uintMatrixCount>> *, bool>, 6 > indel_margins;
	for(uintInDelType type=0; type<2; ++type){
		for( uintBaseCall last_call=6; last_call--; ){
			params.push_back( {kIPFInDels,type,0,0,0,last_call} );
		}
	}
	current_param_ = 0;

	// Run the predetermined parameters in defined number of threads
	if( num_threads > params.size() ){
		num_threads = params.size();
	}
	thread t[num_threads];
	for(auto i = num_threads; i--; ){
		t[i] = std::thread(IPFThread, std::ref(*this), std::cref(stats), std::ref(params), max_iterations, precision_aim/100); // Change precision aim from % into factor
	}
	for(auto i = num_threads; i--; ){
		t[i].join();
	}

	if( error_during_fitting_ ){
		DeleteFile( output ); // Remove the probabilities file if one existed previously so it is clear that we encountered an error and do not accidentally continue with the old file
		return false;
	}

	if(max_iterations){
		if( precision_improved_ ){
			if( string("") == input){
				printInfo << "Writing probability estimates to " << output << std::endl;
			}
			else{
				printInfo << "Updating probability estimates in " << output << std::endl;
			}

			if( !this->Save(output) ){
				return false;
			}
		}
		else{
			printInfo << "Precision aim already met by loaded probabilities. Probabilities won't be updated." << std::endl;
		}
	}
	else{
		printInfo << "No iterations allowed. Probabilities could not be improved." << std::endl;
	}

	return true;
}
