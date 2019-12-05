#ifndef SEQQUALITYSTATS_HPP
#define SEQQUALITYSTATS_HPP

#include <cmath>
#include <stdexcept>
#include <stdint.h>

#include "reportingUtils.hpp"

#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{

	template<typename T> class SeqQualityStats{
	private:
		// Boost serialization
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			 ar & qualities_;
		}
	public:
		Vect<T> qualities_;

		unsigned char mean_;
		unsigned char minimum_;
		unsigned char first_quartile_; // Lower median of lower median-excluded(only for uneven numbers) half
		unsigned char median_; // Upper median
		unsigned char third_quartile_; // Upper median of upper median-excluded(only for uneven numbers) half
		unsigned char maximum_;

		unsigned char probability_mean_; // Mean of the error probabilities expressed as quality

		SeqQualityStats():
			mean_(0),
			minimum_(0),
			first_quartile_(0),
			median_(0),
			third_quartile_(0),
			maximum_(0),
			probability_mean_(0){}
		SeqQualityStats( decltype( qualities_.capacity() ) reserved ):
			SeqQualityStats(){
			qualities_.reserve(reserved);
		}
		template<typename U> SeqQualityStats<T>& operator=(const std::vector<U>& x){
			qualities_ = x;
			return *this;
		}

		// Forward functions to the Vect qualities_
		typename std::vector<T>::iterator begin(){ return qualities_.begin(); }
		typename std::vector<T>::const_iterator begin() const{ return qualities_.begin(); }
		typename std::vector<T>::iterator end(){ return qualities_.end(); }
		typename std::vector<T>::const_iterator end() const{ return qualities_.end(); }

		decltype( qualities_.size() ) size() const{ return qualities_.size(); }
		decltype( qualities_.max_size() ) max_size() const{ return qualities_.max_size(); }
		decltype( qualities_.capacity() ) capacity() const{ return qualities_.capacity(); }
		decltype( qualities_.empty() ) empty() const{ return qualities_.empty(); }
		void reserve(decltype( qualities_.size() ) n){ qualities_.reserve(n); }

		decltype( qualities_.from() ) from() const{ return qualities_.from(); }
		decltype( qualities_.to() ) to() const{ return qualities_.to(); }
		void Shift(decltype( qualities_.size() ) shift_by ){ qualities_.Shift(shift_by); }
		void SetOffset(decltype( qualities_.size() ) offset){ qualities_.SetOffset(offset); }
		void Shrink(){ qualities_.Shrink(); }

		T &front(){ return qualities_.front(); }
		const T &front() const{ return qualities_.front(); }
		T &back(){ return qualities_.back(); }
		const T &back() const{ return qualities_.back(); }
		T &operator[](decltype( qualities_.size() ) n){ return qualities_[n]; }
		const T &operator[](decltype( qualities_.size() ) n) const{ return qualities_[n]; }
		T &at(decltype( qualities_.size() ) n){ return qualities_.at(n); }
		const T &at(decltype( qualities_.size() ) n) const{ return qualities_.at(n); }
		const std::pair< typename std::vector<T>::size_type, typename std::vector<T> > &stdQualities() const{ return qualities_.std(); }

		void PushBack(const T& val){ qualities_.PushBack(val); }
		void PopBack(){ qualities_.PopBack(); }
		template<typename U>  void Acquire(std::vector<U>& x){ qualities_.Acquire(x); }

		SeqQualityStats<T>& operator+=(const SeqQualityStats<T>& right){
			this->qualities_ += right.qualities_;
			return *this;
		}
		template<typename U> SeqQualityStats<T>& operator+=(const std::vector<U>& right){
			this->qualities_ += right;
			return *this;
		}

		// Own class function
		void Clear(){ // Clears qualities_ and the statistics
			qualities_.Clear();
			mean_ = 0;
			minimum_ = 0;
			first_quartile_ = 0;
			median_ = 0;
			third_quartile_ = 0;
			maximum_ = 0;
			probability_mean_= 0;
		}

		void Calculate(bool all=true){
			uint64_t num_values(0);
			uint64_t sum(0);

			for(auto qual = qualities_.from(); qual < qualities_.to(); ++qual){
				num_values += qualities_.at(qual);
				sum += qualities_.at(qual)*qual;
			}

			if( 0 == num_values ){
				mean_ = 0;
				minimum_ = 0;
				first_quartile_ = 0;
				median_ = 0;
				third_quartile_ = 0;
				maximum_ = 0;
			}
			else if( 1 == num_values ){ // Quartiles not defined as median-excluded list is empty: Set them to the only existing element
				mean_ = sum;
				minimum_ = sum;
				first_quartile_ = sum;
				median_ = sum;
				third_quartile_ = sum;
				maximum_ = sum;
			}
			else{
				mean_ = static_cast<unsigned char>((sum+num_values/2)/num_values); // rounded normally

				if(all){
					uint64_t n_first_quartile( (num_values+2)/4 ); // Lower number in case of hitting the middle
					uint64_t n_median( num_values/2 + 1 ); // Upper number in case of hitting the middle
					uint64_t n_third_quartile( (num_values*3+5)/4 ); // Upper number in case of hitting the middle

					uint64_t num_values2( num_values );
					auto qual = qualities_.to();
					while( 0 == qualities_.at(--qual) ); // Ignore all zeros at the end, we don't need to check here for [qual] existing because we are guaranteed to have two non-zero elements in the vector
					maximum_ = static_cast<unsigned char>(qual); // It is not qualities_.to()-1, because we do not shrink before and therefore might have empty bins at the top
					while( 0 < num_values2 ){
						minimum_ = qual;

						if(num_values2 >= n_first_quartile){
							first_quartile_ = qual;

							if(num_values2 >= n_median){
								median_ = qual;

								if(num_values2 >= n_third_quartile){
									third_quartile_ = qual;
								}
							}
						}

						// Subtract after the comparison, because it's going backwards
						num_values2 -= qualities_.at(qual--);
					}
				}
			}
		}

		void CalculateProbabilityMean(){
			uint64_t num_values(0);
			double sum(0.0);

			for(auto qual = qualities_.from(); qual < qualities_.to(); ++qual){
				num_values += qualities_.at(qual);
				sum += qualities_.at(qual) * std::pow( 10.0, -static_cast<double>(qual)/10.0 );
			}

			if( 0 == num_values ){
				probability_mean_ = 0;
			}
			else{
				probability_mean_ = static_cast<unsigned char>( std::round( -10.0*std::log10(sum/num_values) ) );
			}
		}
	};

	template<typename T> bool operator==( decltype( SeqQualityStats<T>().size() ) val, const SeqQualityStats<T>& stats){ // Compares the size_t value to the quality vector
		return val == stats.qualities_;
	}

	template<typename T> Vect<T>& operator+=(Vect<T>& left, const SeqQualityStats<T>& right){
		for( auto i=right.to(); right.from() < i--; ){
			left[i] += right.at(i);
		}
		return left;
	}
}
	
#endif // SEQQUALITYSTATS_HPP
