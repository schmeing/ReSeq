#ifndef VECT_HPP
#define VECT_HPP

#include <algorithm>
#include <limits>
#include <ostream>
#include <stdint.h>
#include <vector>

#include "gtest/gtest.h"

#include <boost/serialization/vector.hpp>

#include "reportingUtils.hpp"

#include "utilities.h"

namespace reseq{
	template<typename T> class Vect{ // std::vector that removes 0s at the beginning and has an offset variable instead and automatically increases size if called outside range
	private:
		std::pair< typename std::vector<T>::size_type, std::vector<T> > vec_; // Offset + Vector

		static const T dummy_; // In case a reference to an empty element has to be returned

		inline void Ensure( typename std::vector<T>::size_type n ){ // Ensure that element with id n exists
			// Make sure offset is not to high
			if( vec_.first > n ){
				vec_.second.insert( vec_.second.begin(), vec_.first-n, 0 );
				vec_.first = n;
			}

			// Make sure size of vector is long enough
			if( vec_.second.size() <= n-vec_.first ){
				vec_.second.resize( n-vec_.first+1, 0 );
			}
		}

		// Boost serialization
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			 ar & vec_;
		}

		// Google test
		FRIEND_TEST(VectTest, BasicFunctionality);
		FRIEND_TEST(VectTest, CopyAndClear);

	public:
		// Constructors
		Vect(){}
		Vect( typename std::vector<T>::size_type reserved ){
			this->reserve(reserved);
		}
		Vect( const typename std::vector<T>::size_type offset, const typename std::vector<T> &values ):vec_{offset, values}{}
		Vect(const Vect& x):vec_(x.vec_){}
		Vect& operator=(const Vect& x){ vec_ = x.vec_; }
		Vect& operator=(const std::vector<T>& x){
			vec_.first = 0;
			vec_.second = x;
		}
		template<typename U> Vect& operator=(const std::vector<utilities::VectorAtomic<U>>& x){
			Clear();

			// Find first element and set offset
			typename std::vector<U>::size_type first(0);
			while(first < x.size() && 0 == x[first]){
				++first;
			}

			if(first < x.size()){
				vec_.first = first;

				// Find last element and resize vector
				auto last = x.size();
				while( 0 == x[--last] );

				vec_.second.resize(last+1-first);
				for(auto i=first; i <= last; ++i){
					vec_.second[i-first] = x[i];
				}
			}
		}
		template<typename U> Vect& operator=(const std::vector<U>& x){
			Clear();

			// Find first element and set offset
			typename std::vector<U>::size_type first(0);
			while(first < x.size() && 0 == x[first].size()){
				++first;
			}

			if(first < x.size()){
				vec_.first = first;

				// Find last element and resize vector
				auto last = x.size();
				while( 0 == x[--last].size() );

				vec_.second.resize(last+1-first);
				for(auto i=first; i <= last; ++i){
					vec_.second[i-first] = x[i];
				}
			}
		}

		// Iterators
		typename std::vector<T>::iterator begin(){ return vec_.second.begin(); }
		typename std::vector<T>::const_iterator begin() const{ return vec_.second.begin(); }
		typename std::vector<T>::iterator end(){ return vec_.second.end(); }
		typename std::vector<T>::const_iterator end() const{ return vec_.second.end(); }

		// Capacity
		typename std::vector<T>::size_type size() const{ return vec_.second.size(); }
		typename std::vector<T>::size_type max_size() const{ return vec_.second.max_size(); }
		typename std::vector<T>::size_type capacity() const{ return vec_.second.capacity(); }
		bool empty() const{ return vec_.second.empty(); }
		void reserve(typename std::vector<T>::size_type n){ vec_.second.reserve(n); }

		typename std::vector<T>::size_type from() const{ return vec_.first; } // Returns the id of the first element
		typename std::vector<T>::size_type to() const{ return vec_.first + vec_.second.size(); } // Returns the id of the last element + 1
		void Shift( int64_t shift_by ){ // Shift all ids by shift_by
			if(-shift_by > static_cast<int64_t>(vec_.first)){ // In case we would get below 0 offset
				vec_.first = 0;
			}
			else{
				vec_.first += shift_by;
			}
		}
		void SetOffset(typename std::vector<T>::size_type offset){ // Sets the offset to the given value if that doesn't conflict with values already in the array. In case of a conflict it sets the offset as high as possible
			if( vec_.second.empty() ){
				vec_.first = offset;
			}
			else{
				if( offset > vec_.first){
					typename std::vector<T>::size_type leading_zeros = 0;
					while( leading_zeros < vec_.second.size() && 0 == vec_.second[leading_zeros] && leading_zeros < offset-vec_.first){ // Check if all entries from current offset to new offset are zero and if not how many leading zeros can be removed
						++leading_zeros;
					}
					vec_.second.erase( vec_.second.begin(), vec_.second.begin()+leading_zeros );
					vec_.first += leading_zeros;
				}
				else{
					vec_.second.insert( vec_.second.begin(), vec_.first-offset, 0 );
					vec_.first = offset;
				}
			}
		}
		void SoftShrink(){ // Shrink vector to minimum size: Remove leading and trailing zeros without freeing reserved space
			// Remove trailing zeros
			typename std::vector<T>::size_type last_non_zero = vec_.second.size();
			while( last_non_zero-- && 0 == vec_.second[last_non_zero] ); // Beware that 0 <= last_non_zero does not work as size_type is most likely unsigned
			vec_.second.resize(last_non_zero+1);

			// Remove leading zeros
			typename std::vector<T>::size_type leading_zeros = 0;
			while( leading_zeros < vec_.second.size() && 0 == vec_.second[leading_zeros] ){
				++leading_zeros;
			}
			vec_.second.erase( vec_.second.begin(), vec_.second.begin()+leading_zeros );
			vec_.first += leading_zeros;

		}
		void Shrink(){ // Shrink vector to minimum size: Remove leading and trailing zeros + clear all reserved space that is not used
			SoftShrink();

			// Clear all reserved space that is not used
			vec_.second.shrink_to_fit();
		}

		// Element access
		T &front(){ return vec_.second.front(); }
		const T &front() const{ return vec_.second.front(); }
		T &back(){ return vec_.second.back(); }
		const T &back() const{ return vec_.second.back(); }

		T &operator[](typename std::vector<T>::size_type n){ // Ensure that the called element exists (if not create it) and return a reference to it
			Ensure(n);
			return vec_.second[ n-vec_.first ]; // Subtract offset
		}
		const T &operator[](typename std::vector<T>::size_type n) const{
			if( this->from() <= n && n < this->to() ){
				return vec_.second[ n-vec_.first ]; // Subtract offset
			}
			else{
				return dummy_;
			}
		}
		T &at(typename std::vector<T>::size_type n){
			if( this->from() > n ){
				throw std::out_of_range( "Index lower than threshold called" );
			}
			else if(this->to() <= n){
				throw std::out_of_range( "Index higher than last element called" );
			}
			else{
				return vec_.second[ n-vec_.first ]; // Subtract offset
			}
		}
		const T &at(typename std::vector<T>::size_type n) const{
			if( this->from() > n ){
				throw std::out_of_range( "Index lower than threshold called" );
			}
			else if(this->to() <= n){
				throw std::out_of_range( "Index higher than last element called" );
			}
			else{
				return vec_.second[ n-vec_.first ]; // Subtract offset
			}
		}
		const std::pair< typename std::vector<T>::size_type, typename std::vector<T> > &std() const{ return vec_; }

		const T &Max() const{
			auto max_it = std::max_element(vec_.second.begin(), vec_.second.end());
			return *max_it;
		}

		// Modifiers
		void PushBack(const T& val){ vec_.second.push_back(val); }
		void PopBack(){ vec_.second.pop_back(); }
		void PopBack(typename std::vector<T>::size_type num_elements){ vec_.second.resize( vec_.second.size() - num_elements ); }
		void Set( typename std::vector<T>::size_type from, typename std::vector<T>::size_type to, const T& value ){
			// Count backwards so the full memory is allocated directly as we start with the highest index
			for( auto i = to; from < i--; ){ //-- included in second term so it also works for 0==from given unsigned integers
				(*this)[i] = value;
			}
		}

		void Clear(){ // Clears the vector and resets the offset
			vec_.second.clear();
			vec_.first = 0;
		}

		template<typename U>  void Acquire(std::vector<U>& x){
			*this = x;
			x.clear();
			x.shrink_to_fit();
		}

		//Operators
		bool operator==(const Vect<T>& right) const{
			decltype(this->from()) from;
			if( this->from() > right.from()){
				from = right.from();
			}
			else{
				from = this->from();
			}

			decltype(this->to()) to;
			if( this->to() > right.to()){
				to = right.to();
			}
			else{
				to = this->to();
			}

			for( auto i=from; i<to; ++i ){
				if( !( (*this)[i] == right[i]) ){
					return false;
				}
			}

			return true;
		}
		Vect<T>& operator/=(T right){
			for( auto i=vec_.second.size(); i--; ){
				vec_.second[i] /= right;
			}
			return *this;
		}
		Vect<T>& operator+=(const Vect<T>& right){
			for( auto i=right.to(); right.from() < i--; ){
				(*this)[i] += right[i];
			}
			return *this;
		}
	};

	template<typename T> bool operator==( typename std::vector<T>::size_type val, const Vect<T>& vect){ // Compares the size_t value to the size
		return val == vect.size();
	}

	template<typename T> std::ostream& operator<<(std::ostream& os, const Vect<T>& vect){
		for( auto x=vect.from(); x < vect.to(); ++x ){
			os << x << ": " << vect[x] << '\n';
		}
		return os;
	}

	template<typename T> const T Vect<T>::dummy_ = 0;

	// GetLimitsSecondDimension
	template<template<typename> class T, typename U> void GetLimitsSecondDimension( const T<U> &vect2d, typename std::vector<U>::size_type &from, typename std::vector<U>::size_type &to ){
		from = std::numeric_limits< typename std::vector<U>::size_type >::max();
		to = 0;
		for( auto &vectLowerDim : vect2d){
			if(vectLowerDim.size()){
				reseq::utilities::SetToMin(from, vectLowerDim.from());
				reseq::utilities::SetToMax(to, vectLowerDim.to());
			}
		}
	}

	template<template<typename> class T, typename U> void GetIndicesFirstDimension( const T<U> &vect2d, std::vector<uint64_t> &indices ){
		indices.clear();
		indices.reserve(vect2d.size());
		for( auto i=vect2d.from(); i < vect2d.to(); ++i){
			if( vect2d[i].size() ){
				indices.push_back(i);
			}
		}
	}

	template<template<typename> class T, typename U> void GetIndicesSecondDimension( const T<U> &vect2d, std::vector<uint64_t> &indices ){
		// Count the indices from the second dimension in a 1d vect
		Vect<uint64_t> index_count;

		for( auto &vectLowerDim : vect2d ){
			for( auto i=vectLowerDim.to(); i-- > vectLowerDim.from(); ){
				if( vectLowerDim[i] ){
					++index_count[i];
				}
			}
		}

		// Use 1d vect to find the existing indices
		indices.clear();
		indices.reserve(index_count.size());
		for( auto i=index_count.from(); i < index_count.to(); ++i){
			if( index_count[i] ){
				indices.push_back(i);
			}
		}
	}

	// ShiftVect
	template<typename T> void Shift2d( T &vect2d, int64_t outer_shift, int64_t inner_shift ){
		for( auto &vect1d: vect2d){
			vect1d.Shift(inner_shift);
		}
		vect2d.Shift(outer_shift);
	}

	// ShrinkVect
	template<template<typename> class T> inline void ShrinkVect( T<signed char> &vect ){
		vect.Shrink();
	}
	template<template<typename> class T> inline void ShrinkVect( T<unsigned char> &vect ){
		vect.Shrink();
	}
	template<template<typename> class T> inline void ShrinkVect( T<uint16_t> &vect ){
		vect.Shrink();
	}
	template<template<typename> class T> inline void ShrinkVect( T<uint64_t> &vect ){
		vect.Shrink();
	}
	template<template<typename> class T> inline void ShrinkVect( T<double> &vect ){
		vect.Shrink();
	}
	template<typename T> void ShrinkVect( T &vect ){
		for( auto &vectLowerDim: vect){
			ShrinkVect(vectLowerDim);
		}
		vect.Shrink();
	}

	// SetVect
	template<template<typename> class T, template<typename> class U, typename V> void SetVect2d(
			T<U<V>> &vect,
			typename std::vector<U<V>>::size_type outer_from,
			typename std::vector<U<V>>::size_type outer_to,
			typename std::vector<V>::size_type inner_from,
			typename std::vector<V>::size_type inner_to,
			V value ){
		for( auto i = outer_to; outer_from < i--; ){
			vect[i].Set(inner_from, inner_to, value);
		}
	}

	// ClearSecondDimension
	template<typename T> void ClearSecondDimension( T &vect2d ){
		for( auto &vectLowerDim : vect2d){
			vectLowerDim.Clear();
		}
	}

	// SumVect
	template<typename T> uint64_t SumVectHelper( const T &vect ){
		uint64_t sum = 0;
		for( auto val: vect){
			sum += val;
		}
		return sum;
	}

	template<typename T> uint64_t SumVectHelperAtomic( const T &vect ){
		uint64_t sum = 0;
		for( const auto &val: vect){
			sum += val;
		}
		return sum;
	}

	template<typename T> double SumVectHelperDouble( const T &vect ){
		double sum = 0;
		for( auto val: vect){
			sum += val;
		}
		return sum;
	}

	template<template<typename> class T> inline uint64_t SumVect( const T<uint16_t> &vect ){
		return SumVectHelper(vect);
	}

	template<template<typename> class T> inline uint64_t SumVect( const T<uint64_t> &vect ){
		return SumVectHelper(vect);
	}

	template<template<typename> class T> inline uint64_t SumVect( const T<utilities::VectorAtomic<uint64_t>> &vect ){
		return SumVectHelperAtomic(vect);
	}

	template<template<typename> class T> inline double SumVectD( const T<double> &vect ){
		return SumVectHelperDouble(vect);
	}

	template<template<class, class> class T> inline uint64_t SumVect( const T<uint64_t, std::allocator<uint64_t>> &vect ){ // For std containers
		return SumVectHelper(vect);
	}

	template<template<class, class> class T> inline uint64_t SumVect( const T<utilities::VectorAtomic<uint64_t>, std::allocator<utilities::VectorAtomic<uint64_t>>> &vect ){ // For std containers
		return SumVectHelperAtomic(vect);
	}

	template<template<class, class> class T> inline double SumVectD( const T<double, std::allocator<double>> &vect ){ // For std containers
		return SumVectHelperDouble(vect);
	}

	template<std::size_t N> uint64_t SumVect( const std::array<uint64_t, N> &vect ){
		uint64_t sum = 0;
		for( auto i=vect.size(); i--; ){
			sum += vect[i];
		}
		return sum;
	}

	template<std::size_t N> double SumVectD( const std::array<double, N> &vect ){
		double sum = 0;
		for( auto i=vect.size(); i--; ){
			sum += vect[i];
		}
		return sum;
	}

	template<std::size_t N> uint64_t SumVect( const std::array<utilities::VectorAtomic<uint64_t>, N> &vect ){
		uint64_t sum = 0;
		for( auto i=vect.size(); i--; ){
			sum += vect[i];
		}
		return sum;
	}

	template<typename T> uint64_t SumVect( const T &vect ){
		uint64_t sum = 0;
		for( auto &vectLowerDim: vect){
			sum += SumVect(vectLowerDim);
		}
		return sum;
	}

	template<typename T, std::size_t N> uint64_t SumVect( const std::array<T, N> &vect ){
		uint64_t sum = 0;
		for( auto i=vect.size(); i--; ){
			sum += SumVect(vect[i]);
		}
		return sum;
	}

	template<typename T, std::size_t N> double SumVectD( const std::array<T, N> &vect ){
		double sum = 0;
		for( auto i=vect.size(); i--; ){
			sum += SumVect(vect[i]);
		}
		return sum;
	}

	template<typename T> double SumVectD( const T &vect ){
		double sum = 0;
		for( auto &vectLowerDim: vect){
			sum += SumVectD(vectLowerDim);
		}
		return sum;
	}
}

#endif // VECT_HPP
