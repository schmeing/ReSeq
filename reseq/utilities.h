#ifndef UTILITIES_H
#define UTILITIES_H

#include <atomic>
#include <random>
#include <stdexcept>

#include <boost/filesystem.hpp>

#include "reportingUtils.hpp"

#include <seqan/sequence.h>

//https://stackoverflow.com/questions/3599160/how-to-suppress-unused-parameter-warnings-in-c
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

namespace reseq{
	namespace utilities{
		inline void CreateDir(const char *file){
			auto file_path = boost::filesystem::path(file).parent_path();
			if( !file_path.empty() ){
				boost::filesystem::create_directories(file_path);
			}
		}

		inline uint64_t Divide(uint64_t nom, uint64_t den){ // Calculates the division: nom/den with proper rounding to int
			return (nom+den/2)/den;
		}

		inline unsigned char Percent(uint64_t nom, uint64_t den){ // Calculates the percentage: nom*100/den with proper rounding to int
			return static_cast<unsigned char>( Divide(nom*100,den) );
		}

		inline unsigned char SafePercent(uint64_t nom, uint64_t den){
			if(den){
				return Percent(nom,den);
			}
			else{
				return 50;
			}
		}

		inline uint64_t MeanWithRoundingToFirst(uint32_t first, uint32_t second){ // Mean between first and second and rounding in the direction of first
			return (first+second+(first>second?1:0))/2;
		}

		template<typename T, typename U> inline void SetToMin(T& current_min, const U potential_min){
			if(potential_min < current_min ){
				current_min = potential_min;
			}
		}

		template<typename T, typename U> inline void SetToMax(T& current_max, const U potential_max){
			if(potential_max > current_max ){
				current_max = potential_max;
			}
		}

		template <typename T> inline int Sign(T val) {
		    return (T(0) < val) - (val < T(0));
		}

		inline double InvLogit2(const double bias){
			return 2/(1+exp(-bias));
		}

		template<typename T, uint64_t N, typename S> void GetDominantLastX( T &dom_base, std::array<uint64_t,N> &seq_content, uint16_t lastx, const S &seq, uint32_t cur_pos ){
			if(cur_pos){
				// Fill seq_content_reference_last5
				for( uint32_t pos = (lastx<cur_pos?cur_pos-lastx:0); pos < cur_pos; ++pos){
					++seq_content.at( seq[pos] );
				}

				// Get count of most appearing base
				uint32_t max_content(0);
				for( uint16_t nucleotide=4; nucleotide--; ){ // Leave out N for this
					SetToMax(max_content, seq_content.at(nucleotide));
				}

				if(0 == max_content){
					// Only N's in lastx
					dom_base = 4;
				}
				else{
					// Get the base that is closest to the current base which matches max_content (in case there are multiple bases with same number of appearances)
					uint32_t pos=cur_pos;
					while(max_content != seq_content.at( seq[--pos] ));
					dom_base = seq[pos];
				}
			}
		}

		template<typename T, uint64_t N, typename U, typename S> void SetDominantLastX( T &dom_base, std::array<uint64_t,N> &seq_content, U base, uint16_t lastx, const S &seq, uint32_t pos ){
			// Add new base
			++seq_content.at(base);
			// Remove old base if we reached lastx bases, so that the sum of seq_content is actually lastx
			if(lastx <= pos) --seq_content.at(static_cast<seqan::Dna5>(seq[pos-lastx]));
			// If we added an N, we just ignore it
			if(4 > base){
				// In case we don't have a defined base yet or the just added base has equal or more appearances than the current dominant one it's easy
				if(dom_base >= N || seq_content.at(base) >= seq_content.at(dom_base)){
					dom_base = base;
				}
				else if(lastx <= pos && dom_base == seq[pos-lastx]){ // In case seq content is not complete yet only the easy case above can happen
					// Otherwise check if there are two most abundant bases
					bool equal_content(false);
					for( uint16_t b=N; b--; ){
						if( b != dom_base && seq_content.at(b) >= seq_content.at(dom_base)){
							equal_content = true;
						}
					}
					// If two bases are most abundant take the one closer to pos
					if(equal_content){
						for( uint16_t p=1; p<lastx; ++p){
							if( seq_content.at(seq[pos-p]) >= seq_content.at(dom_base) ){
								dom_base = seq[pos-p];
								break;
							}
						}
					}
				}
			}
		}

		template<typename T> struct VectorAtomic{
			typename std::atomic<T> value_;

			VectorAtomic():
				value_(0){
			}

			VectorAtomic(const VectorAtomic &UNUSED(right)){
				throw std::runtime_error( "This function should never be called. Did you resize an object with a non-zero length?" );
			}

			operator T() const{ return value_; }

			template<typename U> inline VectorAtomic<T> &operator=(const U &value){
				value_ = value;
				return *this;
			}

			template<typename U> bool operator==(const U& comp) const{
				return value_ == comp;
			}

			inline VectorAtomic<T> &operator++(){
				++value_;
				return *this;
			}

			inline T operator++(int){
				return value_++;
			}

			inline VectorAtomic<T> &operator--(){
				--value_;
				return *this;
			}

			inline T operator--(int){
				return value_--;
			}
		};
		template<typename T, typename U> bool operator==( U lhs, const VectorAtomic<T>& rhs){ // Compares the size_t value to the size
			return rhs == lhs;
		}
		template<typename T> void Acquire(std::vector<T> &receiver, std::vector<VectorAtomic<T>> &donator){
			receiver.resize(donator.size());
			for(auto i=donator.size(); i--; ){
				receiver.at(i) = donator.at(i);
			}
			donator.clear();
			donator.shrink_to_fit();
		}

		template<typename T> void SetDimensions(std::vector<T> &vec, size_t dim1, size_t dim2){
			vec.resize(dim1);
			for(auto i=dim1; i--; ){
				vec.at(i).resize(dim2);
			}
		}

		template<typename T> void SetDimensions(std::vector<T> &vec, size_t dim1, size_t dim2, size_t dim3){
			vec.resize(dim1);
			for(auto i=dim1; i--; ){
				SetDimensions(vec.at(i), dim2, dim3);
			}
		}

		constexpr uint64_t IntPow(uint16_t base, uint16_t power){
			uint64_t result = 1;
			for(auto i=power; i--; ){
				result *= base;
			}
			return result;
		}

		inline unsigned int TrueRandom(){
			std::random_device rd;
			return rd();
		}
	}
}

#endif // UTILITIES_H
