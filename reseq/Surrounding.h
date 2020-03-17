#ifndef SURROUNDING_H
#define SURROUNDING_H

#include <array>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "utilities.hpp"
#include "SurroundingBase.hpp"

namespace reseq{
	class Surrounding : public SurroundingBase<3, 10, 10, int32_t>{
		friend class SurroundingCountAtomic;
		friend class SurroundingCount;
		friend class SurroundingBias;

		// Google test
		friend class FragmentDistributionStatsTest;
		friend class ReferenceTest;
		friend class SimulatorTest;
		friend class SurroundingTest;

	public:
		void ChangeSurroundingBase( uintSurPos pos, seqan::Dna new_base );
		void DeleteSurroundingBaseShiftingOnRightSide( uintSurPos pos, seqan::Dna new_end_base );
		void DeleteSurroundingBaseShiftingOnLeftSide( uintSurPos pos, seqan::Dna new_end_base );
		void InsertSurroundingBasesShiftingOnRightSide( uintSurPos pos, seqan::DnaString new_bases );
		void InsertSurroundingBasesShiftingOnLeftSide( uintSurPos pos, seqan::DnaString new_bases );
	};

	class SurroundingCountAtomic{
	private:
		std::array< std::vector<utilities::VectorAtomic<uintFragCount>>, Surrounding::kNumBlocks> counts_; // count_[BlockNumber][(262144*nucAt0+...+(256*nucAt5)+(64*nucAt6)+(16*nucAt7)+(4*nucAt8)+(nucAt9)] = #fragments (always facing inwards, Block0 is outside fragment, Block1 first in fragment ..., value 262144 is first base in block/fragment(block1), 256 is 5th base in block, ...)

		friend class SurroundingCount;
	public:
		SurroundingCountAtomic(){
			for(auto block = Surrounding::kNumBlocks; block--;){
				counts_.at(block).resize(Surrounding::Size());
			}
		}

		void Count(const Surrounding &sur){
			for( auto block = Surrounding::kNumBlocks; block--; ){
				if(0 <= sur.sur_.at(block)){
					++(counts_.at(block).at( sur.sur_.at(block) ));
				}
			}
		}
	};

	class SurroundingCount{
	private:
		std::array< std::vector<uintFragCount>, Surrounding::kNumBlocks> counts_;

		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & counts_;
		}

		// Google test
		friend class FragmentDistributionStatsTest;
	public:
		SurroundingCount(){
			for(auto block = Surrounding::kNumBlocks; block--;){
				counts_.at(block).resize(Surrounding::Size());
			}
		}

		void Acquire(SurroundingCountAtomic &tmp_counts){
			for(auto block = Surrounding::kNumBlocks; block--; ){
				utilities::Acquire(counts_.at(block), tmp_counts.counts_.at(block));
			}
		}
	};

	class SurroundingBias{
	private:
		std::array<std::vector<double>, Surrounding::kNumBlocks> bias_;

		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & bias_;
		}

		// Google test
		friend class FragmentDistributionStatsTest;
		friend class ReferenceTest;
		friend class SurroundingTest;
	public:
		SurroundingBias(){
			for(auto block = Surrounding::kNumBlocks; block--;){
				bias_.at(block).resize(Surrounding::Size(), 0.0);
			}
		}

		void SetUniform(){
			for(auto block = Surrounding::kNumBlocks; block--;){
				bias_.at(block).clear();
				bias_.at(block).resize(Surrounding::Size(), 1.0);
			}
		}

		void CombinePositions(const std::array<double, 4*Surrounding::Length()> &separated);
		void SeparatePositions(std::array<double, 4*Surrounding::Length()> &separated) const;

		double Bias(const Surrounding &sur) const{
			double bias = 0.0;
			for(auto block = Surrounding::kNumBlocks; block--;){
				bias += bias_.at(block).at( sur.sur_.at(block) );
			}
			return utilities::InvLogit2(bias);
		}
	};

#ifndef SWIG // This part is not needed for the python plotting
	template <uintSurPos K> class Kmer : public SurroundingBase<1, K, 0, int64_t>{
		template<uintSurPos T> friend class KmerCount;
	};

	template <uintSurPos K> class KmerCount{
	public:
		std::vector<uintNucCount> counts_;

		void Prepare(){
			counts_.clear();
			counts_.resize(Kmer<K>::Size(), 0);
		}

		void Clear(){
			counts_.clear();
			counts_.shrink_to_fit();
		}

		void Count( const Kmer<K> &kmer ){
			if(0 <= kmer.sur_.at(0)){
				++(counts_.at( kmer.sur_.at(0) ));
			}
		}

		void Revert( const Kmer<K> &kmer ){
			if(0 <= kmer.sur_.at(0)){
				--(counts_.at( kmer.sur_.at(0) ));
			}
		}

		typename Kmer<K>::intType CountSequenceForward(const seqan::IupacString &seq, uintReadLen start_pos){
			Kmer<K> kmer;
			if(start_pos+kmer.Length() <= length(seq)){
				kmer.ForwardWithN(seq, start_pos);
				auto start_kmer = kmer.sur_.at(0);
				Count(kmer);

				for(uintSeqLen pos=start_pos+1; pos+kmer.Length() <= length(seq); ++pos){
					kmer.UpdateForwardWithN(seq, pos);
					if(kmer.sur_.at(0) == start_kmer){
						// If we hit the start kmer again, we are probably in a segment with short repeats, so we remove counts again and stop
						for(uintSeqLen pos2=start_pos+1; pos2 < pos; ++pos2){
							kmer.UpdateForwardWithN(seq, pos2);
							Revert(kmer);
						}
						return -1;
					}
					Count(kmer);
				}

				return start_kmer; // Return first k-mer in sequence
			}

			return -1;
		}

		typename Kmer<K>::intType CountSequenceReverse(const seqan::IupacString &seq, uintReadLen end_pos){
			Kmer<K> kmer;
			if(kmer.Length() <= end_pos){
				kmer.ReverseWithN(seq, kmer.Length()-1);
				Count(kmer);

				for(uintSeqLen pos=kmer.Length(); pos < end_pos; ++pos){
					kmer.UpdateReverseWithN(seq, pos);
					Count(kmer);
				}

				return kmer.sur_.at(0); // Return first k-mer in reversed sequence
			}

			return -1;
		}

		void CountSequenceCanonical(const seqan::Dna5String &seq){
			Kmer<K> kmer_forward, kmer_reverse;
			kmer_forward.ForwardWithN(seq, 0);
			Count(kmer_forward);
			kmer_reverse.ReverseWithN(seq, kmer_forward.Length()-1);
			if(kmer_forward.sur_.at(0) != kmer_reverse.sur_.at(0)){
				Count(kmer_reverse);
			}

			for(uintSeqLen pos=1; pos <= length(seq)-kmer_forward.Length(); ++pos){
				kmer_forward.UpdateForwardWithN(seq, pos);
				Count(kmer_forward);
				kmer_reverse.UpdateReverseWithN(seq, pos + kmer_forward.Length()-1);
				if(kmer_forward.sur_.at(0) != kmer_reverse.sur_.at(0)){
					// If the reverse of the kmer is the kmer itself, don't add it twice
					Count(kmer_reverse);
				}
			}
		}
	};
#endif //SWIG
}

#endif // SURROUNDING_H
