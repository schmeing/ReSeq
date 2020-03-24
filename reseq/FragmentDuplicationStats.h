#ifndef FRAGMENTDUPLICATIONSTATS_H
#define FRAGMENTDUPLICATIONSTATS_H

#include <algorithm>
#include <deque>
#include <limits>
#include <list>
#include <cmath>
#include <set>
#include <stdint.h>
#include <vector>

#include <seqan/bam_io.h>

#include "Reference.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	class FragmentDuplicationStats{
	public:
		static const uintDupCount kMaxDuplication = 50; // Maximum duplication for direct plotting and plotted dispersion fit

	private:
		// Temporary variables
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_duplication_number_; // tmp_duplication_number_[#FragmentsWithSamePosition] = #FragmentSites

		// Collected variables for plotting
		Vect<uintFragCount> duplication_number_; // duplication_number_[#FragmentsWithSamePosition] = #Fragments
		
		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & duplication_number_;
		}

		// Google test
		friend class FragmentDuplicationStatsTest;
		FRIEND_TEST(FragmentDuplicationStatsTest, DispersionCalculation);
		FRIEND_TEST(FragmentDistributionStatsTest, BiasBinningAndFragmentCounts);
	public:
		// Getter functions
		const Vect<uintFragCount> &DuplicationNumber() const{ return duplication_number_; }

		// Main functions
		inline void PrepareTmpDuplicationVector(){
			tmp_duplication_number_.resize(kMaxDuplication+2);
		}

		inline void AddDuplication(uintDupCount dup){
			if(dup > kMaxDuplication){
				++tmp_duplication_number_.at(kMaxDuplication+1);
			}
			else{
				++tmp_duplication_number_.at(dup);
			}
		}
		inline void AddSites(std::vector<FragmentSite> &sites){
			for(auto &site : sites){
				if(0.0 != site.bias_){
					AddDuplication(site.count_forward_);
					AddDuplication(site.count_reverse_);
				}
			}
		}
		void AddDuplicates( std::vector<uintSeqLen> &fragment_positions );
		void FinalizeDuplicationVector();
	};
}

#endif // FRAGMENTDUPLICATIONSTATS_H
