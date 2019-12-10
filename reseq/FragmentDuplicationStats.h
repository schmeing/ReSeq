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
		std::vector<std::array<std::array<utilities::VectorAtomic<uintFragCount>, kMaxDuplication+2>, 101>> tmp_duplication_number_by_gc_insert_length_; // tmp_duplication_number_by_gc_insert_length_[InsertLength][GC][#Fragments] = #sites

		// Collected variables for plotting
		Vect<Vect<Vect<uintFragCount>>> duplication_number_by_gc_insert_length_; // duplication_number_by_gc_insert_length_[GC][InsertLength][#Fragments] = #Sites

		// Calculated variables for plotting
		Vect<uintFragCount> duplication_number_; // duplication_number_[#ReadsWithSamePosition] = #fragments
		std::vector<double> mean_list_;
		std::vector<double> dispersion_list_;
		
		// Helper functions
		inline void SumUpDuplicationNumber();

		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & duplication_number_by_gc_insert_length_;
		}

		// Google test
		friend class FragmentDuplicationStatsTest;
		FRIEND_TEST(FragmentDuplicationStatsTest, DispersionCalculation);
		FRIEND_TEST(FragmentDistributionStatsTest, BiasBinningAndFragmentCounts);
	public:
		// Getter functions
		const std::vector<double> &DispersionList() const{ return dispersion_list_; }
		const std::vector<double> &MeanList() const{ return mean_list_; }

		const Vect<uintFragCount> &DuplicationNumber() const{ return duplication_number_; }

		// Main functions
		inline void PrepareTmpDuplicationVector(uintSeqLen maximum_insert_length){
			tmp_duplication_number_by_gc_insert_length_.resize(maximum_insert_length+1);
		}
		inline void AddDuplication(uintSeqLen insert_length, uintPercent gc, uintDupCount dup){
			if(dup > kMaxDuplication){
				++tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(kMaxDuplication+1);
			}
			else{
				++tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(dup);
			}
		}
		inline void AddSites(std::vector<FragmentSite> &sites, uintSeqLen insert_length){
			for(auto &site : sites){
				if(0.0 != site.bias_){
					AddDuplication(insert_length, site.gc_, site.count_forward_);
					AddDuplication(insert_length, site.gc_, site.count_reverse_);
				}
			}
		}
		void AddDuplicates( std::vector<uintSeqLen> &fragment_positions, uintRefSeqId ref_seq_id, uintSeqLen insert_length, uintSeqLen base_pos, const Reference &reference );
		void FinalizeDuplicationVector(const std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc);

		void CalculateDispersionPlot();

		void PreparePlotting();
		void PrepareTesting();
	};
}

#endif // FRAGMENTDUPLICATIONSTATS_H
