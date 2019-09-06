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
#include "utilities.h"
#include "Vect.hpp"

namespace reseq{
	class FragmentDuplicationStats{
	public:
		static const uint16_t max_duplication = 50;

	private:
		// Temporary variables
		std::vector<std::array<std::array<utilities::VectorAtomic<uint64_t>, max_duplication+2>, 101>> tmp_duplication_number_by_gc_insert_length_; // tmp_duplication_number_by_gc_insert_length_[InsertLength][GC][#Fragments] = #sites

		// Collected variables for plotting
		Vect<Vect<Vect<uint64_t>>> duplication_number_by_gc_insert_length_; // duplication_number_by_gc_insert_length_[GC][InsertLength][#Fragments] = #Sites

		// Calculated variables for plotting
		Vect<uint64_t> duplication_number_; // duplication_number_[#ReadsWithSamePosition] = #fragments
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

		const Vect<uint64_t> &DuplicationNumber() const{ return duplication_number_; }

		// Main functions
		inline void PrepareTmpDuplicationVector(std::vector<uint32_t>::size_type maximum_insert_length){
			tmp_duplication_number_by_gc_insert_length_.resize(maximum_insert_length+1);
		}
		inline void AddDuplication(uint32_t insert_length, uint16_t gc, uint16_t dup){
			if(dup > max_duplication){
				++tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(max_duplication+1);
			}
			else{
				++tmp_duplication_number_by_gc_insert_length_.at(insert_length).at(gc).at(dup);
			}
		}
		inline void AddSites(std::vector<FragmentSite> &sites, uint32_t insert_length){
			for(auto &site : sites){
				if(0.0 != site.bias_){
					AddDuplication(insert_length, site.gc_, site.count_forward_);
					AddDuplication(insert_length, site.gc_, site.count_reverse_);
				}
			}
		}
		void AddDuplicates( std::vector<uint32_t> &fragment_positions, uint32_t ref_seq_id, uint32_t insert_length, const Reference &reference );
		void FinalizeDuplicationVector(const std::vector<std::vector<utilities::VectorAtomic<uint64_t>>> &site_count_by_insert_length_gc);

		void CalculateDispersionPlot();

		void PreparePlotting();
		void PrepareTesting();
	};
}

#endif // FRAGMENTDUPLICATIONSTATS_H
