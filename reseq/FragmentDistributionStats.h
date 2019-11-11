#ifndef FRAGMENTDISTRIBUTIONSTATS_H
#define FRAGMENTDISTRIBUTIONSTATS_H

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <mutex>
#include <random>
#include <stdint.h>
#include <vector>

#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#include "nlopt.hpp"

#include "FragmentDuplicationStats.h"
#include "Reference.h"
#include "utilities.hpp"
#include "Vect.hpp"

namespace reseq{
	enum RefSeqBiasSimulation{
		kKeep,
		kNo,
		kDraw,
		kFile,
		kError
	};

	struct BiasCalculationParams{
		uintRefSeqId ref_seq_id;
		uintSeqLen fragment_length;
	};

	class BiasCalculationParamsSplitSeqs{
	public:
		uintRefSeqBin ref_seq_bin_;
		uintSeqLen fragment_length_;

		bool bias_calculation_;

		BiasCalculationParamsSplitSeqs():
			ref_seq_bin_(0),
			fragment_length_(0),
			bias_calculation_(false)
		{}

		BiasCalculationParamsSplitSeqs(uintRefSeqBin ref_seq_bin, uintSeqLen fragment_length, bool bias_calculation):
			ref_seq_bin_(ref_seq_bin),
			fragment_length_(fragment_length),
			bias_calculation_(bias_calculation)
		{}

		void Set(uintRefSeqBin ref_seq_bin, uintSeqLen fragment_length, bool bias_calculation){
			ref_seq_bin_ = ref_seq_bin;
			fragment_length_ = fragment_length;
			bias_calculation_ = bias_calculation;
		}

		void Clear(uintRefSeqBin ref_seq_bin){
			ref_seq_bin_ = ref_seq_bin; // It is important to have the correct one in here, because it may define the vector that is cleared and we don't want to clear an active one
			fragment_length_ = 0; // Signal that nothing has to be done
			bias_calculation_ = false;
		}
	};

	class BiasCalculationVectors{
	private:
		static uintPercent OptimizeSplineIdToKnot(uintNumFits id, const uintPercent range);
		static uintPercentShift OptimizeSplineIdToShift(uintNumFits id, const uintPercent range);

		// Used only for paper output (kParameterInfoFile, kDispersionInfoFile)
		static std::mutex file_mutex_;
	public:
		// Definitions
		static constexpr double kStopCriterion = 1e-6; // Minimum relative change in chi2 to finish fit
		static const uintNumFits kMaxLikelihoodCalculations = 100;
		static const uint16_t kSplinePrecisionFactor = 10; // Spline fit is fast so we can require a higher precision

		static const uintPercent kGCSplineDf = 6; // Number of knots for gc spline
		static const uintPercent kMaxKnotShift = 20; // Maximum shift for a single knot from one fit to the next for greedy knot adjustment
		static const uintPercent kPercentGCSitesForNormalization = 80; // Use this percent of gc values with highest number of sites for normalization to avoid effects from overfitted gc with low number of sites

		static const uintNumFits kNumFitsInsertLength = 30; // Number of insert length values with the highest number of counts fitted
		static constexpr double kNXXRefSeqs = 80.0; // Number of reference sequences to fit (so Nxx is reached)

		static constexpr double kLowerBound = -1e25;
		static constexpr double kUpperBound = 1e25;
		static constexpr double kBaseValue = 1e-25;

		// Definitions for paper (extra output)
		static const bool kSurMult = false; // This is only for the fit for the comparison plots, the rest of the program is using the better sum without a toggle, so do not switch this to true and use the simulator

		static constexpr const char *kParameterInfoFile = NULL; //"maxlike_fit.csv";
		static constexpr const char *kDispersionInfoFile = NULL; //"dispersion_fit.csv";
		static const uintDupCount kMaxDuplications = 100; // Maximum number of duplications used for dispersion fit reported in kDispersionInfoFile
		static const uintSeqLen kDispersionBinSize = 200000; // Bin size used to get sample mean and for dispersion fit reported in kDispersionInfoFile (2*sites, due to forward+reverse)

		// Variables used for fit
		// Input
		std::vector<FragmentSite> sites_;

		uintFragCount total_counts_;
		uintFragCount total_sites_;

		std::array<uintFragCount, 101> gc_count_;
		std::array<uintFragCount, 4*Surrounding::Length()> sur_count_;

		std::array<uintFragCount, 101> gc_sites_;
		std::array<uintFragCount, 4*Surrounding::Length()> sur_sites_;

		// Results
		std::array<uintPercent, kGCSplineDf> gc_knots_;
		std::array<double, 101> gc_bias_;

		std::array<double, 4*Surrounding::Length()> sur_bias_;

		std::vector<double> dispersion_;

		bool converged_;
		double loglike_;

		// Temporary calculation variables
		double loglike_pois_base_;

		std::array<double, 101> gc_bias_no_logit_;
		std::array<double, 101> gc_bias_grad_;
		std::array<std::pair<double, double>, 101> gc_bias_sum_;
		std::array<std::array<std::pair<double, double>, 4*Surrounding::Length()>, 101> grad_gc_bias_sum_;
		std::array<double, 4*Surrounding::Length()> sur_grad_;

		std::array<std::array<std::array<double, kGCSplineDf>, kGCSplineDf-1>, 3> lin_comb_gc_splines_;

		// Optimizer variables
		std::vector<double> fit_pars_;
		std::vector<double> bounds_;
		nlopt::opt optimizer_poisson_;
		nlopt::opt optimizer_nbinom_;

		std::vector<double> gc_spline_pars_;
		nlopt::opt optimizer_spline_;

		// Variables only used for paper output
		uintRefSeqBin ref_seq_bin_;
		uintSeqLen insert_length_;

		uintNumFits func_calls_;

		uintSeqLen start_site_dispersion_fit_;
		uintSeqLen end_site_dispersion_fit_;
		double use_sample_mean_;
		std::array<uintFragCount, kMaxDuplications+2> duplication_count_part_;
		uintNumFits func_calls_const_disp_;

		// Functions for fit
		BiasCalculationVectors():
			optimizer_poisson_(nlopt::LD_LBFGS,4*Surrounding::Length() + (!kSurMult?1:0) ),
			optimizer_nbinom_(nlopt::LD_LBFGS, 4*Surrounding::Length() + (!kSurMult?1:0) + kGCSplineDf + 1 + 2),
			optimizer_spline_(nlopt::LD_LBFGS, kGCSplineDf+1)
		{
			dispersion_.reserve(2);

			fit_pars_.reserve(4*Surrounding::Length() + (!kSurMult?1:0) + kGCSplineDf + 1 + 2);
			bounds_.reserve(fit_pars_.capacity());
			gc_spline_pars_.resize(1+kGCSplineDf);

			optimizer_poisson_.set_ftol_rel(kStopCriterion);
			optimizer_poisson_.set_maxeval(kMaxLikelihoodCalculations);

			optimizer_nbinom_.set_ftol_rel(kStopCriterion);
			optimizer_nbinom_.set_maxeval(kMaxLikelihoodCalculations);

			optimizer_spline_.set_ftol_rel(kStopCriterion / kSplinePrecisionFactor);
			optimizer_spline_.set_maxeval(kMaxLikelihoodCalculations * kSplinePrecisionFactor);

			bounds_.resize(gc_spline_pars_.size(), kUpperBound);
			optimizer_spline_.set_upper_bounds(bounds_);

			bounds_.clear();
			bounds_.resize(gc_spline_pars_.size(), kLowerBound);
			bounds_.at(0) = kBaseValue;
			optimizer_spline_.set_lower_bounds(bounds_);
		}

		void AddCountsFromSite(const FragmentSite &site, std::array<uintFragCount, 101> &gc_count, std::array<uintFragCount, 4*Surrounding::Length()> &sur_count);
		void GetCounts();
		void RemoveUnnecessarySites();
		void AddBaseValue(uintDupCount count);
		void GetLogLikeBase();
		void DeactivateZeroCounts(std::vector<double> &par, double deactive_value, uintNumFits sur_shift);

		void NormGC();
		void NormSurroundings(const std::vector<double> &x);
		void UnnormSurroundingGradients(std::vector<double> &grad, const std::vector<double> &x);

		double SurBiasAtSite(const std::pair<double, double>& bias);
		std::pair<double, double> SurBiasAtSiteSplit(const FragmentSite& site);
		double SurBiasAtSite(const FragmentSite& site);

		void DefineStartingKnots();
		void PrepareSplines();
		void GetSplineCoefficients(double &a, double &b, double &c, double &d, uintPercent k, const std::vector<double> &spline_pars);
		void CalculateSpline(const std::vector<double> &spline_pars);
		void CalculateSplineGrad(const std::vector<double> &spline_pars, std::vector<double> &grad, uintNumFits grad_offset);
		static double LogLikeGcSpline(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
		double OptimizeSpline();
		double OptimizeSpline(uintPercent knot, uintPercentShift shift);
		void GetGCSpline();

		static double LogLikelihoodPoisson(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
		static double GetDispersion(double bias, double a, double b);
		void LogLikelihoodNbinomSite(double &loglike, double &a, double &b, double &ga, double &gb, std::vector<double> &grad, const FragmentSite &site);
		static double LogLikelihoodNbinom(const std::vector<double> &x, std::vector<double> &grad, void* f_data);

		void PostProcessFit();

		// Functions for paper output
		static double LogLikelihoodConstDispersion(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
		void WriteOutParameterInfo(const char *context);
		void WriteOutDispersionInfo(const char *context);
	};

	class FragmentDistributionStats{
	public:
		class ThreadData{
		private:
			std::vector<uintSeqLen> num_sites_per_insert_length_;
			std::vector<std::pair<uintSeqLen, std::pair<uintRefSeqBin, uintSeqLen>>> bias_calc_tmp_params_;
			BiasCalculationVectors bias_calc_vects_;

			friend class FragmentDistributionStats;

			// Google test
			friend class FragmentDistributionStatsTest;
		public:
			ThreadData( uintSeqLen maximum_insert_length, uintSeqLen max_seq_bin_len ){
				bias_calc_tmp_params_.reserve( maximum_insert_length );
				bias_calc_vects_.sites_.reserve(max_seq_bin_len);
			}
		};
	private:
		// Definitions
		const uintSeqLen kOutskirtRange = 20; // 20 bases before and after each fragment are reported
		uintPercent surrounding_range_; // Get it from Reference class
		static const uint16_t kMaxBinsQueuedForBiasCalc = 5; // Defines length of the parameter vector which is used to feed calculation threads
		const uintErrorCount kMaxErrorsShownPerFile = 50;

		// User parameter
		uintSeqLen max_ref_seq_bin_length_;

		// Flags
		bool calculate_bias_; // Only deactivated to speed up tests

		// Temporary variables
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_abundance_;
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_insert_lengths_;
		std::vector<utilities::VectorAtomic<uintFragCount>> tmp_gc_fragment_content_;
		SurroundingCountAtomic tmp_fragment_surroundings_;

		std::array<std::array<std::vector<utilities::VectorAtomic<uintFragCount>>, 4>, 2> tmp_outskirt_content_;

		std::vector<std::vector<std::pair<uintRefLenCalc, uintSeqLen>>> fragment_sites_by_ref_seq_bin_; // fragment_sites_by_ref_seq_bin_[ReferenceSequenceBinId][UniqueId] = {startPosition*2 + 0/1(forward/reverse),fragment_length}
		std::vector<utilities::VectorAtomic<uintFragCount>> fragment_sites_by_ref_seq_bin_cur_id_; // fragment_sites_by_ref_seq_bin_cur_id_[ReferenceSequenceBinId] = CurrentUniqueId

		std::vector<std::vector<std::vector<uintRefLenCalc>>> fragment_sites_by_ref_seq_bin_by_insert_length_; // fragment_sites_by_ref_seq_bin_by_insert_length_[ReferenceSequenceBinId][FragmentLength][UniqueId] = startPosition*2 + 0/1(forward/reverse)
		std::atomic<uintRefSeqBin> num_handled_reference_sequence_bins_; // num_handled_reference_sequence_bins_ = #ofAlreadyHandledReferenceSequences = IdOfFirstUnhandledReferenceSequence

		std::vector<bool> ref_seq_in_nxx_;
		std::vector<uintRefSeqBin> ref_seq_start_bin_;
		std::array<std::vector<double>, 101> tmp_gc_bias_; // tmp_gc_bias_[GC][#Fit]
		std::array<std::vector<double>, 4*Surrounding::Length()> tmp_sur_bias_; // tmp_sur_bias_[SurBase][#Fit]
		std::array<std::vector<double>, 2> tmp_dispersion_parameters_; // tmp_dispersion_parameters_[dispPar][#Fit]

		std::array<std::atomic_flag, kMaxBinsQueuedForBiasCalc> claimed_bias_bins_;
		std::array<std::atomic<uintNumFits>, kMaxBinsQueuedForBiasCalc> current_bias_param_;
		std::array<std::atomic<uintNumFits>, kMaxBinsQueuedForBiasCalc> finished_bias_calcs_;
		std::atomic<uintNumFits> current_bias_result_;
		std::vector<BiasCalculationParamsSplitSeqs> bias_calc_params_;
		std::atomic<uintNumFits> params_left_for_calculation_;
		std::atomic<uintNumFits> params_fitted_;

		// Collected variables for bias calculation
		std::vector<uintFragCount> abundance_; // abundance_[referenceID] = #numberOfPairsMapToIt
		Vect<uintFragCount> insert_lengths_; // insert_lengths_[length] = #pairs
		Vect<uintFragCount> gc_fragment_content_; // gc_fragment_content_[gcContent(%)] = #fragments
		SurroundingCount fragment_surroundings_;

		// Collected variables for dispersion calculation
		Vect<Vect<uintFragCount>> site_count_; // site_count_[GCcontent(%)][FragmentLength] = #SitesInReference

		// Collected variables for plotting
		std::array<std::array<Vect<uintFragCount>, 4>, 2> outskirt_content_; // outskirt_content_[forward/reverse][refBase][position] = #(reads with given reference content at given position before or after read)

		// Calculated bias variables
		std::vector<double> ref_seq_bias_; // ref_seq_bias_[referenceID] = #numberOfPairsMapToIt/#PossibilitiesInReference
		Vect<double> insert_lengths_bias_; // insert_lengths_bias_[length] = #pairs/#PossibilitiesInReference
		Vect<double> gc_fragment_content_bias_; // gc_fragment_content_bias_[gcContent(%)] = #fragments/#
		SurroundingBias fragment_surroundings_bias_;

		std::array<double, 2> dispersion_parameters_;

		// Calculated variables for plotting
		std::array<std::vector<double>, 4> fragment_surrounding_bias_by_base_;

		// Helper functions
		uintSeqLen RefSeqSplitLength(uintRefSeqId ref_seq_id, const Reference &reference){
			return reference.SequenceLength(ref_seq_id)/(reference.SequenceLength(ref_seq_id)/max_ref_seq_bin_length_+1);
		}
		uintRefSeqId GetRefSeqId(uintRefSeqBin ref_seq_bin){
			// We shouldn't have many super long reference sequences, so that ref_seq_bin will be close to ref_seq_id and this should be fairly efficient
			uintRefSeqId ref_seq_id = std::min(ref_seq_bin, static_cast<uintRefSeqBin>(ref_seq_start_bin_.size())-1);
			while(ref_seq_bin < ref_seq_start_bin_.at(ref_seq_id) ){
				--ref_seq_id;
			}
			return ref_seq_id;
		}
		inline void IncreaseErrorCounter(uintErrorCount &errors);
		void PrepareBiasCalculation( const Reference &ref, uintSeqLen maximum_insert_length, const std::vector<uintFragCount> &reads_per_ref_seq_bin );
		void SortFragmentSites(uintRefSeqBin ref_seq_bin, std::vector<uintSeqLen> &num_sites_per_insert_length);
		void UpdateBiasCalculationParams(uintRefSeqBin ref_seq_bin, uint32_t queue_spot, std::vector<std::pair<uintSeqLen, std::pair<uintRefSeqBin, uintSeqLen>>> &tmp_params, std::mutex &print_mutex );
		void FillParams(std::vector<BiasCalculationParams> &params, const Reference &reference) const;

		void CountDuplicates(FragmentDuplicationStats &duplications, const BiasCalculationParamsSplitSeqs &params, const Reference &reference);
		void AddFragmentsToSites(std::vector<FragmentSite> &sites, const std::vector<uintRefLenCalc> &fragment_positions, uintSeqLen min_dist_to_ref_seq_ends);
		void CalculateBiasByBin(BiasCalculationVectors &tmp_calc, const Reference &reference, FragmentDuplicationStats &duplications, uintRefSeqBin ref_seq_bin, uintSeqLen insert_length);
		void AcquireBiases(const BiasCalculationVectors &calc, std::mutex &print_mutex);
		bool StoreBias();
		void CalculateInsertLengthAndRefSeqBias(const Reference &reference, uintNumThreads num_threads, std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc);
		void AddNewBiasCalculations(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, std::mutex &print_mutex);
		void ExecuteBiasCalculations( const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, std::mutex &print_mutex );
		void HandleReferenceSequencesUntil(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);
		static void BiasSumThread( const FragmentDistributionStats &self, const Reference &reference, const std::vector<BiasCalculationParams> &params, std::atomic<uintNumFits> &current_param, std::vector<double> &insert_length_sum, std::vector<double> &ref_seq_sum, std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc, std::mutex &result_mutex );

		static void BiasNormalizationThread( const FragmentDistributionStats &self, const Reference &reference, const std::vector<BiasCalculationParams> &params, std::atomic<uintNumFits> &current_param, double &norm, std::mutex &result_mutex, double &max_bias );
		
		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & abundance_;
			ar & insert_lengths_;
			ar & gc_fragment_content_;
			ar & fragment_surroundings_;

			ar & site_count_;

			ar & outskirt_content_;

			ar & ref_seq_bias_;
			ar & insert_lengths_bias_;
			ar & gc_fragment_content_bias_;
			ar & fragment_surroundings_bias_;

			ar & dispersion_parameters_;
		}

		// Google test
		friend class FragmentDistributionStatsTest;
		FRIEND_TEST(FragmentDistributionStatsTest, UpdateRefSeqBias);
		friend class SimulatorTest;
	public:
		FragmentDistributionStats():
			max_ref_seq_bin_length_(0), // Will be set later as it is a program parameter
			calculate_bias_(true), // Not calculated only to speed up tests
			num_handled_reference_sequence_bins_(0),
			current_bias_result_(0),
			dispersion_parameters_({0, 1e-100}) // Initialize with poisson to avoid writing random values into file if saving without fitting
			{
		}

		// Getter functions
		const std::vector<uintFragCount> &Abundance() const{ return abundance_; }
		const Vect<uintFragCount> &InsertLengths() const{ return insert_lengths_; }
		const Vect<uintFragCount> &GCFragmentContent() const{ return gc_fragment_content_; }

		const Vect<Vect<uintFragCount>> &SiteCount() const{ return site_count_; }

		const std::vector<double> &RefSeqBias() const{ return ref_seq_bias_; }
		const Vect<double> &InsertLengthsBias() const{ return insert_lengths_bias_; }
		const Vect<double> &GCFragmentContentBias() const{ return gc_fragment_content_bias_; }
		
		const Vect<uintFragCount> &OutskirtContent(uintTempSeq direction, uintBaseCall nucleotide) const{ return outskirt_content_.at(direction).at(nucleotide); }
		const std::vector<double> &FragmentSurroundingBiasByBase(uintBaseCall nucleotide) const{ return fragment_surrounding_bias_by_base_.at(nucleotide); }

		// Setter functions
		void AddAbundance(uintRefSeqId ref_seq_id){ ++tmp_abundance_.at(ref_seq_id); }
		void AddInsertLengths(uintSeqLen length){ ++tmp_insert_lengths_.at(length); }
		void AddGCContent(uintPercent gc){ ++tmp_gc_fragment_content_.at(gc); }
		void AddFragmentSite(uintRefSeqId ref_seq_id, uintSeqLen length, uintSeqLen position, uintTempSeq template_segment, const Reference &reference){
			auto ref_seq_bin = GetRefSeqBin(ref_seq_id,position,reference);
			fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin)++) = {(static_cast<uintRefLenCalc>(position)<<1)+template_segment, length};
		}

		void AddThreadData(ThreadData &thread){
			thread.num_sites_per_insert_length_.clear();
			thread.num_sites_per_insert_length_.shrink_to_fit();
		}

		void ActivateBiasCalculation(){ calculate_bias_ = true; }
		void DeactivateBiasCalculation(){ calculate_bias_ = false; }

		void SetUniformBias();

		// Main functions
		uintRefSeqBin CreateRefBins( const Reference &ref, uintSeqLen max_ref_seq_bin_size );
		uintRefSeqBin GetRefSeqBin(uintRefSeqId ref_seq_id, uintSeqLen position, const Reference &reference){
			return ref_seq_start_bin_.at(ref_seq_id) + position/RefSeqSplitLength(ref_seq_id, reference);
		}
		uintSeqLen MaxRefSeqBinLength(const Reference &reference){
			uintSeqLen max(0);
			for(uintRefSeqId seq=0; seq < reference.NumberSequences(); ++seq){
				utilities::SetToMax(max, RefSeqSplitLength(seq, reference));
			}
			return max;
		}

		void Prepare( const Reference &ref, uintSeqLen maximum_insert_length, const std::vector<uintFragCount> &reads_per_ref_seq_bin );
		void FillInOutskirtContent( const Reference &reference, const seqan::BamAlignmentRecord &record_start, uintSeqLen fragment_start_pos, uintSeqLen fragment_end_pos );

		void HandleReferenceSequencesUntil(uintRefSeqId still_needed_reference_sequence, uintSeqLen still_needed_position, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);
		void FinishThreads(ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);

		void Finalize();
		void Shrink();
		bool FinalizeBiasCalculation(const Reference &reference, uintNumThreads num_threads, FragmentDuplicationStats &duplications);
		bool UpdateRefSeqBias(RefSeqBiasSimulation model, const std::string &bias_file, const Reference &ref, std::mt19937_64 &rgen);
		double CalculateBiasNormalization(double &max_bias, const Reference &reference, uintNumThreads num_threads, uintFragCount total_reads) const;
		double CalculateNonZeroThreshold(double bias_normalization, double max_bias) const;

		static uintDupCount NegativeBinomial(double p, double r, double probability_chosen);
		double Dispersion(double bias) const{ return BiasCalculationVectors::GetDispersion( bias, dispersion_parameters_.at(0), dispersion_parameters_.at(1) ); }
		uintDupCount GetFragmentCounts(double bias_normalization, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintPercent gc, const Surrounding &fragment_start, const Surrounding &fragment_end, double probability_chosen, double non_zero_threshold=0.0) const;

		void PreparePlotting();
	};
}

#endif // FRAGMENTDISTRIBUTIONSTATS_H
