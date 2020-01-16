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

		static constexpr double kDistanceWeightReductionFactor = 0.5;
		static const uintNumFits kNumFitsInsertLength = 30; // Number of insert length values with the highest number of counts fitted
		static constexpr double kNXXRefSeqs = 80.0; // Number of reference sequences to fit (so Nxx is reached)

		static constexpr double kLowerBound = -1e25;
		static constexpr double kUpperBound = 1e25;
		static constexpr double kBaseValue = 1e-25;

		// Definitions for paper (extra output)
		static const bool kSurMult = false; // This is only for the fit for the comparison plots, the rest of the program is using the better sum without a toggle, so do not switch this to true and use the simulator

		static constexpr const char *kParameterInfoFile = NULL; // "maxlike_fit.csv";
		static constexpr const char *kDispersionInfoFile = NULL; // "dispersion_fit2.csv";
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

		std::array<double, 101> gc_weights_;

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
		void CalculateGCWeights();
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
			std::vector<uintFragCount> tmp_frag_count_;

			uintRefSeqBin last_bin_;
			uintSeqLen last_added_region_;
			uintSeqLen correction_;

			friend class FragmentDistributionStats;

			// Google test
			friend class FragmentDistributionStatsTest;
		public:
			ThreadData( uintSeqLen maximum_insert_length, uintSeqLen max_seq_bin_len, uintSeqLen start_exclusion_length ):
				last_bin_(0),
				last_added_region_(0),
				correction_(start_exclusion_length){
				bias_calc_tmp_params_.reserve( maximum_insert_length );
				bias_calc_vects_.sites_.reserve(max_seq_bin_len);
				tmp_frag_count_.reserve(max_seq_bin_len+maximum_insert_length+1);
			}
		};
	private:
		// Definitions
		const uintSeqLen kOutskirtRange = 20; // 20 bases before and after each fragment are reported
		static const uint16_t kMaxBinsQueuedForBiasCalc = 5; // Defines length of the parameter vector which is used to feed calculation threads
		const uintErrorCount kMaxErrorsShownPerFile = 50;
		const double kPrecisionAimRefSeqFragLengthFit = 1.0001;
		const uintNumFits kMaxIterationsRefSeqFragLengthFit = 200;
		const double kMinLowQFractionForExclusion = 0.1;
		const uintFragCount kMinHighQToSeparateLowQ = 10;
		const uintFragCount kMinSites = 25000;
		const uintFragCount kMinCounts = 100;
		const uintFragCount kMaxSitesPerCount = 10000;

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

		std::vector<std::vector<std::pair<uintSeqLen, uintSeqLen>>> fragment_sites_by_ref_seq_bin_; // fragment_sites_by_ref_seq_bin_[ReferenceSequenceBinId][UniqueId] = {startPositionCorrectedForExcludedRegionsRelativeToBinStart*2 + 0/1(forward/reverse),fragment_length}
		std::vector<utilities::VectorAtomic<uintFragCount>> fragment_sites_by_ref_seq_bin_cur_id_; // fragment_sites_by_ref_seq_bin_cur_id_[ReferenceSequenceBinId] = CurrentUniqueId

		std::vector<std::vector<std::vector<uintSeqLen>>> fragment_sites_by_ref_seq_bin_by_insert_length_; // fragment_sites_by_ref_seq_bin_by_insert_length_[ReferenceSequenceBinId][FragmentLength][UniqueId] = startPositionCorrectedForExcludedRegionsRelativeToBinStart*2 + 0/1(forward/reverse)
		std::atomic<uintRefSeqBin> num_handled_reference_sequence_bins_; // num_handled_reference_sequence_bins_ = #ofAlreadyHandledReferenceSequences = IdOfFirstUnhandledReferenceSequence

		std::vector<std::vector<uintSeqLen>> lowq_site_start_by_ref_seq_bin_; // lowq_site_start_by_ref_seq_bin_[ReferenceSequenceBinId][UniqueId] = startPositionCorrectedForExcludedRegionsRelativeToBinStart
		std::vector<utilities::VectorAtomic<uintFragCount>> lowq_site_start_by_ref_seq_bin_cur_id_; // lowq_site_start_by_ref_seq_bin_cur_id_[ReferenceSequenceBinId] = CurrentUniqueId
		std::vector<std::vector<uintSeqLen>> lowq_site_end_by_ref_seq_bin_; // lowq_site_end_by_ref_seq_bin_[ReferenceSequenceBinId][UniqueId] = endPositionCorrectedForExcludedRegionsRelativeToBinStart
		std::vector<utilities::VectorAtomic<uintFragCount>> lowq_site_end_by_ref_seq_bin_cur_id_; // lowq_site_end_by_ref_seq_bin_cur_id_[ReferenceSequenceBinId] = CurrentUniqueId
		std::vector<std::vector<std::pair<uintSeqLen,uintSeqLen>>> lowq_site_start_exclusion_;
		std::vector<std::vector<std::pair<uintSeqLen,uintSeqLen>>> lowq_site_end_exclusion_;

		std::atomic<uintRefLenCalc> excluded_lowq_positions_;
		std::atomic<uintRefLenCalc> excluded_lowq_regions_;
		std::vector<utilities::VectorAtomic<uintFragCount>> corrected_abundance_; // corrected_abundance_[RefSeqId] = AbundanceCorrectedForSitesWithLowQualityReads
		std::vector<utilities::VectorAtomic<uintFragCount>> filtered_sites_; // filtered_sites_[RefSeqId] = SitesWithoutLowQualityReads
		std::vector<utilities::VectorAtomic<uintSeqLen>> fragment_lengths_used_for_correction_; // fragment_lengths_used_for_correction_[RefSeqId] = FragmentLengthsThatWereUsedForLowQualityCorrection

		std::vector<bool> ref_seq_in_nxx_;
		std::vector<uintRefSeqBin> ref_seq_start_bin_; // ref_seq_start_bin_[RefSeqId] = First RefSeqBin
		std::vector<std::pair<uintRefSeqId,uintSeqLen>> ref_seq_bin_def_; // ref_seq_bin_def_[RefSeqBin] = {RefSeqId,StartPos}
		std::array<std::vector<std::pair<double,double>>, 101> tmp_gc_bias_; // tmp_gc_bias_[GC][#Fit] = {FittedBiasValue, WeightOfFit}
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
			auto corrected_length = reference.SequenceLength(ref_seq_id)-reference.SumExcludedBases(ref_seq_id);
			return utilities::Divide(corrected_length, utilities::DivideAndCeil(corrected_length, max_ref_seq_bin_length_));
		}
		uintRefSeqId GetRefSeqId(uintRefSeqBin ref_seq_bin){
			return ref_seq_bin_def_.at(ref_seq_bin).first;
		}
		inline void IncreaseErrorCounter(uintErrorCount &errors);
		void PrepareBiasCalculation( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const std::vector<uintFragCount> &reads_per_frag_len_bin, const std::vector<uintFragCount> &lowq_reads_per_frag_len_bin );
		inline void CheckLowQExclusion(std::vector<std::pair<uintSeqLen,uintSeqLen>> &lowq_site_exclusion, uintSeqLen excluded_pos, uintFragCount lowq_count, const std::vector<uintFragCount> &tmp_frag_count );
		void CheckLowQExclusions(uintRefSeqBin ref_seq_bin, std::vector<uintFragCount> &tmp_frag_count, const Reference &reference);
		void SortFragmentSites(uintRefSeqBin ref_seq_bin, std::vector<uintSeqLen> &num_sites_per_insert_length);
		void UpdateBiasCalculationParams(uintRefSeqBin ref_seq_bin, uint32_t queue_spot, std::vector<std::pair<uintSeqLen, std::pair<uintRefSeqBin, uintSeqLen>>> &tmp_params, std::mutex &print_mutex );
		void FillParams(std::vector<BiasCalculationParams> &params, const Reference &reference) const;
		void FillParamsSimulation(std::vector<BiasCalculationParams> &params) const;

		void CountEndExclusion(uintFragCount &lowq_sites, uintSeqLen &cur_end_exclusion, uintSeqLen last_cut, uintSeqLen new_cut, uintRefSeqBin ref_seq_bin );
		void CountDuplicates(FragmentDuplicationStats &duplications, const BiasCalculationParamsSplitSeqs &params, const Reference &reference);
		void AddFragmentsToSites(std::vector<FragmentSite> &sites, const std::vector<uintSeqLen> &fragment_positions);
		uintSeqLen ExcludedLowQBases(const std::vector<std::pair<uintSeqLen,uintSeqLen>> &exlcuded_regions){
			uintSeqLen num_excluded(0);
			for(auto &exclusion : exlcuded_regions){
				num_excluded += exclusion.second - exclusion.first;
			}
			return num_excluded;
		}

		uintSeqLen ExcludeLowQualitySites(std::vector<FragmentSite> &sites, uintRefSeqBin ref_seq_bin, uintSeqLen insert_length);
		void CalculateBiasByBin(BiasCalculationVectors &tmp_calc, const Reference &reference, FragmentDuplicationStats &duplications, uintRefSeqBin ref_seq_bin, uintSeqLen insert_length);
		void AcquireBiases(const BiasCalculationVectors &calc, std::mutex &print_mutex);
		bool StoreBias();
		double MedianFragmentCoverage(const Reference &ref) const;
		void CalculateInsertLengthAndRefSeqBias(const Reference &reference, uintNumThreads num_threads, std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc);
		void ReplaceUncertainCorrectedAbundanceWithMedian(const Reference &ref);
		void AddNewBiasCalculations(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, std::mutex &print_mutex, const Reference &reference);
		void ExecuteBiasCalculations( const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, std::mutex &print_mutex );
		void HandleReferenceSequencesUntil(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);
		static void BiasSumThread( const FragmentDistributionStats &self, const Reference &reference, const std::vector<BiasCalculationParams> &params, std::atomic<uintNumFits> &current_param, std::vector<std::vector<double>> &bias_sum, std::vector<std::vector<utilities::VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc, std::mutex &print_mutex );

		uintRefSeqId SplitCoverageGroups(std::vector<uintRefSeqId> &coverage_groups) const;
		static void BiasNormalizationThread( const FragmentDistributionStats &self, const Reference &reference, const std::vector<BiasCalculationParams> &params, std::atomic<uintNumFits> &current_param, double &norm, std::mutex &result_mutex, const std::vector<uintRefSeqId> &coverage_groups, std::vector<std::vector<double>> &max_bias );
		double CalculateNonZeroThreshold(double bias_normalization, double max_bias) const;
		
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
		void AddFragmentSite(uintRefSeqId ref_seq_id, uintSeqLen length, uintSeqLen position, uintTempSeq template_segment, const Reference &reference, ThreadData &thread){
			auto bin_coords = GetRefSeqBin(ref_seq_id, position, reference, thread);
			fragment_sites_by_ref_seq_bin_.at(bin_coords.first).at(fragment_sites_by_ref_seq_bin_cur_id_.at(bin_coords.first)++) = {(bin_coords.second<<1)+template_segment, length};
		}
		void AddLowQSiteStart(uintRefSeqId ref_seq_id, uintSeqLen position, const Reference &reference, ThreadData &thread){
			auto bin_coords = GetRefSeqBin(ref_seq_id, position, reference, thread);
			lowq_site_start_by_ref_seq_bin_.at(bin_coords.first).at(lowq_site_start_by_ref_seq_bin_cur_id_.at(bin_coords.first)++) = bin_coords.second;
		}
		void AddLowQSiteEnd(uintRefSeqId ref_seq_id, uintSeqLen position, const Reference &reference, ThreadData &thread){
			auto bin_coords = GetRefSeqBin(ref_seq_id, position, reference, thread);
			lowq_site_end_by_ref_seq_bin_.at(bin_coords.first).at(lowq_site_end_by_ref_seq_bin_cur_id_.at(bin_coords.first)++) = bin_coords.second;
			if(bin_coords.second < tmp_insert_lengths_.size()-1 && bin_coords.first && ref_seq_bin_def_.at(bin_coords.first).first == ref_seq_bin_def_.at(bin_coords.first-1).first ){
				// Add position also to the previous bin as fragments starting there potentially end on this position
				lowq_site_end_by_ref_seq_bin_.at(bin_coords.first-1).at(lowq_site_end_by_ref_seq_bin_cur_id_.at(bin_coords.first-1)++) = bin_coords.second + RefSeqSplitLength(ref_seq_id, reference);
			}
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
		uintRefSeqBin GetRefSeqBin(uintRefSeqId ref_seq_id, uintSeqLen position, uintRefSeqBin last_bin){
			uintRefSeqBin new_bin(last_bin);
			if(last_bin < ref_seq_start_bin_.at(ref_seq_id)){
				new_bin = ref_seq_start_bin_.at(ref_seq_id);
			}
			else if(last_bin >= ref_seq_start_bin_.at(ref_seq_id+1)){
				new_bin = ref_seq_start_bin_.at(ref_seq_id+1) - 1;
			}

			while(new_bin+1 < ref_seq_start_bin_.at(ref_seq_id+1) && ref_seq_bin_def_.at(new_bin+1).second <= position){
				++new_bin;
			}
			while(new_bin > ref_seq_start_bin_.at(ref_seq_id) && ref_seq_bin_def_.at(new_bin).second > position){
				--new_bin;
			}

			return new_bin;
		}

		std::pair<uintRefSeqBin, uintSeqLen> GetRefSeqBin(uintRefSeqId ref_seq_id, uintSeqLen position, const Reference &reference, ThreadData &thread){
			auto new_bin = GetRefSeqBin(ref_seq_id, position, thread.last_bin_);

			if(thread.last_bin_ < new_bin){
				if(ref_seq_bin_def_.at(new_bin).first != ref_seq_bin_def_.at(thread.last_bin_).first ){
					thread.last_added_region_ = 0;
					thread.correction_ = (new_bin - ref_seq_start_bin_.at(ref_seq_id)) * RefSeqSplitLength(ref_seq_id, reference) + reference.StartExclusion(ref_seq_id);
				}
				else{
					// Start positions at start of this bin and not the last
					thread.correction_ += (new_bin - thread.last_bin_) * RefSeqSplitLength(ref_seq_id, reference);
				}
			}
			else if(thread.last_bin_ > new_bin){
				if(ref_seq_bin_def_.at(new_bin).first != ref_seq_bin_def_.at(thread.last_bin_).first ){
					thread.last_added_region_ = reference.NumExcludedRegions(ref_seq_id)-1;
					thread.correction_ = (new_bin - ref_seq_start_bin_.at(ref_seq_id)) * RefSeqSplitLength(ref_seq_id, reference) + reference.SumExcludedBases(ref_seq_id);
				}
				else{
					// Start positions at start of this bin and not the last
					thread.correction_ -= (thread.last_bin_ - new_bin) * RefSeqSplitLength(ref_seq_id, reference);
				}
			}
			thread.last_bin_ = new_bin;

			reference.CorrectPositionRemovingExcludedRegions(thread.correction_, thread.last_added_region_, ref_seq_id, position);

			return {new_bin, position-thread.correction_};
		}
		uintSeqLen MaxRefSeqBinLength(const Reference &reference){
			uintSeqLen max(0);
			for(uintRefSeqId seq=0; seq < reference.NumberSequences(); ++seq){
				if(!reference.ReferenceSequenceExcluded(seq)){
					utilities::SetToMax(max, RefSeqSplitLength(seq, reference) + ref_seq_start_bin_.at(seq+1)-ref_seq_start_bin_.at(seq)); // Due to rounding the last bin of a reference sequence can have more than RefSeqSplitLength bases. The number of bin is an overestimate for this, but this little increase compared to RefSeqSplitLength doesn't justify more thinking for a better upper limit
				}
			}
			return max;
		}

		void Prepare( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const std::vector<uintFragCount> &reads_per_frag_len_bin, const std::vector<uintFragCount> &lowq_reads_per_frag_len_bin );
		void FillInOutskirtContent( const Reference &reference, const seqan::BamAlignmentRecord &record_start, uintSeqLen fragment_start_pos, uintSeqLen fragment_end_pos );

		void HandleReferenceSequencesUntil(uintRefSeqId still_needed_reference_sequence, uintSeqLen still_needed_position, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);
		void FinishThreads(ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, std::mutex &print_mutex);

		void Finalize();
		bool FinalizeBiasCalculation(const Reference &reference, uintNumThreads num_threads, FragmentDuplicationStats &duplications);
		double CorrectedCoverage(const Reference &ref, uintReadLen average_read_len);
		bool UpdateRefSeqBias(RefSeqBiasSimulation model, const std::string &bias_file, const Reference &ref, std::mt19937_64 &rgen);
		double CalculateBiasNormalization(std::vector<uintRefSeqId> &coverage_groups, std::vector<std::vector<double>> &non_zero_thresholds, const Reference &reference, uintNumThreads num_threads, uintFragCount total_reads) const;

		static uintDupCount NegativeBinomial(double p, double r, double probability_chosen);
		double Dispersion(double bias) const{ return BiasCalculationVectors::GetDispersion( bias, dispersion_parameters_.at(0), dispersion_parameters_.at(1) ); }
		uintDupCount GetFragmentCounts(double bias_normalization, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintPercent gc, const Surrounding &fragment_start, const Surrounding &fragment_end, double probability_chosen) const;

		void PreparePlotting();
	};
}

#endif // FRAGMENTDISTRIBUTIONSTATS_H
