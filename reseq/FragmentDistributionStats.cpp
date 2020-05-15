#include "FragmentDistributionStats.h"
using reseq::FragmentDistributionStats;
using reseq::BiasCalculationVectors;
using reseq::InsertLengthSpline;

//include <algorithm>
using std::max;
using std::max_element;
using std::min_element;
using std::min;
using std::sort;
//include <array>
using std::array;
//include <atomic>
using std::atomic;
#include <chrono>
#include <exception>
using std::exception;
#include <fstream>
using std::ifstream;
using std::getline;
using std::ofstream;
#include <functional> // std::ref, std::cref
#include <iterator>
using std::distance;
#include <limits>
using std::numeric_limits;
#include <cmath>
using std::abs;
using std::exp;
using std::isnan;
using std::log;
using std::pow;
using std::nan;
using std::sqrt;
//include <mutex>
using std::lock_guard;
using std::mutex;
//include <random>
using std::mt19937_64;
using std::uniform_int_distribution;
#include <set>
using std::set;
#include <string>
using std::stod;
using std::string;
#include <sstream>
using std::stringstream;
#include <thread>
using std::thread;
#include <unordered_map>
using std::unordered_map;
#include <utility>
using std::pair;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::Dna5;

//include "utilities.hpp"
using reseq::utilities::Complement;
using reseq::utilities::at;
using reseq::utilities::Divide;
using reseq::utilities::DivideAndCeil;
using reseq::utilities::IsN;
using reseq::utilities::Percent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;
using reseq::utilities::Sign;
using reseq::utilities::InvLogit2;
using reseq::utilities::VectorAtomic;

constexpr double BiasCalculationVectors::kLowerBound;
constexpr double BiasCalculationVectors::kUpperBound;
constexpr double BiasCalculationVectors::kBaseValue;
mutex BiasCalculationVectors::file_mutex_;

void BiasCalculationVectors::AddCountsFromSite(const FragmentSite &site, array<uintFragCount, 101> &gc_count, array<uintFragCount, 4*Surrounding::Length()> &sur_count){
	if(0 < site.count_forward_ + site.count_reverse_){
		gc_count.at(site.gc_) += site.count_forward_ + site.count_reverse_;

		for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
			sur_count.at(4*sur_pos + site.start_surrounding_.BaseAt(sur_pos)) += site.count_forward_ + site.count_reverse_;
			sur_count.at(4*sur_pos + site.end_surrounding_.BaseAt(sur_pos)) += site.count_forward_ + site.count_reverse_;
		}
	}
}

void BiasCalculationVectors::GetCounts(){
	gc_count_.fill(0);
	sur_count_.fill(0);

	for(auto &site : sites_){
		if(0.0 != site.bias_){
			AddCountsFromSite(site, gc_count_, sur_count_);
		}
	}

	total_counts_ = SumVect(gc_count_);
}

void BiasCalculationVectors::RemoveUnnecessarySites(){
	gc_sites_.fill(0);
	sur_sites_.fill(0);

	uintSeqLen needed_sites(0);
	bool zero_site;
	for(uintSeqLen cur_site = 0; cur_site < sites_.size(); ++cur_site){
		if(0.0 != sites_.at(cur_site).bias_){ // Filter out sides where at least one surrounding couldn't be determined
			zero_site = false; // sites with zero for gc are ok as we use splines there

			for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
				if( !sur_count_.at(4*sur_pos + sites_.at(cur_site).start_surrounding_.BaseAt(sur_pos)) ){
					zero_site = true;
					break;
				}
				if( !sur_count_.at(4*sur_pos + sites_.at(cur_site).end_surrounding_.BaseAt(sur_pos)) ){
					zero_site = true;
					break;
				}
			}

			if(!zero_site){
				sites_.at(needed_sites++) = sites_.at(cur_site);

				gc_sites_.at(sites_.at(cur_site).gc_) += 2;

				for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
					sur_sites_.at(4*sur_pos + sites_.at(cur_site).start_surrounding_.BaseAt(sur_pos)) += 2;
					sur_sites_.at(4*sur_pos + sites_.at(cur_site).end_surrounding_.BaseAt(sur_pos)) += 2;
				}
			}
		}
	}
	sites_.resize(needed_sites);

	total_sites_ = 2*needed_sites; // Forward + Reverse
}

void BiasCalculationVectors::CalculateGCWeights(){
	gc_weights_.fill(0.0);
	for(uintPercent gc=0; gc < gc_sites_.size(); ++gc){
		if(1 < gc_sites_.at(gc)){
			double center_weight = sqrt(gc_sites_.at(gc))-1;
			gc_weights_.at(gc) += center_weight;

			// Reducing weight by kDistanceWeightReductionFactor every gc percent we go away from the sites
			double weight = center_weight;
			for(auto low_gc = gc; low_gc--; ){
				weight *= kDistanceWeightReductionFactor;
				gc_weights_.at(low_gc) += weight;
			}

			weight = center_weight;
			for(auto high_gc = gc+1; high_gc < gc_sites_.size(); ++high_gc){
				weight *= kDistanceWeightReductionFactor;
				gc_weights_.at(high_gc) += weight;
			}
		}
	}

	// Normalize total weight to 1
	double weight_sum(0.0);
	for(uintPercent gc=0; gc < gc_weights_.size(); ++gc){
		weight_sum += gc_weights_.at(gc);
	}
	for(uintPercent gc=0; gc < gc_weights_.size(); ++gc){
		gc_weights_.at(gc) /= weight_sum;
	}
}

void BiasCalculationVectors::AddBaseValue(uintDupCount count){
	if(1 < count){
		double k_fac = 2;
		for(uintDupCount k = 3; k <= count; ++k){
			k_fac *= k;
		}
		loglike_pois_base_ -= log(k_fac);
	}
}

void BiasCalculationVectors::GetLogLikeBase(){
	loglike_pois_base_ = 0;
	for(auto &site : sites_){
		AddBaseValue(site.count_forward_);
		AddBaseValue(site.count_reverse_);
	}
}

void BiasCalculationVectors::DeactivateZeroCounts(vector<double> &par, double deactive_value, uintNumFits sur_shift){
	for( uintNumFits sur=sur_count_.size(); sur--; ){
		if(!sur_count_.at(sur)){
			par.at(sur+sur_shift) = deactive_value;
		}
	}
}

void BiasCalculationVectors::NormGC(){
	// Sort gc by number of sites
	array<pair<uintFragCount, uintPercent>, 101> sites;
	for(uintPercent gc = gc_sites_.size(); gc--; ){
		sites.at(gc) = {gc_sites_.at(gc), gc};
	}
	sort(sites.begin(), sites.end());

	// Sum the biases up to kPercentGCSitesForNormalization of the sites starting with the gc with most sites
	uintFragCount sum_sites = 0;
	double sum_bias = 0.0;
	uintPercent num_bias = 0;
	for(uintPercent gci = gc_sites_.size(); gci-- && sum_sites < Divide(total_sites_*kPercentGCSitesForNormalization,static_cast<uintFragCount>(100)); ){
		sum_sites += sites.at(gci).first;
		sum_bias += gc_bias_.at( sites.at(gci).second );
		++num_bias;
	}

	// Normalize using the biases of the gc values having kPercentGCSitesForNormalization of the sites to take only well determined values
	for(uintPercent gc = gc_sites_.size(); gc--; ){
		gc_bias_.at(gc) *= num_bias/sum_bias;
	}
}

void BiasCalculationVectors::NormSurroundings(const vector<double> &x){
	if(kSurMult){
		for(auto sur_pos = 0; sur_pos < Surrounding::Length(); ++sur_pos){
			double sur_sum = 0.0;
			uintBaseCall valid_sur = 0;
			for(auto base = 0; base < 4; ++base){
				sur_sum += x.at(4*sur_pos+base);
				if(sur_count_.at(4*sur_pos+base)){
					++valid_sur;
				}
			}
			for(auto base = 0; base < 4; ++base){
				sur_bias_.at(4*sur_pos+base) = valid_sur*x.at(4*sur_pos+base)/sur_sum;
			}
		}
	}
	else{
		for(auto sur_pos = 0; sur_pos < Surrounding::Length(); ++sur_pos){
			auto from = 1+sur_pos*4;
			auto to = 1+(sur_pos+1)*4;

			double sur_sum(0.0);
			uintBaseCall valid_sur = 0;
			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur-1)){
					sur_sum += x.at(sur);
					++valid_sur;
				}
			}
			sur_sum /= valid_sur;

			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur-1)){
					sur_bias_.at(sur-1) = x.at(sur) - sur_sum;
				}
			}
		}

		for(auto sur = 0; sur < 4; ++sur){
			sur_bias_.at(sur) += x.at(0);
		}
	}
}

void BiasCalculationVectors::UnnormSurroundingGradients(vector<double> &grad, const vector<double> &x){
	if(kSurMult){
		for(auto sur_pos = 0; sur_pos < Surrounding::Length(); ++sur_pos){
			auto from = sur_pos*4;
			auto to = (sur_pos+1)*4;

			double grad_sum(0.0), bias_sum(0.0);
			uintBaseCall valid_sur = 0;
			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur)){
					grad_sum += sur_bias_.at(sur) * sur_grad_.at(sur);
					bias_sum += x.at(sur);
					++valid_sur;
				}
			}

			for(auto sur = from; sur < to; ++sur){
				grad.at(sur) = (valid_sur*sur_grad_.at(sur)-grad_sum)/bias_sum;
			}
		}
	}
	else{
		for(auto sur_pos = 0; sur_pos < Surrounding::Length(); ++sur_pos){
			auto from = sur_pos*4;
			auto to = (sur_pos+1)*4;

			double grad_mean(0.0);
			uintBaseCall valid_sur = 0;
			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur)){
					grad_mean += sur_grad_.at(sur);
					++valid_sur;
				}
			}
			grad_mean /= valid_sur;

			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur)){
					grad.at(sur+1) = sur_grad_.at(sur)-grad_mean;
				}
				else{
					grad.at(sur+1) = 0.0;
				}
			}
		}

		grad.at(0) = 0.0;
		for(uintBaseCall sur = 0; sur < 4; ++sur){
			grad.at(0) += sur_grad_.at(sur);
		}
	}
}

double BiasCalculationVectors::SurBiasAtSite(const pair<double, double>& bias){
	if(kSurMult){
		return static_cast<double>(total_counts_) / total_sites_ * bias.first * bias.second;
	}
	else{
		return static_cast<double>(total_counts_) / total_sites_ * bias.first * bias.second + kBaseValue;
	}
}

pair<double, double> BiasCalculationVectors::SurBiasAtSiteSplit(const FragmentSite& site){
	pair<double, double> bias;
	if(kSurMult){
		bias = {1.0, 1.0};
	}
	else{
		bias = {0.0, 0.0};
	}

	for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
		if(kSurMult){
			bias.first *= sur_bias_.at(4*sur_pos + site.start_surrounding_.BaseAt(sur_pos));
			bias.second *= sur_bias_.at(4*sur_pos + site.end_surrounding_.BaseAt(sur_pos));
		}
		else{
			bias.first += sur_bias_.at(4*sur_pos + site.start_surrounding_.BaseAt(sur_pos));
			bias.second += sur_bias_.at(4*sur_pos + site.end_surrounding_.BaseAt(sur_pos));
		}
	}

	if(!kSurMult){
		bias.first = InvLogit2(bias.first);
		bias.second = InvLogit2(bias.second);
	}

	return bias;
}

double BiasCalculationVectors::SurBiasAtSite(const FragmentSite& site){
	auto bias = SurBiasAtSiteSplit(site);
	return SurBiasAtSite(bias);
}

void BiasCalculationVectors::DefineStartingKnots(){
	// Define knots for gc
	uintFragCount sum(0);
	uintPercent k(0);
	for(uintPercent gc=0; gc < gc_sites_.size(); ++gc){
		if( gc_count_.at(gc) ){
			// Get an equal spacing of sites (First site where the cumulated sum is at least the k-th quartile)
			sum += gc_sites_.at(gc);
			gc_knots_.at(k) = gc;
			if(sum >= Divide(k*total_sites_,static_cast<uintFragCount>(kGCSplineDf)-1)){
				if( ++k < gc_knots_.size() ){
					gc_knots_.at(k) = 102; // Set it invalid
				}
			}
		}
	}

	if( k < gc_knots_.size() && 102 > gc_knots_.at(k)){
		++k; // gc_knots_.at(k) is the last gc with a count, so don't overwrite that one in the filling procedure coming next
	}

	// In case there are not enough gc bins with counts
	if( k < gc_knots_.size() ){
		// Fill rest of array with invalid gc, to not jump over an unfilled bin that accidentally happens to have the right gc
		for(auto k_fill=k; k_fill<gc_knots_.size(); ++k_fill){
			gc_knots_.at(k_fill) = 102;
		}

		// Start from the first knot and assign unused gc bins with zero counts
		uintPercent next_k = 1;
		uintPercent gc = gc_knots_.at(0) + 1;
		while(k < gc_knots_.size() && gc < 101){
			if(gc == gc_knots_.at(next_k)){
				// GC bin is already a knot, go to next
				++gc;
				++next_k;
			}
			else{
				gc_knots_.at(k++) = gc++;
			}
		}

		// In case we reached GC 101
		gc = gc_knots_.at(0);
		while(k < gc_knots_.size()){
			// Has to be lower than the lowest knot as we start from that one
			gc_knots_.at(k++) = --gc;
		}

		sort(gc_knots_.begin(), gc_knots_.end());
	}
}

template<typename T, typename S, typename U, typename V, typename W> void BiasCalculationVectors::PrepareSplines(array<T, 3> &lin_comb_splines, S &x, W &l, U &c, U &z, V &mu, V &h, T &beta){
	// https://en.wikipedia.org/wiki/Spline_(mathematics) "Computation of Natural Cubic Splines" slightly adapted to get the linear combination of the spline parameters a
	auto &b = lin_comb_splines.at(0); // b[i][a_index]: i=[0, n-1]
	auto &d = lin_comb_splines.at(2); // d[i][a_index]: i=[0, n-1]

	for(uintPercent i=0; i < h.size(); ++i){
		h.at(i) = x.at(i+1) - x.at(i);
	}

	for(uintPercent k=1; k < beta.size(); ++k){
		//beta.at(k).fill(0.0);
		beta.at(k).at(k+1) = 3/h.at(k);
		beta.at(k).at(k) = -3/h.at(k) - 3/h.at(k-1);
		beta.at(k).at(k-1) = 3/h.at(k-1);
	}

	l.at(0) = 0.0;
	mu.at(0) = 0.0;
	//z.at(0).fill(0.0);

	for(uintPercent k=1; k < beta.size(); ++k){
		l.at(k) = 2*(x.at(k+1)-x.at(k-1)) - h.at(k-1)*mu.at(k-1);
		mu.at(k) = h.at(k)/l.at(k);
		for(uintPercent ai=0; ai<x.size(); ++ai){
			z.at(k).at(ai) = (beta.at(k).at(ai) - h.at(k-1) * z.at(k-1).at(ai))/l.at(k);
		}
	}

	l.at(l.size()-1) = 1.0;
	//z.at(z.size()-1).fill(0.0);
	//c.at(c.size()-1).fill(0.0);

	for(uintPercent i=h.size(); i--;){
		for(uintPercent ai=0; ai<x.size(); ++ai){
			c.at(i).at(ai) = z.at(i).at(ai) - mu.at(i)*c.at(i+1).at(ai);
			b.at(i).at(ai) = -h.at(i)*(c.at(i+1).at(ai) + 2*c.at(i).at(ai))/3;
			d.at(i).at(ai) = (c.at(i+1).at(ai) - c.at(i).at(ai))/3/h.at(i);
		}

		b.at(i).at(i+1) += 1/h.at(i);
		b.at(i).at(i) -= 1/h.at(i);
	}

	// Store values from c in lin_comb_gc_splines_ (as c needs to be one larger for the calculation it is not a reference like b and d)
	for(uintPercent i=h.size(); i--;){
		for(uintPercent ai=0; ai<x.size(); ++ai){
			lin_comb_splines.at(1).at(i).at(ai) = c.at(i).at(ai);
		}
	}
}

void BiasCalculationVectors::PrepareSplines(){
	array<double, kGCSplineDf> l; // j=[0, n]
	array<array<double, kGCSplineDf>, kGCSplineDf> c, z; // c[j][a_index]: j=[0, n]
	array<double, kGCSplineDf-1> mu, h; // i=[0, n-1]
	array<array<double, kGCSplineDf>, kGCSplineDf-1> beta; // beta[k][a_index]: k=[1, n-1]

	for(uintPercent k=1; k < beta.size(); ++k){
		beta.at(k).fill(0.0);
	}
	z.at(0).fill(0.0);
	z.at(z.size()-1).fill(0.0);
	c.at(c.size()-1).fill(0.0);

	PrepareSplines(lin_comb_gc_splines_, gc_knots_, l, c, z, mu, h, beta);
}

template<typename T> void BiasCalculationVectors::GetSplineCoefficients(double &a, double &b, double &c, double &d, uintSeqLen k, const vector<double> &spline_pars, const array<T, 3> &lin_comb_splines){
	a = spline_pars.at(k+1);
	b = 0.0;
	c = 0.0;
	d = 0.0;

	for(uintPercent ai=1; ai < spline_pars.size(); ++ai){
		b += spline_pars.at(ai) * lin_comb_splines.at(0).at(k).at(ai-1);
		c += spline_pars.at(ai) * lin_comb_splines.at(1).at(k).at(ai-1);
		d += spline_pars.at(ai) * lin_comb_splines.at(2).at(k).at(ai-1);
	}
}

void BiasCalculationVectors::GetSplineCoefficients(double &a, double &b, double &c, double &d, uintPercent k, const vector<double> &spline_pars){
	GetSplineCoefficients(a, b, c, d, k, spline_pars, lin_comb_gc_splines_);
}

void BiasCalculationVectors::CalculateSpline(const vector<double> &spline_pars){
	// Constant bias values for gc below any observed sites: Extrapolate the linear function at the lowest measured GC by half a GC unit and use this value
	double norm = spline_pars.at(0);
	double a, b, c, d;
	GetSplineCoefficients(a, b, c, d, 0, spline_pars);
	double lower_const_no_logit = a - 0.5*b;
	double lower_const;
	if( kGCExp ){
		lower_const = norm*exp(lower_const_no_logit) + kBaseValue;
	}
	else{
		lower_const = norm*InvLogit2(lower_const_no_logit) + kBaseValue;
	}
	for(uintPercent gc=0; gc < gc_knots_.at(0); ++gc){
		gc_bias_no_logit_.at(gc) = lower_const_no_logit;
		gc_bias_.at(gc) = lower_const;
	}

	// Fill observed range of gc with spline values
	gc_bias_no_logit_.at( gc_knots_.at(0) ) = a;
	if( kGCExp ){
		gc_bias_.at( gc_knots_.at(0) ) = norm*exp(a) + kBaseValue;
	}
	else{
		gc_bias_.at( gc_knots_.at(0) ) = norm*InvLogit2(a) + kBaseValue;
	}

	uintPercent k=0;
	for(uintPercent gc=gc_knots_.at(0)+1; gc < gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		if(gc == gc_knots_.at(k+1)){
			// Reached new spline part
			GetSplineCoefficients(a, b, c, d, ++k, spline_pars);
			gc_bias_no_logit_.at(gc) = a;
			if( kGCExp ){
				gc_bias_.at(gc) = norm*exp(a) + kBaseValue;
			}
			else{
				gc_bias_.at(gc) = norm*InvLogit2(a) + kBaseValue;
			}
		}
		else{
			// Continue on current spline part
			uintPercent cur_gc = gc-gc_knots_.at(k);
			gc_bias_no_logit_.at(gc) = a + b*cur_gc + c*cur_gc*cur_gc + d*cur_gc*cur_gc*cur_gc;
			if( kGCExp ){
				gc_bias_.at(gc) = norm*exp( gc_bias_no_logit_.at(gc) ) + kBaseValue;
			}
			else{
				gc_bias_.at(gc) = norm*InvLogit2( gc_bias_no_logit_.at(gc) ) + kBaseValue;
			}
		}
	}
	++k;

	gc_bias_no_logit_.at( gc_knots_.at(k) ) = spline_pars.at(k+1);
	if( kGCExp ){
		gc_bias_.at( gc_knots_.at(k) ) = norm*exp(spline_pars.at(k+1)) + kBaseValue;
	}
	else{
		gc_bias_.at( gc_knots_.at(k) ) = norm*InvLogit2(spline_pars.at(k+1)) + kBaseValue;
	}

	// Constant bias values for gc above any observed sites: Extrapolate the linear function at the highest measured GC by half a GC unit and use this value
	uintPercent cur_gc = gc_knots_.at(k)-gc_knots_.at(k-1);
	double upper_const_no_logit = spline_pars.at(k+1) + 0.5*(b + c*cur_gc); // d=0 since second derivative must be 0 (natural spline)
	double upper_const;
	if( kGCExp ){
		upper_const = norm*exp( upper_const_no_logit ) + kBaseValue;
	}
	else{
		upper_const = norm*InvLogit2( upper_const_no_logit ) + kBaseValue;
	}
	for(uintPercent gc=gc_knots_.at(k)+1; gc < gc_bias_.size(); ++gc){
		gc_bias_no_logit_.at(gc) = upper_const_no_logit;
		gc_bias_.at(gc) = upper_const;
	}
}

void BiasCalculationVectors::CalculateSplineGrad(const vector<double> &spline_pars, vector<double> &grad, uintNumFits grad_offset){
	for( uintPercent gc=0; gc < gc_bias_grad_.size(); ++gc ){
		grad.at(grad_offset) += gc_bias_grad_.at(gc) / spline_pars.at(0);
		if( !kGCExp ){
			gc_bias_grad_.at(gc) *= InvLogit2(-gc_bias_no_logit_.at(gc))/2;
		}
	}

	// Constant bias values for gc below any observed sites
	for( uintPercent gc=0; gc < gc_knots_.at(0); ++gc ){
		grad.at(grad_offset+1) += gc_bias_grad_.at(gc);

		for(uintPercent ai=0; ai < gc_spline_pars_.size()-1; ++ai){
			grad.at(grad_offset+ai+1) -= 0.5*lin_comb_gc_splines_.at(0).at(0).at(ai) * gc_bias_grad_.at(gc);
		}
	}

	// Observed range of gc
	uintPercentShift k=-1;
	for( uintPercent gc=gc_knots_.at(0); gc <= gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		if(gc == gc_knots_.at(k+1)){
			// Reached new spline part
			++k;
		}
		else{
			// Continue on current spline part
			uintPercent cur_gc = gc-gc_knots_.at(k);
			for(uintPercent ai=0; ai < gc_spline_pars_.size()-1; ++ai){
				grad.at(grad_offset+ai+1) += (lin_comb_gc_splines_.at(0).at(k).at(ai)*cur_gc + lin_comb_gc_splines_.at(1).at(k).at(ai)*cur_gc*cur_gc + lin_comb_gc_splines_.at(2).at(k).at(ai)*cur_gc*cur_gc*cur_gc) * gc_bias_grad_.at(gc);
			}
		}
		grad.at(grad_offset+k+1) += gc_bias_grad_.at(gc);
	}

	// Constant bias values for gc above any observed sites
	uintPercent cur_gc = gc_knots_.at(k)-gc_knots_.at(k-1);
	for( uintPercent gc=gc_knots_.at( gc_knots_.size()-1 )+1; gc < gc_count_.size(); ++gc ){
		grad.at(grad_offset+k+1) += gc_bias_grad_.at(gc);
		for(uintPercent ai=0; ai < gc_spline_pars_.size()-1; ++ai){
			grad.at(grad_offset+ai+1) += 0.5*(lin_comb_gc_splines_.at(0).at(k-1).at(ai) + lin_comb_gc_splines_.at(1).at(k-1).at(ai)*cur_gc) * gc_bias_grad_.at(gc);
		}
	}
}

double BiasCalculationVectors::LogLikeGcSpline(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_;

	calc.CalculateSpline(x);

	double loglike(calc.loglike_pois_base_);
	for( uintPercent gc=0; gc < calc.gc_count_.size(); ++gc){
		loglike += calc.gc_bias_sum_.at(gc).first + calc.gc_count_.at(gc)*log(calc.gc_bias_.at(gc)) - calc.gc_bias_.at(gc)*calc.gc_bias_sum_.at(gc).second;
	}

	if( grad.size() ){
		for( auto &g : grad ){
			g = 0.0;
		}

		for( uintPercent gc=0; gc < calc.gc_count_.size(); ++gc){
			//double grad_gc = gc_count_.at(gc)/gc_bias_.at(gc) - gc_bias_sum_.at(gc).second;
			calc.gc_bias_grad_.at(gc) = (calc.gc_count_.at(gc) - calc.gc_bias_sum_.at(gc).second*calc.gc_bias_.at(gc));
		}

		calc.CalculateSplineGrad(x, grad, 0);
	}

	return loglike;
}

double BiasCalculationVectors::OptimizeSpline(){
	PrepareSplines();

	double max_val = 0.0; // This will be our normalization variable
	for( uintPercent k=gc_spline_pars_.size()-1; k--; ){
		SetToMax(max_val, gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second);
	}
	if(!kGCExp){
		max_val /= 1.6; // The maximum normalized value should be 1.6, so that we stay in the mostly linear part of the 2*InvLogit
	}
	gc_spline_pars_.at(0) = max_val;

	for( uintPercent k=gc_spline_pars_.size()-1; k--; ){
		if( gc_count_.at( gc_knots_.at(k) ) ){
			double val = gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second;

			// Set gc bias to set the gradient to zero for this gc
			if(kGCExp){
				gc_spline_pars_.at(k+1) = log(val/max_val);
			}
			else{
				gc_spline_pars_.at(k+1) = log(val/(2*max_val-val));
			}

			if(kUpperBound < gc_spline_pars_.at(k+1)){
				gc_spline_pars_.at(k+1) = kUpperBound;
			}
			else if(kLowerBound > gc_spline_pars_.at(k+1)){
				gc_spline_pars_.at(k+1) = kLowerBound;
			}
		}
		else{
			gc_spline_pars_.at(k+1) = kLowerBound;
		}
	}

	func_calls_ = 0;
	double loglike;

	converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == optimizer_spline_.optimize(gc_spline_pars_, loglike) ){
			func_calls_ += 10000;
			converged_ = false;
		}
	}
	catch(const std::exception& e){
		func_calls_ += 20000;
		converged_ = false;
	}

	return loglike;
}

double BiasCalculationVectors::OptimizeSpline(uintPercent knot, uintPercentShift shift){
	// knot is assumed to be an inner knot
	if( (shift < 0 && gc_knots_.at(knot-1) >= gc_knots_.at(knot)+shift) || (shift > 0 && gc_knots_.at(knot+1) <= gc_knots_.at(knot)+shift) || (0 == gc_count_.at( gc_knots_.at(knot) + shift )) ){
		return -1 * numeric_limits<double>::infinity();
	}
	else{
		gc_knots_.at(knot) += shift;
		double loglike = OptimizeSpline();
		gc_knots_.at(knot) -= shift;
		return loglike;
	}
}

reseq::uintPercent BiasCalculationVectors::OptimizeSplineIdToKnot(uintNumFits id, const uintPercent range){
	return (id+2*range-1)/(2*range);
}

reseq::uintPercentShift BiasCalculationVectors::OptimizeSplineIdToShift(uintNumFits id, const uintPercent range){
	return static_cast<uintPercentShift>((id-1)%(2*range)) - range + ((id-1)%(2*range) < range?0:1);
}

void BiasCalculationVectors::GetGCSpline(){
	DefineStartingKnots();

	NormSurroundings(fit_pars_);

	// Get sums over all sites needed for calculation
	gc_bias_sum_.fill({0.0,0.0});
	for(auto &site : sites_){
		auto bias = SurBiasAtSite(site);

		if(site.count_forward_+site.count_reverse_){
			gc_bias_sum_.at(site.gc_).first += (site.count_forward_+site.count_reverse_)*log(bias);
		}
		gc_bias_sum_.at(site.gc_).second += bias; // Multiplication with 2 later after sum is complete
	}

	for( uintPercent gc=gc_knots_.at(0); gc <= gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		gc_bias_sum_.at(gc).second *= 2; // Multiplication with 2 here to save one multiplication per site
	}

	if(kParameterInfoFile || kDispersionInfoFile){
		// Get "perfect" GC bias by setting gradient to zero for each gc bin
		for( uintPercent gc=0; gc < gc_count_.size(); ++gc){
			if(gc_count_.at(gc)){
				gc_bias_.at(gc) = gc_count_.at(gc)/gc_bias_sum_.at(gc).second;
			}
			else{
				gc_bias_.at(gc) = 0;
			}
		}

		if(kDispersionInfoFile){
			WriteOutDispersionInfo("poisson");
		}

		if(BiasCalculationVectors::kParameterInfoFile){
			WriteOutParameterInfo("poisson");
		}
	}

	// Optimize gc spline
	// Knots are adjusted by comparing the likelihoods to the ones with all inner knots in turn moved by 1 to range up or down and taking the best likelihood in a greedy fashion
	// Iterate until likelihood cannot be improved by moving an inner knot by 1 to range gc
	array<double, (kGCSplineDf-2)*2*kMaxKnotShift+1 > loglikes;
	uintNumFits bestlike(1);
	loglikes.at(0) = OptimizeSpline();
	loglikes.at(1) = OptimizeSpline(1, -static_cast<uintPercentShift>(kMaxKnotShift));
	while(bestlike != 0){ // loglikes[0] is current knot position
		// Calculate loglikelihoods with inner knots in turn moved by one up or down
		for(uintNumFits cur_like=loglikes.size(); --cur_like; ){ // Leave out 0
			if(cur_like != bestlike){ // Otherwise we still have the value from last iteration
				loglikes.at(cur_like) = OptimizeSpline( OptimizeSplineIdToKnot(cur_like, kMaxKnotShift), OptimizeSplineIdToShift(cur_like, kMaxKnotShift) );
			}
		}

		bestlike = distance( loglikes.begin(), max_element(loglikes.begin(), loglikes.end()) );
		if(bestlike != 0){
			// Save the two likelihoods that will be identical in the next round in the proper spots
			uintNumFits flipped_bestlike = bestlike - 2*(static_cast<uintPercentShift>((bestlike-1)%(2*kMaxKnotShift)) - kMaxKnotShift) - 1;
			loglikes.at(flipped_bestlike) = loglikes.at(0);
			loglikes.at(0) = loglikes.at(bestlike);
			// Apply the knot movement
			gc_knots_.at( OptimizeSplineIdToKnot(bestlike, kMaxKnotShift) ) += OptimizeSplineIdToShift(bestlike, kMaxKnotShift);
			// Store at which new position the original likelihood is, so that we don't recalculate it
			bestlike = flipped_bestlike;
		}
	}

	// Get gc bias for best spline (We require a higher precision here as the surrounding gradient assumes that the gc gradients are exactly zero, so we don't want to deviate much)
	loglike_ = OptimizeSpline();

	if(kParameterInfoFile || kDispersionInfoFile){
		CalculateSpline(gc_spline_pars_);

		if(kDispersionInfoFile){
			WriteOutDispersionInfo("gcspline");
		}

		if(BiasCalculationVectors::kParameterInfoFile){
			WriteOutParameterInfo("gcspline");
		}
	}
}

double BiasCalculationVectors::LogLikelihoodPoisson(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_;

	calc.NormSurroundings(x);

	if( grad.size() ){
		if(BiasCalculationVectors::kSurMult){
			for(uintNumFits sur=calc.sur_count_.size(); sur--;){
				calc.sur_grad_.at(sur) = calc.sur_count_.at(sur)/calc.sur_bias_.at(sur);
			}
		}
		else{
			calc.sur_grad_.fill(0.0);
		}
	}

	// Get sums over all sites needed for calculation
	calc.gc_bias_sum_.fill({0.0,0.0});
	for( auto gc=calc.grad_gc_bias_sum_.size(); gc--; ){
		calc.grad_gc_bias_sum_.at(gc).fill({0.0, 0.0});
	}

	for(auto &site : calc.sites_){
		auto bias_split = calc.SurBiasAtSiteSplit(site);
		auto bias = calc.SurBiasAtSite(bias_split);

		if(site.count_forward_+site.count_reverse_){
			calc.gc_bias_sum_.at(site.gc_).first += (site.count_forward_+site.count_reverse_)*log(bias);
		}
		calc.gc_bias_sum_.at(site.gc_).second += bias; // Multiplication with 2 later after sum is complete

		if( grad.size() ){
			auto cur_grad_start = 1-bias_split.first/2;
			auto cur_grad_end = 1-bias_split.second/2;

			for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
				auto i_start = 4*sur_pos + site.start_surrounding_.BaseAt(sur_pos);
				auto i_end = 4*sur_pos + site.end_surrounding_.BaseAt(sur_pos);

				if(BiasCalculationVectors::kSurMult){
					calc.grad_gc_bias_sum_.at(site.gc_).at(i_start).first += bias/calc.sur_bias_.at(i_start);
					calc.grad_gc_bias_sum_.at(site.gc_).at(i_end).first += bias/calc.sur_bias_.at(i_end);
				}
				else{
					calc.grad_gc_bias_sum_.at(site.gc_).at(i_start).first += bias*cur_grad_start;
					calc.grad_gc_bias_sum_.at(site.gc_).at(i_end).first += bias*cur_grad_end;

					if(site.count_forward_+site.count_reverse_){
						calc.sur_grad_.at(i_start) += (site.count_forward_+site.count_reverse_)*cur_grad_start;
						calc.sur_grad_.at(i_end) += (site.count_forward_+site.count_reverse_)*cur_grad_end;
					}
				}
			}
		}
	}

	for( uintPercent gc=calc.gc_count_.size(); gc--; ){
		calc.gc_bias_sum_.at(gc).second *= 2; // Multiplication with 2 here to save one multiplication per site
	}

	// Calculate gc bias by setting gradient to zero for every gc bin
	double loglike = calc.loglike_pois_base_;
	for( uintPercent gc=calc.gc_count_.size(); gc--; ){
		if(calc.gc_count_.at(gc)){
			calc.gc_bias_.at(gc) = calc.gc_count_.at(gc)/calc.gc_bias_sum_.at(gc).second;
			loglike += calc.gc_bias_sum_.at(gc).first + calc.gc_count_.at(gc)*log(calc.gc_bias_.at(gc)) - calc.gc_count_.at(gc);
		}
		else{
			calc.gc_bias_.at(gc) = 0;
		}
	}

	// Calculate gradients for surroundings based on the gc bias
	if( grad.size() ){
		for( uintPercent gc=calc.gc_count_.size(); gc--; ){
			for(auto sur=calc.sur_bias_.size(); sur--;){
				calc.grad_gc_bias_sum_.at(gc).at(sur).first *= 2; // Multiplication with 2 here to save two multiplication per site
				calc.sur_grad_.at(sur) -= calc.gc_bias_.at(gc) * calc.grad_gc_bias_sum_.at(gc).at(sur).first;
			}
		}

		calc.UnnormSurroundingGradients(grad, x);
	}

	return loglike;
}

double BiasCalculationVectors::GetDispersion(double bias, double a, double b){
	double r = bias/(a+b*bias);
	if(r > bias*1e10){
		// r is so much larger than the bias that r + bias is having precision issues, so shift it back into precision (it anyways is approximately poisson at this point)
		r = bias*1e10;
	}
	return r;
}

void BiasCalculationVectors::LogLikelihoodNbinomSite(double &loglike, double &a, double &b, double &ga, double &gb, vector<double> &grad, const FragmentSite &site){
	auto bias_split = SurBiasAtSiteSplit(site);
	double bias = SurBiasAtSite(bias_split) * gc_bias_.at(site.gc_);
	double r = GetDispersion(bias, a, b);
	double grad_term1;

	grad_term1 = r*(site.count_forward_ + site.count_reverse_ - 2*bias)/(bias+r);

	ga += r * grad_term1/bias;
	gb += r * grad_term1;

	grad_term1 *= (1-a*r/bias);

	if( site.count_forward_ + site.count_reverse_ ){
		loglike += (site.count_forward_ + site.count_reverse_) * log(bias/(bias+r));
	}

	double grad_term2 = 2*log(1/(bias/r+1)); // 2* because of forward+reverse; Be careful to handle r->inf without numerical instabilities, should result in r*grad_term2 -> -2*bias
	loglike += r*grad_term2;

	for(uintDupCount i=1; i <= site.count_forward_; ++i){
		grad_term2 += 1/(r+(i-1)); // Be careful to handle r->0 (brackets around i-1, so that r+1-1=r and not 0)
		loglike += log( (r+(i-1))/i );
	}
	for(uintDupCount i=1; i <= site.count_reverse_; ++i){
		grad_term2 += 1/(r+(i-1));
		loglike += log( (r+(i-1))/i );
	}

	grad_term2 *= r*r;

	ga -= grad_term2/bias;
	gb -= grad_term2;

	grad_term2 *= a/bias;

	grad_term1 += grad_term2;

	gc_bias_grad_.at(site.gc_) += grad_term1;

	if(grad.size()){
		auto cur_grad_start = 1-bias_split.first/2;
		auto cur_grad_end = 1-bias_split.second/2;

		for(uintSurPos sur_pos = Surrounding::Length(); sur_pos--; ){
			if(kSurMult){
				sur_grad_.at(4*sur_pos + site.start_surrounding_.BaseAt(sur_pos)) += grad_term1; // Division by bias later
				sur_grad_.at(4*sur_pos + site.end_surrounding_.BaseAt(sur_pos)) += grad_term1;
			}
			else{
				sur_grad_.at(4*sur_pos + site.start_surrounding_.BaseAt(sur_pos)) += grad_term1*cur_grad_start;
				sur_grad_.at(4*sur_pos + site.end_surrounding_.BaseAt(sur_pos)) += grad_term1*cur_grad_end;
			}
		}
	}
}

double BiasCalculationVectors::LogLikelihoodNbinom(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_;

	calc.NormSurroundings(x);

	uintNumFits spline_par_offset(calc.sur_count_.size());
	if(!BiasCalculationVectors::kSurMult){
		spline_par_offset += 1;
	}

	for(uintPercent k=0; k<calc.gc_spline_pars_.size(); ++k){
		calc.gc_spline_pars_.at(k) = x.at( spline_par_offset + k );
	}
	calc.CalculateSpline(calc.gc_spline_pars_);

	double loglike(0.0), a( x.at(x.size()-2) ), b( x.at(x.size()-1) ), ga(0.0), gb(0.0);
	if( grad.size() ){
		for(auto sur=spline_par_offset; sur--;){
			grad.at(sur) = 0.0;
		}
		for(auto k=spline_par_offset; k<grad.size(); ++k){
			grad.at(k) = 0.0;
		}

		calc.sur_grad_.fill(0.0);
	}

	// Get sums over all sites needed for calculation
	calc.gc_bias_grad_.fill(0.0);
	for(auto &site : calc.sites_){
		calc.LogLikelihoodNbinomSite(loglike, a, b, ga, gb, grad, site);
	}

	if( grad.size() ){
		// Calculate gradient of spline parameter from gradient of individual gc bins
		calc.CalculateSplineGrad(x, grad, spline_par_offset);

		// Calculate gradient of unnormalized surroundings from normalized surroundings
		if(BiasCalculationVectors::kSurMult){
			for(auto sur=calc.sur_count_.size(); sur--;){
				calc.sur_grad_.at(sur) /= calc.sur_bias_.at(sur); // Division by bias here instead of in the site loop
			}
		}
		calc.UnnormSurroundingGradients(grad, x);

		grad.at( grad.size()-2 ) = ga;
		grad.at( grad.size()-1 ) = gb;
	}

	return loglike;
}

void BiasCalculationVectors::PostProcessFit(){
	// Normalize sur biases to 4
	NormSurroundings(fit_pars_);

	// Get correct gc bias values
	for(uintNumFits k=0; k<gc_spline_pars_.size(); ++k){
		gc_spline_pars_.at(k) = fit_pars_.at( sur_count_.size() + (!kSurMult?1:0) + k );
	}
	CalculateSpline(gc_spline_pars_);

	if(kDispersionInfoFile){
		WriteOutDispersionInfo("nbinom");
	}

	dispersion_.resize(2);
	dispersion_.at(0) = fit_pars_.at( fit_pars_.size()-2 );
	dispersion_.at(1) = fit_pars_.at( fit_pars_.size()-1 );

	if(BiasCalculationVectors::kParameterInfoFile){
		WriteOutParameterInfo("nbinom");
	}

	NormGC();
}


double BiasCalculationVectors::LogLikelihoodConstDispersion(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_const_disp_;
	double loglike(0.0), gr(0.0), gr0(0.0), r(x.at(0)), bias; // Calculate log likelihood for negative binomial with constant dispersion

	for( auto site = calc.sites_.begin() + calc.start_site_dispersion_fit_; site != calc.sites_.begin() + calc.end_site_dispersion_fit_; ++site){
		if( isnan(calc.use_sample_mean_)){
			bias = site->bias_;
		}
		else{
			bias = calc.use_sample_mean_;
		}
		gr0 += (2*bias - site->count_forward_ - site->count_reverse_ )/(bias+r);
		gr += log(r/(bias+r));
		loglike += (site->count_forward_ + site->count_reverse_) * log(bias/(bias+r));
	}

	gr *= 2;
	loglike += r*gr;
	gr += gr0;
	for(uintDupCount i=1; i < calc.duplication_count_part_.size(); ++i){
		gr += calc.duplication_count_part_.at(i)/(r+i-1);
		loglike += calc.duplication_count_part_.at(i) * log( (r+i-1)/i );
	}

	if( grad.size() ){
		grad.at(0) = gr;
	}

	return loglike;
}

void BiasCalculationVectors::WriteOutParameterInfo(const char *context){
	stringstream output;

	std::array<double, 101> unnormalized_gc;
	for(auto gc=unnormalized_gc.size(); gc--; ){
		unnormalized_gc.at(gc) = gc_bias_.at(gc);
	}
	NormGC();

	output << context << ", " << ref_seq_bin_ << ", " << insert_length_ << ", " << total_counts_ << ", " << total_sites_ << ", " << func_calls_ << ", " << loglike_ << ", ";

	if(2 == dispersion_.size()){
		output << dispersion_.at(0) << ", " << dispersion_.at(1);
	}
	else{
		output << 0 << ", " << 1e-100;
	}

	for(auto gc = 0; gc < gc_bias_.size(); ++gc){
		output << ", " << unnormalized_gc.at(gc) << ", " << gc_bias_.at(gc) << ", " << gc_count_.at(gc) << ", " << gc_sites_.at(gc);
	}
	for(auto k = 0; k < gc_knots_.size(); ++k){
		output << ", " << static_cast<uintPercentPrint>(gc_knots_.at(k));
	}
	for(auto sur = 0; sur < sur_bias_.size(); ++sur){
		output << ", " << sur_bias_.at(sur) << ", " << sur_count_.at(sur) << ", " << sur_sites_.at(sur);
	}
	output << std::endl;

	ofstream myfile;
	lock_guard<mutex> lock(file_mutex_);
	myfile.open(BiasCalculationVectors::kParameterInfoFile, ofstream::out | ofstream::app);
	myfile << output.str();
	myfile.close();
}

void BiasCalculationVectors::WriteOutDispersionInfo(const char *context){
	// Set predicted site bias
	for(auto &site : sites_){
		site.bias_ = SurBiasAtSite(site) * gc_bias_.at(site.gc_);
	}

	// Sort sites from lowest to highest bias
	sort(sites_.begin(), sites_.end(), FragmentSite::LowerBias);

	// Prepare optimizer
	dispersion_.resize(1);
	nlopt::opt optimizer(nlopt::LD_LBFGS, dispersion_.size());
	optimizer.set_ftol_abs(kStopCriterion);
	optimizer.set_maxeval(kMaxLikelihoodCalculations);

	vector<double> bounds;
	bounds.resize(1);
	bounds.at(0) = 1e-200;
	optimizer.set_lower_bounds(bounds);
	bounds.at(0) = 1e200;
	optimizer.set_upper_bounds(bounds);

	optimizer.set_max_objective(BiasCalculationVectors::LogLikelihoodConstDispersion, reinterpret_cast<void*>(this));

	// Prepare fits
	double loglike, sample_mean, prediction_mean;
	uintDupCount max_duplication;
	vector<double> grad;
	stringstream output;

	// Constant number of sites per bin
	start_site_dispersion_fit_ = ((total_sites_)%kDispersionBinSize)/2; // Ignore equally many sites in the beginning and end
	end_site_dispersion_fit_ = start_site_dispersion_fit_ + kDispersionBinSize/2;

	while( end_site_dispersion_fit_ < sites_.size() ){
		// Get statistics
		sample_mean = 0.0;
		prediction_mean = 0.0;
		max_duplication = 0;
		duplication_count_part_.fill(0);

		for( auto s = start_site_dispersion_fit_; s < end_site_dispersion_fit_; ++s){
			sample_mean += sites_.at(s).count_forward_ + sites_.at(s).count_reverse_;
			prediction_mean += sites_.at(s).bias_;

			SetToMax(max_duplication, sites_.at(s).count_forward_);
			SetToMax(max_duplication, sites_.at(s).count_reverse_);
			++duplication_count_part_.at(min(sites_.at(s).count_forward_, static_cast<uintDupCount>(kMaxDuplications+1)));
			++duplication_count_part_.at(min(sites_.at(s).count_reverse_, static_cast<uintDupCount>(kMaxDuplications+1)));
		}
		sample_mean /= 2*(end_site_dispersion_fit_ - start_site_dispersion_fit_);
		prediction_mean /= (end_site_dispersion_fit_ - start_site_dispersion_fit_);

		output << ref_seq_bin_ << ", " << insert_length_ << ", " << context << ", " << prediction_mean << ", " << sample_mean << ", " << 2*(end_site_dispersion_fit_ - start_site_dispersion_fit_) << ", " << max_duplication;
		for(uintDupCount dups=1; dups<4; ++dups){
			output << ", " << duplication_count_part_.at(dups);
		}

		// Predicted means
		use_sample_mean_ = nan("0");
		dispersion_.at(0) = 1.0;

		func_calls_const_disp_ = 0;
		try{
			if( nlopt::MAXEVAL_REACHED == optimizer.optimize(dispersion_, loglike) ){
				func_calls_const_disp_ += 10000;
			}
		}
		catch(const std::exception& e){
			func_calls_const_disp_ += 20000;
		}
		output << ", " << func_calls_const_disp_ << ", " << dispersion_.at(0);

		// Sample mean
		use_sample_mean_ = sample_mean;
		dispersion_.at(0) = 1.0;

		func_calls_const_disp_ = 0;
		try{
			if( nlopt::MAXEVAL_REACHED == optimizer.optimize(dispersion_, loglike) ){
				func_calls_const_disp_ += 10000;
			}
		}
		catch(const std::exception& e){
			func_calls_const_disp_ += 20000;
		}

		output << ", " << func_calls_const_disp_ << ", " << dispersion_.at(0) << std::endl;

		// Prepare next iteration
		start_site_dispersion_fit_ += kDispersionBinSize/2;
		end_site_dispersion_fit_ += kDispersionBinSize/2;
	}

	// Write output to file
	ofstream myfile;
	lock_guard<mutex> lock(file_mutex_);
	myfile.open(kDispersionInfoFile, ofstream::out | ofstream::app);
	myfile << output.str();
	myfile.close();
}

template<class T> void InsertLengthSpline::GetSampledValues(const Vect<double> &insert_lengths_bias, const T &insert_lengths){
	// Get insert_lengths/insert_lengths_bias_ ratios at sample_positions_
	sampled_values_.resize(sample_positions_.size());
	for(uintSeqLen s=0; s < sample_positions_.size(); ++s){
		// sample_positions_ are selected in a way that count and therefore bias cannot be zero
		sampled_values_.at(s) = insert_lengths.at( sample_positions_.at(s) ) / insert_lengths_bias.at( sample_positions_.at(s) ); // Count/Bias ratio
	}

	// Get mean in case the fit does not work (only relevant for normalization fit)
	sample_mean_ = 0.0;
	for(uintSeqLen s=0; s < sample_positions_.size(); ++s){
		sample_mean_ += sampled_values_.at(s);
	}
}

void InsertLengthSpline::DistributeStartingKnots(){
	// At the beginning every sample gets a knot
	knots_.resize( sample_positions_.size() );
	for( uintSeqLen k=0; k < knots_.size(); ++k ){
		knots_.at(k) = sample_positions_.at(k);
	}
}

void InsertLengthSpline::ConfigureOptimizer(){
	// We need the number of knots before we can configure the optimizer, so run after DistributeKnots
	optimizer_ = nlopt::opt(nlopt::LD_LBFGS, knots_.size()+1);

	optimizer_.set_min_objective(Chi2, reinterpret_cast<void*>(this));

	optimizer_.set_ftol_rel(kStopCriterion);
	optimizer_.set_maxeval(kMaxChi2Calculations);

	// Set the bounds for the knots
	vector<double> bounds;
	const double upper_bound = kUpperBound;
	bounds.resize(knots_.size()+1, upper_bound);
	bounds.at(0) = 1.0; // Normalization parameter that we don't use here
	optimizer_.set_upper_bounds(bounds);

	bounds.clear();
	const double lower_bound = kLowerBound;
	bounds.resize(knots_.size()+1, lower_bound);
	bounds.at(0) = 1.0; // Normalization parameter that we don't use here
	optimizer_.set_lower_bounds(bounds);
}

void InsertLengthSpline::PrepareInsertLengthSpline(){
	// Create coefficient vectors for linear combinations used for the spline
	//b[i][a_index]: i=[0, n-1]
	for(uint8_t dim = 3; dim--; ){
		lin_comb_splines_.at(dim).clear();
		lin_comb_splines_.at(dim).resize(knots_.size()-1);
		for(auto i=lin_comb_splines_.at(dim).size(); i--; ){
			lin_comb_splines_.at(dim).at(i).resize(knots_.size(), 0.0);
		}
	}

	// Create variables needed for PrepareSplines
	vector<double> l; // j=[0, n]
	l.resize(knots_.size(), 0.0);

	vector<vector<double>> c, z; // c[j][a_index]: j=[0, n]
	c.resize(knots_.size());
	z.resize(knots_.size());
	for(uintPercent k=0; k < c.size(); ++k){
		c.at(k).resize(knots_.size(), 0.0);
		z.at(k).resize(knots_.size(), 0.0);
	}

	vector<double> mu, h; // i=[0, n-1]
	mu.resize(knots_.size()-1, 0.0);
	h.resize(knots_.size()-1, 0.0);

	vector<vector<double>> beta; // beta[k][a_index]: k=[1, n-1]
	beta.resize(knots_.size()-1);
	for(uintPercent k=1; k < beta.size(); ++k){
		beta.at(k).resize(knots_.size(), 0.0);
	}

	BiasCalculationVectors::PrepareSplines(lin_comb_splines_, knots_, l, c, z, mu, h, beta);
}

void InsertLengthSpline::SetStartingParameters(){
	pars_.resize(knots_.size()+1);
	pars_.at(0) = 1.0; // That is a normalization parameter that we don't need here
	uintSeqLen s = 0;
	for(uintSeqLen k=0; k < knots_.size(); ++k){
		// Go to the next sampling position after the knot (or on the knot if one exists)
		while( knots_.at(k) > sample_positions_.at(s) ){
			++s;
		}

		if( knots_.at(k) == sample_positions_.at(s) ){
			// Knot is on a sampling position
			pars_.at(k+1) = sampled_values_.at(s);
		}
		else{
			// Knot is in between sampling positions: Mean value
			pars_.at(k+1) = 0.5 * (sampled_values_.at(s-1) + sampled_values_.at(s));
		}

		// We need log because we use exp later on the spline result
		if(pars_.at(k+1) > 0.0){
			pars_.at(k+1) = log( pars_.at(k+1) );
		}
		else{
			pars_.at(k+1) = log( 1e-10 );
		}
	}
}

double InsertLengthSpline::Chi2(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
	InsertLengthSpline &self(*reinterpret_cast<InsertLengthSpline *>(f_data));

	double chi2 = 0.0;
	for(auto &g : grad){
		g = 0.0;
	}

	double a, b, c, d;
	uintSeqLen k = 0;
	BiasCalculationVectors::GetSplineCoefficients(a, b, c, d, k, x, self.lin_comb_splines_);

	for(uintSeqLen s=0; s < self.sample_positions_.size()-1; ++s){
		if( self.sample_positions_.at(s) >= self.knots_.at(k+1) ){
			BiasCalculationVectors::GetSplineCoefficients(a, b, c, d, ++k, x, self.lin_comb_splines_);
		}

		uintSeqLen cur_len = self.sample_positions_.at(s) - self.knots_.at(k);
		double estimate = exp(a + b*cur_len + c*cur_len*cur_len + d*cur_len*cur_len*cur_len);

		chi2 += (self.sampled_values_.at(s) - estimate) * (self.sampled_values_.at(s) - estimate);

		if(grad.size()){
			// Update gradient: Coefficient a
			double factor = -2.0 * (self.sampled_values_.at(s) - estimate) * estimate;
			grad.at(k+1) += factor; // k+1 because parameter 0 is the normalization parameter we ignore

			// Update gradient: Coefficients b, c, d
			factor *= cur_len;
			for(uintPercent ai=1; ai < x.size(); ++ai){
				grad.at(ai) += factor * self.lin_comb_splines_.at(0).at(k).at(ai-1); // b
				grad.at(ai) += factor * cur_len * self.lin_comb_splines_.at(1).at(k).at(ai-1); // c
				grad.at(ai) += factor * cur_len * cur_len * self.lin_comb_splines_.at(2).at(k).at(ai-1); // d
			}
		}
	}

	// Handle last sample, that falls on top of last knot
	double estimate = exp( x.at(x.size()-1) );

	chi2 += (self.sampled_values_.at(self.sample_positions_.size()-1) - estimate) * (self.sampled_values_.at(self.sample_positions_.size()-1) - estimate);

	if(grad.size()){
		grad.at( grad.size()-1 ) = -2.0 * (self.sampled_values_.at(self.sample_positions_.size()-1) - estimate) * estimate;
	}

	return chi2;
}

double InsertLengthSpline::OptimizeParameters(){
	PrepareInsertLengthSpline();
	SetStartingParameters();

	double chi2;
	try{
		optimizer_.optimize(pars_, chi2);
	}
	catch(const std::exception& e){
		// Make sure we have the chi2 that fits the pars_
		std::vector<double> grad;
		chi2 = Chi2(pars_, grad, reinterpret_cast<void*>(this));
	}


	return chi2;
}

double InsertLengthSpline::OptimizeParameters(uintSeqLen removed_knot){
	auto tmp_pos = knots_.at(removed_knot);

	knots_.erase(knots_.begin() + removed_knot);

	ConfigureOptimizer();
	double chi2 = OptimizeParameters();

	knots_.insert(knots_.begin() + removed_knot, tmp_pos);

	return chi2;
}

double InsertLengthSpline::OptimizeParameters(uintSeqLen shifted_knot, intSeqShift shifted_positions){
	if( 0 > shifted_positions ){
		if(-shifted_positions > knots_.at(shifted_knot) || knots_.at(shifted_knot-1) >= knots_.at(shifted_knot) + shifted_positions){
			// Knot shifted over the previous knot (or below zero)
			return numeric_limits<double>::max();
		}
	}
	else{
		if(knots_.at(shifted_knot+1) <= knots_.at(shifted_knot) + shifted_positions){
			// Knot shifted over the next knot
			return numeric_limits<double>::max();
		}
	}

	knots_.at(shifted_knot) += shifted_positions;

	double chi2 = OptimizeParameters();

	knots_.at(shifted_knot) -= shifted_positions;

	return chi2;
}

void InsertLengthSpline::ShiftOptimization(){
	ConfigureOptimizer();
	double min_chi2 = OptimizeParameters();
	double new_chi2;

	// Greedily shift inner knots to minimize chi^2 (optimize parameters for every tested knot combination)
	uintSeqLen best_shifted_knot = 0;
	intSeqShift best_shifted_positions = 0;
	intSeqShift better_direction;

	// Loop until convergence
	do{ // while(0 != best_shifted_positions)
		// Adjust knots based on result from last round
		knots_.at(best_shifted_knot) += best_shifted_positions;
		best_shifted_knot = 0;
		best_shifted_positions = 0;

		// Try to shift every inner knot and see if that improves the result
		for(uintSeqLen shifted_knot = 1; shifted_knot < knots_.size()-1; ++shifted_knot){
			for(intSeqShift shifted_positions = kKnotShiftLength; shifted_positions < kKnotMaxShift; shifted_positions += kKnotShiftLength){
				// Test positive and negative shift
				double new_chi2_positive = OptimizeParameters(shifted_knot, shifted_positions);
				double new_chi2_negative = OptimizeParameters(shifted_knot,-shifted_positions);

				// Select better of the two
				if(new_chi2_negative < new_chi2_positive){
					better_direction = -1;
					new_chi2 = new_chi2_negative;
				}
				else{
					better_direction = 1;
					new_chi2 = new_chi2_positive;
				}

				// Compare to previously tested values
				if(new_chi2 < min_chi2){
					best_shifted_knot = shifted_knot;
					best_shifted_positions = better_direction * shifted_positions;
					min_chi2 = new_chi2;
				}
			}
		}
	} while(0 != best_shifted_positions);
}

void InsertLengthSpline::Optimize(){
	// Get chi2 for starting knot distribution
	double min_chi2;
	uintSeqLen least_important_knot;

	// Find optimal number of knots by removing the least important inner knot until we have the desired number of knots
	uintSeqLen next_shift_optimization = knots_.size()/2;
	while(knots_.size() > kNumKnots){
		// Remove least important knot
		least_important_knot = 0;
		min_chi2 = numeric_limits<double>::max();

		for(uintSeqLen removed_knot = 1; removed_knot < knots_.size()-1; ++removed_knot){
			double cur_chi2 = OptimizeParameters(removed_knot);
			if(cur_chi2 < min_chi2){
				min_chi2 = cur_chi2;
				least_important_knot = removed_knot;
			}
		}

		knots_.erase(knots_.begin() + least_important_knot);

		if( knots_.size() == next_shift_optimization ){
			next_shift_optimization = knots_.size()/2;
			ShiftOptimization();
		}
	}

	// Get parameters for best knot distribution
	ShiftOptimization();
	OptimizeParameters();
}

void InsertLengthSpline::FillInWithFittedRatios(Vect<double> &insert_lengths_bias, const Vect<uintFragCount> &insert_lengths){
	// Everything before the first knot is zero
	for(uintSeqLen len=min(static_cast<size_t>(1),insert_lengths_bias.from()); len < knots_.at(0); ++len){
		insert_lengths_bias.at(len) = 0.0;
	}

	// Interpolate with spline
	uintSeqLen k = 0;
	double a, b, c, d;
	for(; k < knots_.size()-1; ++k){
		BiasCalculationVectors::GetSplineCoefficients(a, b, c, d, k, pars_, lin_comb_splines_);

		insert_lengths_bias.at( knots_.at(k) ) = insert_lengths.at( knots_.at(k) ) / exp(a);

		for(uintSeqLen len = knots_.at(k)+1; len < knots_.at(k+1); ++len){
			uintSeqLen cur_len = len-knots_.at(k);
			insert_lengths_bias.at(len) = insert_lengths.at(len) / exp(a + b*cur_len + c*cur_len*cur_len + d*cur_len*cur_len*cur_len);
		}
	}

	insert_lengths_bias.at( knots_.at(k) ) = insert_lengths.at( knots_.at(k) ) / exp( pars_.at(k+1) ); // Set last knot to its value defined by the parameter

	// Linearly continue values for insert length after last knot
	uintSeqLen cur_len = knots_.at(k) - knots_.at(k-1); // Take slope at the last knot
	double slope = (b + c*cur_len); // d=0 since second derivative must be 0 (natural spline)
	for(uintSeqLen len=knots_.at(k)+1; len < insert_lengths_bias.to(); ++len){
		insert_lengths_bias.at(len) = insert_lengths.at(len) / exp(pars_.at(k+1) + (len-knots_.at(k))*slope);
	}
}

void InsertLengthSpline::FillInWithFittedRatios(vector<double> &normalization, const Vect<double> &insert_lengths_bias){
	// Everything before the first knot is zero
	for(uintSeqLen len=1; len < knots_.at(0); ++len){
		normalization.at(len) = 0.0;
	}

	// Interpolate with spline
	uintSeqLen k = 0;
	double a, b, c, d;
	for(; k < knots_.size()-1; ++k){
		BiasCalculationVectors::GetSplineCoefficients(a, b, c, d, k, pars_, lin_comb_splines_);

		normalization.at( knots_.at(k) ) = insert_lengths_bias.at( knots_.at(k) ) * exp(a);

		for(uintSeqLen len = knots_.at(k)+1; len < knots_.at(k+1); ++len){
			uintSeqLen cur_len = len-knots_.at(k);
			normalization.at(len) = insert_lengths_bias.at(len) * exp(a + b*cur_len + c*cur_len*cur_len + d*cur_len*cur_len*cur_len);
		}
	}

	normalization.at( knots_.at(k) ) = insert_lengths_bias.at( knots_.at(k) ) * exp( pars_.at(k+1) ); // Set last knot to its value defined by the parameter

	// Linearly continue values for insert length after last knot
	uintSeqLen cur_len = knots_.at(k) - knots_.at(k-1); // Take slope at the last knot
	double slope = (b + c*cur_len); // d=0 since second derivative must be 0 (natural spline)
	for(uintSeqLen len=knots_.at(k)+1; len < normalization.size(); ++len){
		normalization.at(len) = insert_lengths_bias.at(len) * exp(pars_.at(k+1) + (len-knots_.at(k))*slope);
	}
}

template<class T> bool InsertLengthSpline::CheckFit(const Vect<double> &insert_lengths_bias, const T &insert_lengths){
	// Check if a estimated ratio is way to high or low and remove the sample next to it, that likely is responsible
	// Start at low insert length and stop after a removal of a sample, fluctuations decrease with insert length and therefore the causal samples are more likely to sit at the beginning
	uintSeqLen s(0);
	for(uintSeqLen len=knots_.at(0); len < knots_.at( knots_.size()-1 ) ; ++len){
		// Update sample
		if(len > sample_positions_.at(s) ){
			++s;
		}

		// Check if ratio values are too high or low
		const double outlier_factor = 1.5;
		uintSeqLen sample_to_delete = 0;

		if(sample_positions_.at(s) == len){
			// We are on a sample so only compare to the correct value
			if( insert_lengths.at(len)/insert_lengths_bias.at(len) < sampled_values_.at(s) / outlier_factor || // Factor too low
				insert_lengths.at(len)/insert_lengths_bias.at(len) > sampled_values_.at(s) * outlier_factor ){ // Factor too high

				if(0 == s){
					// We must not remove first sample, so remove second
					sample_to_delete = 1;
				}
				else{
					sample_to_delete = s;
				}
			}
		}
		else{
			// In between original sample so compare to min/max of the two neighbouring samples
			if( insert_lengths.at(len)/insert_lengths_bias.at(len) < min(sampled_values_.at(s-1), sampled_values_.at(s)) / outlier_factor ){
				// Factor too low
				if(1 == s){
					// We must not remove first sample, so remove second
					sample_to_delete = 1;
				}
				else{
					if(sampled_values_.at(s-1) < sampled_values_.at(s)){
						sample_to_delete = s-1;
					}
					else{
						sample_to_delete = s;
					}
				}
			}
			else if( insert_lengths.at(len)/insert_lengths_bias.at(len) > max(sampled_values_.at(s-1), sampled_values_.at(s)) * outlier_factor ){
				// Factor too high
				if(1 == s){
					// We must not remove first sample, so remove second
					sample_to_delete = 1;
				}
				else{
					if(sampled_values_.at(s-1) > sampled_values_.at(s)){
						sample_to_delete = s-1;
					}
					else{
						sample_to_delete = s;
					}
				}
			}
		}

		if(sample_to_delete){
			// Remove sample
			sample_positions_.erase(sample_positions_.begin()+sample_to_delete);
			sampled_values_.erase(sampled_values_.begin()+sample_to_delete);

			return false;
		}
	}

	return true;
}

void InsertLengthSpline::FillInWithMeanRatio(Vect<double> &insert_lengths_bias, const Vect<uintFragCount> &insert_lengths){
	// We normalize anyways afterwards, so dividing by mean ratio is the same as directly taking insert_length
	for(uintSeqLen len=insert_lengths_bias.from(); len < insert_lengths_bias.to(); ++len){
		insert_lengths_bias.at(len) = insert_lengths.at(len);
	}
}

void InsertLengthSpline::FillInWithMeanRatio(std::vector<double> &normalization, const Vect<double> &insert_lengths_bias){
	for(uintSeqLen len=1; len < normalization.size(); ++len){
		normalization.at(len) = sample_mean_ * insert_lengths_bias.at(len);
	}
}

bool InsertLengthSpline::GetSamplePositions(const Vect<uintFragCount> &insert_lengths){
	// Find first sample that is not zero
	uintSeqLen first_sample = max(static_cast<size_t>(1), insert_lengths.from());
	while(first_sample < insert_lengths.to() && 0 == insert_lengths.at(first_sample)){
		++first_sample;
	}

	// Distribute the samples equidistantly
	uintSeqLen num_samples = 0;
	uintSeqLen hit_zero = 0; // position where we hit zero (0 = nowhere)
	for(uintSeqLen len = first_sample; len < insert_lengths.to(); len += kInsertLengthSamplingDistance){
		if(hit_zero){
			if(insert_lengths.at(len) >= 10){
				// Rescue if we suddenly find higher values again
				num_samples += (len-hit_zero)/kInsertLengthSamplingDistance + 1;
				hit_zero = 0;
			}
		}
		else{
			if(insert_lengths.at(len) > 0){
				++num_samples;
			}
			else{
				// End knots when we find the first zero
				hit_zero = len;
			}
		}
	}

	if(2 > num_samples){
		printErr << "Sampling insert lengths did not find at least two usable lengths." << std::endl;
		return false;
	}

	sample_positions_.resize(num_samples);
	sample_positions_.at(0) = first_sample;
	uintSeqLen found_zeros = 0;
	for(uintSeqLen k=1; k<num_samples-found_zeros; ++k){
		sample_positions_.at(k) = sample_positions_.at(k-1) + kInsertLengthSamplingDistance;

		// Exclude samples at fragment lengths with no data
		while(0 == insert_lengths.at( sample_positions_.at(k) )){ // Last sample cannot be a zero, so no further checks necessary
			++found_zeros;
			sample_positions_.at(k) += kInsertLengthSamplingDistance;
		}
	}
	sample_positions_.resize(num_samples-found_zeros);

	return true;
}

bool InsertLengthSpline::FitInsertLengthWithSpline(Vect<double> &insert_lengths_bias, const Vect<uintFragCount> &insert_lengths){
	GetSampledValues(insert_lengths_bias, insert_lengths);

	bool fit_accepted;
	do{
		DistributeStartingKnots();
		Optimize();

		FillInWithFittedRatios(insert_lengths_bias, insert_lengths);

		fit_accepted = CheckFit(insert_lengths_bias, insert_lengths);
		if(sample_positions_.size() < 2){
			printWarn << "Count/Bias ratio could not be fitted for fragment length. Using mean ratio.";
			FillInWithMeanRatio(insert_lengths_bias, insert_lengths);
			return true;
		}
	} while(!fit_accepted);


	return true;
}

bool InsertLengthSpline::InterpolateNormalizationWithSpline(vector<double> &normalization, const Vect<double> &insert_lengths_bias){
	GetSampledValues(insert_lengths_bias, normalization);

	DistributeStartingKnots();

	PrepareInsertLengthSpline();
	SetStartingParameters();

	FillInWithFittedRatios(normalization, insert_lengths_bias);

	return true;
}

inline void FragmentDistributionStats::IncreaseErrorCounter(uintErrorCount &errors){
	if( errors++ == kMaxErrorsShownPerFile ){
		printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
	}
}

void FragmentDistributionStats::PrepareBiasCalculation( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const vector<uintFragCount> &reads_per_frag_len_bin, const vector<uintFragCount> &lowq_reads_per_frag_len_bin ){
	// Bin reference sequence
	auto num_bins = CreateRefBins(ref, max_ref_seq_bin_size);

	vector<uintFragCount> reads_per_ref_seq_bin;
	reads_per_ref_seq_bin.resize(num_bins);
	vector<uintFragCount> lowq_reads_per_ref_seq_bin;
	lowq_reads_per_ref_seq_bin.resize(num_bins);

	uintRefLenCalc ref_seq_start_bin(0);
	uintRefSeqBin ref_seq_bin(0);
	for(uintRefSeqId ref_seq = 0; ref_seq < ref.NumberSequences(); ++ref_seq){
		if( !ref.ReferenceSequenceExcluded(ref_seq) ){
			for(uintSeqLen cur_bin = 0; cur_bin <= ref.SequenceLength(ref_seq)/maximum_insert_length; ++cur_bin){
				ref_seq_bin = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length, ref_seq_bin);
				reads_per_ref_seq_bin.at(ref_seq_bin) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
				lowq_reads_per_ref_seq_bin.at(ref_seq_bin) += lowq_reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);

				// Add it also to all bins the fragment might overlap, as we don' know the exact start position of the fragment
				if(cur_bin < ref.SequenceLength(ref_seq)/maximum_insert_length){
					auto ref_seq_bin2 = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length + maximum_insert_length, ref_seq_bin);
					if(ref_seq_bin != ref_seq_bin2){
						reads_per_ref_seq_bin.at(ref_seq_bin2) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
						lowq_reads_per_ref_seq_bin.at(ref_seq_bin2) += lowq_reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
					}
				}

				if(cur_bin){
					auto ref_seq_bin2 = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length - maximum_insert_length, ref_seq_bin);
					if(ref_seq_bin != ref_seq_bin2){
						reads_per_ref_seq_bin.at(ref_seq_bin2) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
						lowq_reads_per_ref_seq_bin.at(ref_seq_bin2) += lowq_reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
					}
				}
			}
		}
		ref_seq_start_bin += ref.SequenceLength(ref_seq)/maximum_insert_length + 1;
	}

	// Read in vectors to calculate bias from
	fragment_sites_by_ref_seq_bin_.resize(num_bins);
	for( uintRefSeqBin ref_bin=0; ref_bin < fragment_sites_by_ref_seq_bin_.size(); ++ref_bin){
		fragment_sites_by_ref_seq_bin_.at(ref_bin).resize(reads_per_ref_seq_bin.at(ref_bin)/2); // Sites are defined by a read pair, so /2
	}
	fragment_sites_by_ref_seq_bin_cur_id_.resize(num_bins);
	fragment_sites_by_ref_seq_bin_by_insert_length_.resize(num_bins);

	lowq_site_start_by_ref_seq_bin_.resize(num_bins);
	lowq_site_end_by_ref_seq_bin_.resize(num_bins);
	for( uintRefSeqBin ref_bin=0; ref_bin < lowq_site_start_by_ref_seq_bin_.size(); ++ref_bin){
		lowq_site_start_by_ref_seq_bin_.at(ref_bin).resize(lowq_reads_per_ref_seq_bin.at(ref_bin));
		lowq_site_end_by_ref_seq_bin_.at(ref_bin).resize(lowq_reads_per_ref_seq_bin.at(ref_bin));
	}
	lowq_site_start_by_ref_seq_bin_cur_id_.resize(num_bins);
	lowq_site_end_by_ref_seq_bin_cur_id_.resize(num_bins);
	lowq_site_start_exclusion_.resize(num_bins);
	lowq_site_end_exclusion_.resize(num_bins);
	excluded_lowq_positions_ = 0;
	excluded_lowq_regions_ = 0;

	corrected_abundance_.resize(ref.NumberSequences());
	filtered_sites_.resize(ref.NumberSequences());
	fragment_lengths_used_for_correction_.resize(ref.NumberSequences());

	// Prepare bias fitting
	bias_calc_params_.resize( maximum_insert_length * kMaxBinsQueuedForBiasCalc );

	ref.RefSeqsInNxx(ref_seq_in_nxx_, BiasCalculationVectors::kNXXRefSeqs);

	// Temporary bias fit result vectors
	params_fitted_ = 0;
	auto max_fits = ref.NumRefSeqBinsInNxx(ref_seq_in_nxx_, max_ref_seq_bin_length_);
	if(max_fits >= 20){
		num_fits_per_ref_bin_= 5;
	}
	else if(max_fits >= 17){
		num_fits_per_ref_bin_= 6;
	}
	else if(max_fits >= 15){
		num_fits_per_ref_bin_= 7;
	}
	else if(max_fits >= 13){
		num_fits_per_ref_bin_= 8;
	}
	else if(max_fits >= 10){
		num_fits_per_ref_bin_= 9;
	}
	else if(max_fits >= 6){
		num_fits_per_ref_bin_= 10;
	}
	else if(max_fits >= 5){
		num_fits_per_ref_bin_= 12;
	}
	else if(max_fits >= 3){
		num_fits_per_ref_bin_= 15;
	}
	else if(max_fits >= 2){
		num_fits_per_ref_bin_= 20;
	}
	else{
		num_fits_per_ref_bin_= 30;
	}

	max_fits *= num_fits_per_ref_bin_;
	for(uintPercent gc=tmp_gc_bias_.size(); gc--;){
		tmp_gc_bias_.at(gc).resize(max_fits, {nan("0"),0});
	}
	for(uintNumFits sur=tmp_sur_bias_.size(); sur--;){
		tmp_sur_bias_.at(sur).resize(max_fits, nan("0"));
	}
	for(uintNumFits par=tmp_dispersion_parameters_.size(); par--;){
		tmp_dispersion_parameters_.at(par).resize(max_fits, nan("0"));
	}

	// Initialize bias calculation queue information
	params_left_for_calculation_ = ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1) * maximum_insert_length;
	for( uintNumFits bin = 0; bin < current_bias_param_.size(); ++bin){
		current_bias_param_.at(bin) = maximum_insert_length+1; // So no thread tries to calculate on the unfilled bins
		finished_bias_calcs_.at(bin) = 0;
		claimed_bias_bins_.at(bin).clear();
	}

	// Creating files for additional output
	if( BiasCalculationVectors::kParameterInfoFile ){
		ofstream myfile;
		myfile.open(BiasCalculationVectors::kParameterInfoFile);
		myfile << "Fit, RefSeqBin, InsertLength, Counts, Sites, FunctionCalls, LogLikelihood, DispersionA, DispersionB";
		for(auto gc = 0; gc < 101; ++gc){
			myfile << ", " << "RawGCbias" << gc << ", " << "GCbias" << gc << ", " << "GCcount" << gc << ", " << "GCsites" << gc;
		}
		for(auto k = 0; k < BiasCalculationVectors::kGCSplineDf; ++k){
			myfile << ", " << "GCknot" << k;
		}
		for(auto sur = 0; sur < 4*Surrounding::Length(); ++sur){
			char base;
			switch(sur%4){
			case 0:
				base='A';
				break;
			case 1:
				base='C';
				break;
			case 2:
				base='G';
				break;
			case 3:
				base='T';
				break;
			}
			myfile << ", " << "SurBias" << sur/4 << base << ", " << "SurCount" << sur/4 << base << ", " << "SurSites" << sur/4 << base;
		}
		myfile << std::endl;
		myfile.close();
	}

	if( BiasCalculationVectors::kDispersionInfoFile ){
		ofstream myfile;
		myfile.open(BiasCalculationVectors::kDispersionInfoFile);
		myfile << "ref_seq_bin, insert_length, mean_fit, prediction_mean, sample_mean, num_sites, max_duplication, counts1, counts2, counts3, calls_prediction, dispersion_prediction, calls_sample, dispersion_sample" << std::endl;
		myfile.close();
	}
}

inline void FragmentDistributionStats::CheckLowQExclusion(vector<pair<uintSeqLen,uintSeqLen>> &lowq_site_exclusion, uintSeqLen excluded_pos, uintFragCount lowq_count, const vector<uintFragCount> &tmp_frag_count ){
	if(lowq_count > kMinLowQFractionForExclusion*tmp_frag_count.at(excluded_pos) - 1e-10){
		// Low quality is not just spurious, so exclude this position
		if(0 == lowq_site_exclusion.size()){
			lowq_site_exclusion.emplace_back( excluded_pos, excluded_pos+1 );
		}
		else{
			uintFragCount separating_high_q(0);
			for( auto last_pos = lowq_site_exclusion.back().second; last_pos < excluded_pos && separating_high_q < kMinHighQToSeparateLowQ; ++last_pos ){
				separating_high_q += tmp_frag_count.at(last_pos);
			}

			if(separating_high_q < kMinHighQToSeparateLowQ){
				// Connect the two regions
				if(1 == lowq_site_exclusion.back().second - lowq_site_exclusion.back().first ){
					// Extend the the region in the beginning
					separating_high_q = 0;
					while(lowq_site_exclusion.back().first && separating_high_q + tmp_frag_count.at(lowq_site_exclusion.back().first-1) < kMinHighQToSeparateLowQ/2 ){
						separating_high_q += tmp_frag_count.at(--lowq_site_exclusion.back().first);
					}
				}
				lowq_site_exclusion.back().second = excluded_pos+1;
			}
			else{
				if(1 == lowq_site_exclusion.back().second - lowq_site_exclusion.back().first ){
					// Replace the single random event
					lowq_site_exclusion.back().first = excluded_pos;
					lowq_site_exclusion.back().second = excluded_pos+1;
				}
				else{
					// Extend last region at the end
					separating_high_q = 0;
					while(lowq_site_exclusion.back().second < tmp_frag_count.size() && separating_high_q + tmp_frag_count.at(lowq_site_exclusion.back().second) < kMinHighQToSeparateLowQ/2 ){
						separating_high_q += tmp_frag_count.at(lowq_site_exclusion.back().second++);
					}

					// Insert new position
					lowq_site_exclusion.emplace_back( excluded_pos, excluded_pos+1 );
				}
			}
		}
	}
}

void FragmentDistributionStats::CheckLowQExclusions(uintRefSeqBin ref_seq_bin, vector<uintFragCount> &tmp_frag_count, const Reference &reference){
	// Handle low quality site starts
	auto ref_seq = ref_seq_bin_def_.at(ref_seq_bin).first;
	auto tmp_size = RefSeqSplitLength(ref_seq, reference) + ref_seq_start_bin_.at(ref_seq+1)-ref_seq_start_bin_.at(ref_seq);
	tmp_frag_count.clear();
	tmp_frag_count.resize(tmp_size, 0);

	for( uintFragCount site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
		++tmp_frag_count.at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).first/2);
	}

	lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).resize( lowq_site_start_by_ref_seq_bin_cur_id_.at(ref_seq_bin) );
	if( lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).size() ){
		lowq_site_start_exclusion_.at(ref_seq_bin).reserve( lowq_site_start_by_ref_seq_bin_cur_id_.at(ref_seq_bin) );
		sort(lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).begin(), lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).end());
		auto excluded_pos = lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).at(0);
		uintFragCount lowq_count(1);

		for( uintFragCount site_id=1; site_id < lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).size(); ++site_id ){
			if(lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).at(site_id) == excluded_pos){
				++lowq_count;
			}
			else{
				CheckLowQExclusion(lowq_site_start_exclusion_.at(ref_seq_bin), excluded_pos, lowq_count, tmp_frag_count );

				excluded_pos = lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).at(site_id);
				lowq_count = 1;
			}
		}

		CheckLowQExclusion(lowq_site_start_exclusion_.at(ref_seq_bin), excluded_pos, lowq_count, tmp_frag_count );

		if(1 == lowq_site_start_exclusion_.at(ref_seq_bin).back().second - lowq_site_start_exclusion_.at(ref_seq_bin).back().first ){
			// Remove the single random event
			lowq_site_start_exclusion_.at(ref_seq_bin).pop_back();
		}
		else{
			// Extend last region at the end
			uintFragCount separating_high_q = 0;
			while(lowq_site_start_exclusion_.at(ref_seq_bin).back().second < tmp_frag_count.size() && separating_high_q + tmp_frag_count.at(lowq_site_start_exclusion_.at(ref_seq_bin).back().second) < kMinHighQToSeparateLowQ/2 ){
				separating_high_q += tmp_frag_count.at(lowq_site_start_exclusion_.at(ref_seq_bin).back().second++);
			}
		}

		lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).clear();
		lowq_site_start_by_ref_seq_bin_.at(ref_seq_bin).shrink_to_fit();
		lowq_site_start_exclusion_.at(ref_seq_bin).shrink_to_fit();
	}

	// Handle low quality site ends
	tmp_frag_count.clear();
	tmp_frag_count.resize(tmp_size+tmp_insert_lengths_.size()+1, 0);
	for( uintFragCount site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
		++tmp_frag_count.at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).first/2 + fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).second);
	}

	lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).resize( lowq_site_end_by_ref_seq_bin_cur_id_.at(ref_seq_bin) );
	if( lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).size() ){
		lowq_site_end_exclusion_.at(ref_seq_bin).reserve( lowq_site_end_by_ref_seq_bin_cur_id_.at(ref_seq_bin) );
		sort(lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).begin(), lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).end());
		auto excluded_pos = lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).at(0);
		uintFragCount lowq_count(1);

		for( uintFragCount site_id=1; site_id < lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).size(); ++site_id ){
			if(lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).at(site_id) == excluded_pos){
				++lowq_count;
			}
			else{
				CheckLowQExclusion(lowq_site_end_exclusion_.at(ref_seq_bin), excluded_pos, lowq_count, tmp_frag_count );

				excluded_pos = lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).at(site_id);
				lowq_count = 1;
			}
		}

		CheckLowQExclusion(lowq_site_end_exclusion_.at(ref_seq_bin), excluded_pos, lowq_count, tmp_frag_count );

		if(1 == lowq_site_end_exclusion_.at(ref_seq_bin).back().second - lowq_site_end_exclusion_.at(ref_seq_bin).back().first ){
			// Remove the single random event
			lowq_site_end_exclusion_.at(ref_seq_bin).pop_back();
		}
		else{
			// Extend last region at the end
			uintFragCount separating_high_q = 0;
			while(lowq_site_end_exclusion_.at(ref_seq_bin).back().second < tmp_frag_count.size() && separating_high_q + tmp_frag_count.at(lowq_site_end_exclusion_.at(ref_seq_bin).back().second) < kMinHighQToSeparateLowQ/2 ){
				separating_high_q += tmp_frag_count.at(lowq_site_end_exclusion_.at(ref_seq_bin).back().second++);
			}
		}

		lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).clear();
		lowq_site_end_by_ref_seq_bin_.at(ref_seq_bin).shrink_to_fit();
		lowq_site_end_exclusion_.at(ref_seq_bin).shrink_to_fit();
	}

	// Use start as we add positions close to bin ends twice to end
	excluded_lowq_positions_ += ExcludedLowQBases( lowq_site_start_exclusion_.at(ref_seq_bin) );
	excluded_lowq_regions_ += lowq_site_start_exclusion_.at(ref_seq_bin).size();
}

void FragmentDistributionStats::SortFragmentSites(uintRefSeqBin ref_seq_bin, vector<uintSeqLen> &num_sites_per_insert_length){
	// Sort sites into insert_length bins
	// Start by getting size of each bin
	num_sites_per_insert_length.clear();
	num_sites_per_insert_length.resize(tmp_insert_lengths_.size(), 0);
	for( uintFragCount site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
		++num_sites_per_insert_length.at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).second);
	}

	uintSeqLen max_used_length = tmp_insert_lengths_.size();
	while(--max_used_length && 0 == num_sites_per_insert_length.at(max_used_length) ); // Ignore empty bins at the end
	if(max_used_length){
		// Reserve needed space
		fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).resize(max_used_length+1);
		for(uintSeqLen insert_length = 1; insert_length < max_used_length+1; ++insert_length ){
			fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).reserve(num_sites_per_insert_length.at(insert_length));
		}

		// Sort into bins
		for( uintFragCount site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
			fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).second).push_back(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).first);
		}
	}

	// Free space that is not needed anymore
	fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).clear();
	fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).shrink_to_fit();
}

void FragmentDistributionStats::UpdateBiasCalculationParams(uintRefSeqBin ref_seq_bin, uint32_t queue_spot, vector<pair<uintSeqLen, pair<uintRefSeqBin, uintSeqLen>>> &tmp_params, mutex &print_mutex){
	auto qbin_size=bias_calc_params_.size()/kMaxBinsQueuedForBiasCalc;
	if(fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size()){
		// Select ref_seq, insert_length bins to use for bias calculation
		if(ref_seq_in_nxx_.at( GetRefSeqId(ref_seq_bin) )){
			// Sort them by counts and only take num_fits_per_ref_bin_ from top kNumFitsInsertLength
			tmp_params.clear();
			for(uintSeqLen insert_length=1; insert_length < fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); ++insert_length){
				if( fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size() ){
					tmp_params.push_back({fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size(), {ref_seq_bin, insert_length}});
				}
				else{
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Clear(ref_seq_bin); // If we won't set anything to it later we have to clear it from it's previous content
				}
			}

			sort(tmp_params.begin(), tmp_params.end(), std::greater<pair<uintSeqLen, pair<uintRefSeqBin, uintSeqLen>>>());

			for( uintNumFits n=0; n < tmp_params.size(); ++n ){
				bool calculate_param = false;
				if(n < BiasCalculationVectors::kNumFitsInsertLength && (0 == n || 1 == n*num_fits_per_ref_bin_/BiasCalculationVectors::kNumFitsInsertLength - (n-1)*num_fits_per_ref_bin_/BiasCalculationVectors::kNumFitsInsertLength)){
					calculate_param = calculate_bias_;
				}

				bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, calculate_param);
			}

			if(tmp_params.size() < num_fits_per_ref_bin_){
				params_fitted_ += num_fits_per_ref_bin_-tmp_params.size(); // Do not have to be fitted, so are already done

				lock_guard<mutex> lock(print_mutex);
				printInfo << "Finished " << static_cast<uintPercentPrint>(Percent(params_fitted_, static_cast<uintNumFits>(tmp_gc_bias_.at(0).size()))) << "% of the bias fits." << std::endl;
			}
		}
		else{
			for(uintSeqLen insert_length=1; insert_length < fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); ++insert_length){
				if( fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size() ){
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Set(ref_seq_bin, insert_length, false);
				}
				else{
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Clear(ref_seq_bin); // If we won't set anything to it later we have to clear it from it's previous content
				}
			}
		}

		// Clear all the other insert_length, so no old values from other ref bins are in there
		for(uintSeqLen insert_length=fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); insert_length <= qbin_size; ++insert_length){
			bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Clear(ref_seq_bin);
		}
		// Reset queue bin
		finished_bias_calcs_.at(queue_spot) = 0;
		current_bias_param_.at(queue_spot) = 0; // This allows other threads to do bias calculations with the values from bias_calc_params_ for this queue bin, so do this as the very last step
	}
	else{
		params_left_for_calculation_ -= qbin_size;
		claimed_bias_bins_.at(queue_spot).clear(); // Let go of claim as reference sequence bin does not need to be handled as it is empty
	}
}

void FragmentDistributionStats::FillParams(vector<BiasCalculationParams> &params, const Reference &reference, const vector<uintSeqLen> &sample_frag_length) const{
	params.reserve( sample_frag_length.size() * abundance_.size() );

	for( uintRefSeqId ref_id=abundance_.size(); ref_id--; ){
		if(!reference.ReferenceSequenceExcluded(ref_id)){
			for( auto frag_length : sample_frag_length ){
				// Add parameters for later calculation
				params.push_back({ref_id, frag_length});
			}
		}
	}

	params.shrink_to_fit();
}

void FragmentDistributionStats::FillParamsSimulation(vector<BiasCalculationParams> &params, const Reference &reference, const vector<uintSeqLen> &sample_frag_length) const{
	params.reserve( insert_lengths_.size() * abundance_.size() );

	for( uintRefSeqId ref_id=ref_seq_bias_.size(); ref_id--; ){
		if(0.0 != ref_seq_bias_.at(ref_id)){
			for( auto frag_length : sample_frag_length ){
				if(frag_length <= reference.SequenceLength(ref_id)){
					// Add parameters for later calculation
					params.push_back({ref_id, frag_length});
				}
			}
		}
	}

	params.shrink_to_fit();
}

void FragmentDistributionStats::CountEndExclusion(uintFragCount &lowq_sites, uintSeqLen &cur_end_exclusion, uintSeqLen last_cut, uintSeqLen new_cut, uintRefSeqBin ref_seq_bin ){
	while( cur_end_exclusion < lowq_site_end_exclusion_.at(ref_seq_bin).size() && last_cut > lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).second ){
		++cur_end_exclusion; // Ignore everything where we cannot end, because it is before fragment length, or where we are overlapping with start exclusion
	}
	if( cur_end_exclusion < lowq_site_end_exclusion_.at(ref_seq_bin).size() && last_cut > lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).first ){
		if( lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).second <= new_cut ){
			lowq_sites += lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).second - last_cut;
			++cur_end_exclusion;
		}
		else{
			lowq_sites += new_cut - last_cut;
		}
	}

	while( cur_end_exclusion < lowq_site_end_exclusion_.at(ref_seq_bin).size() && lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).second <= new_cut ){
		lowq_sites += lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).second - lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).first;
		++cur_end_exclusion;
	}

	if( cur_end_exclusion < lowq_site_end_exclusion_.at(ref_seq_bin).size() && lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).first < new_cut && last_cut <= lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).first ){
		lowq_sites += new_cut - lowq_site_end_exclusion_.at(ref_seq_bin).at(cur_end_exclusion).first;
	}
}

void FragmentDistributionStats::CountDuplicates(FragmentDuplicationStats &duplications, const BiasCalculationParamsSplitSeqs &params, const Reference &reference){
	duplications.AddDuplicates( fragment_sites_by_ref_seq_bin_by_insert_length_.at(params.ref_seq_bin_).at(params.fragment_length_) );

	// fragment_sites_by_ref_seq_bin_by_insert_length_ already sorted from AddDuplicates call
	uintFragCount highq_counts( fragment_sites_by_ref_seq_bin_by_insert_length_.at(params.ref_seq_bin_).at(params.fragment_length_).size() );
	uintSeqLen cur_start_exclusion(0);
	uintSeqLen cur_end_exclusion(0);
	for( auto pos : fragment_sites_by_ref_seq_bin_by_insert_length_.at(params.ref_seq_bin_).at(params.fragment_length_) ){
		while(cur_start_exclusion < lowq_site_start_exclusion_.at(params.ref_seq_bin_).size() && pos >= lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).second ){
			++cur_start_exclusion;
		}
		while(cur_end_exclusion < lowq_site_end_exclusion_.at(params.ref_seq_bin_).size() && pos+params.fragment_length_ >= lowq_site_end_exclusion_.at(params.ref_seq_bin_).at(cur_end_exclusion).second ){
			++cur_end_exclusion;
		}

		if( ( cur_start_exclusion < lowq_site_start_exclusion_.at(params.ref_seq_bin_).size() && pos >= lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first ) ||
			( cur_end_exclusion < lowq_site_end_exclusion_.at(params.ref_seq_bin_).size() && pos+params.fragment_length_ >= lowq_site_end_exclusion_.at(params.ref_seq_bin_).at(cur_end_exclusion).first ) ){
			--highq_counts;
		}
	}

	uintSeqLen bin_end;
	if(params.ref_seq_bin_+1 < ref_seq_bin_def_.size() && ref_seq_bin_def_.at(params.ref_seq_bin_).first == ref_seq_bin_def_.at(params.ref_seq_bin_+1).first){
		bin_end = ref_seq_bin_def_.at(params.ref_seq_bin_+1).second;
	}
	else{
		bin_end = reference.EndExclusion(ref_seq_bin_def_.at(params.ref_seq_bin_).first) - params.fragment_length_+1;
	}
	uintSeqLen correction, current_region;
	if(0 == ref_seq_bin_def_.at(params.ref_seq_bin_).second){
		correction = reference.StartExclusion( ref_seq_bin_def_.at(params.ref_seq_bin_).first );
		current_region = 0;
	}
	else{
		correction = 0;
		current_region = reference.NextExclusionRegion(ref_seq_bin_def_.at(params.ref_seq_bin_).first, ref_seq_bin_def_.at(params.ref_seq_bin_).second) - 1; // The -1 is important as this function returns the region after the position, but CorrectPositionRemovingExcludedRegions needs the one before the position
	}

	reference.CorrectPositionRemovingExcludedRegions(correction, current_region, ref_seq_bin_def_.at(params.ref_seq_bin_).first, bin_end);

	uintFragCount highq_sites = bin_end - ref_seq_bin_def_.at(params.ref_seq_bin_).second;
	highq_sites -= correction;

	uintFragCount lowq_sites(0);
	cur_start_exclusion = 0;
	cur_end_exclusion = 0;
	uintSeqLen last_cut = params.fragment_length_;
	while(cur_start_exclusion < lowq_site_start_exclusion_.at(params.ref_seq_bin_).size() && highq_sites >= lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).second){
		CountEndExclusion(lowq_sites, cur_end_exclusion, last_cut, lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first + params.fragment_length_, params.ref_seq_bin_ );

		lowq_sites += lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).second - lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first;
		last_cut = lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).second + params.fragment_length_;
		++cur_start_exclusion;
	}
	if(cur_start_exclusion < lowq_site_start_exclusion_.at(params.ref_seq_bin_).size() && highq_sites > lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first){
		CountEndExclusion(lowq_sites, cur_end_exclusion, last_cut, lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first + params.fragment_length_, params.ref_seq_bin_ );

		lowq_sites += highq_sites - lowq_site_start_exclusion_.at(params.ref_seq_bin_).at(cur_start_exclusion).first;
	}
	else{
		CountEndExclusion(lowq_sites, cur_end_exclusion, last_cut, highq_sites + params.fragment_length_, params.ref_seq_bin_ );
	}
	highq_sites -= lowq_sites;

	if( highq_sites ){
		auto end_pos = reference.SequenceLength(ref_seq_bin_def_.at(params.ref_seq_bin_).first)-params.fragment_length_+1; // bin_end above removes the EndExclusion
		if(params.ref_seq_bin_+1 < ref_seq_bin_def_.size() && ref_seq_bin_def_.at(params.ref_seq_bin_).first == ref_seq_bin_def_.at(params.ref_seq_bin_+1).first){
			end_pos = ref_seq_bin_def_.at(params.ref_seq_bin_+1).second;
		}

		corrected_abundance_.at(ref_seq_bin_def_.at(params.ref_seq_bin_).first) += Divide( highq_counts * (end_pos - ref_seq_bin_def_.at(params.ref_seq_bin_).second), highq_sites );
		filtered_sites_.at(ref_seq_bin_def_.at(params.ref_seq_bin_).first) += highq_sites;
		++fragment_lengths_used_for_correction_.at(ref_seq_bin_def_.at(params.ref_seq_bin_).first);
	}
}

void FragmentDistributionStats::AddFragmentsToSites(vector<FragmentSite> &sites, const vector<uintSeqLen> &fragment_positions){
	for(auto pos : fragment_positions){
		if(pos%2){
			++sites.at(pos/2).count_reverse_;
		}
		else{
			++sites.at(pos/2).count_forward_;
		}
	}
}

reseq::uintSeqLen FragmentDistributionStats::ExcludeLowQualitySites(vector<FragmentSite> &sites, uintRefSeqBin ref_seq_bin, uintSeqLen insert_length){
	for(auto &region : lowq_site_start_exclusion_.at(ref_seq_bin)){
		for( auto pos = region.first; pos < region.second && pos < sites.size(); ++pos ){ // Depending on insert_length the last bases of reference sequence cannot be a start position
			sites.at(pos).bias_ = 0.0;
		}
	}
	for(auto &region : lowq_site_end_exclusion_.at(ref_seq_bin)){
		for( auto pos = max(insert_length, region.first); pos < region.second && pos-insert_length < sites.size(); ++pos ){
			sites.at(pos-insert_length).bias_ = 0.0;
		}
	}

	uintSeqLen excluded_sites(0);
	for(auto &site : sites){
		if(0.0 == site.bias_){
			++excluded_sites;
		}
	}

	return excluded_sites;
}

void FragmentDistributionStats::CalculateBiasByBin(BiasCalculationVectors &tmp_calc, const Reference &reference, FragmentDuplicationStats &duplications, uintRefSeqBin ref_seq_bin, uintSeqLen insert_length){
	// Prepare fits
	auto ref_seq_id = GetRefSeqId(ref_seq_bin);
	uintSeqLen ref_start = ref_seq_bin_def_[ref_seq_bin].second;
	uintSeqLen ref_end = reference.SequenceLength(ref_seq_id);
	if(ref_seq_bin+1 < ref_seq_start_bin_.at(ref_seq_id+1)){ // Not the last bin of the reference
		ref_end = ref_seq_bin_def_[ref_seq_bin+1].second;
	}

	tmp_calc.ref_seq_bin_ = ref_seq_bin;
	tmp_calc.insert_length_ = insert_length;
	tmp_calc.converged_ = false;

	reference.GetFragmentSites(tmp_calc.sites_, ref_seq_id, insert_length, ref_start, ref_end);
	if(0 == tmp_calc.sites_.size()){
		return;
	}
	AddFragmentsToSites(tmp_calc.sites_, fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length) );
	duplications.AddSites(tmp_calc.sites_);
	auto excluded_sites = ExcludeLowQualitySites(tmp_calc.sites_, ref_seq_bin, insert_length);
	tmp_calc.GetCounts();
	if(0 == tmp_calc.total_counts_){
		return;
	}

	// CountsInHighQRegion * (SitesInHighQRegion + SitesInLowQRegion)/SitesInHighQRegion ~ CountsInHighQRegion + CountsInLowQRegion
	auto end_pos = reference.SequenceLength(ref_seq_bin_def_.at(ref_seq_bin).first)-insert_length+1;
	if(ref_seq_bin+1 < ref_seq_bin_def_.size() && ref_seq_bin_def_.at(ref_seq_bin).first == ref_seq_bin_def_.at(ref_seq_bin+1).first){
		end_pos = ref_seq_bin_def_.at(ref_seq_bin+1).second;
	}

	corrected_abundance_.at(ref_seq_bin_def_.at(ref_seq_bin).first) += Divide( tmp_calc.total_counts_ * (end_pos-ref_seq_bin_def_.at(ref_seq_bin).second), tmp_calc.sites_.size()-excluded_sites );
	filtered_sites_.at(ref_seq_bin_def_.at(ref_seq_bin).first) += tmp_calc.sites_.size()-excluded_sites;
	++fragment_lengths_used_for_correction_.at(ref_seq_bin_def_.at(ref_seq_bin).first);

	tmp_calc.RemoveUnnecessarySites();
	if(0 == tmp_calc.total_sites_){
		return;
	}


	if(kMinSites > tmp_calc.total_sites_ || kMinCounts > tmp_calc.total_counts_ || tmp_calc.total_sites_ > kMaxSitesPerCount*tmp_calc.total_counts_){
		return;
	}

	tmp_calc.CalculateGCWeights();
	tmp_calc.GetLogLikeBase();

	// Fit poisson with independent gc bins
	tmp_calc.fit_pars_.clear();
	if(BiasCalculationVectors::kSurMult){
		tmp_calc.fit_pars_.resize(tmp_calc.sur_bias_.size(), 1.0);
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else{
		tmp_calc.fit_pars_.resize(tmp_calc.sur_bias_.size()+1, 0.0);
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::kUpperBound, 1);
		tmp_calc.fit_pars_.at(0) = 1.0;
	}

	tmp_calc.func_calls_ = 0;

	tmp_calc.optimizer_poisson_.set_max_objective(BiasCalculationVectors::LogLikelihoodPoisson, reinterpret_cast<void*>(&tmp_calc));

	double lower_bound;
	if(BiasCalculationVectors::kSurMult){
		lower_bound = BiasCalculationVectors::kBaseValue;
	}
	else{
		lower_bound = BiasCalculationVectors::kLowerBound;
	}
	const double upper_bound = BiasCalculationVectors::kUpperBound;
	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), lower_bound);
	if(BiasCalculationVectors::kSurMult){
		tmp_calc.DeactivateZeroCounts(tmp_calc.bounds_, 0.0, 0);
	}
	tmp_calc.optimizer_poisson_.set_lower_bounds(tmp_calc.bounds_);

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), upper_bound);
	if(BiasCalculationVectors::kSurMult){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else{
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::kUpperBound, 1);
	}
	tmp_calc.optimizer_poisson_.set_upper_bounds(tmp_calc.bounds_);

	tmp_calc.converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == tmp_calc.optimizer_poisson_.optimize(tmp_calc.fit_pars_, tmp_calc.loglike_) ){
			//printWarn << "Optimizer reached maximum number of likelihood calculations for seq " << ref_seq_bin << " and length " << insert_length << std::endl;
			tmp_calc.func_calls_ += 10000;
			tmp_calc.converged_ = false;
		}

	}
	catch(const std::exception& e){
		//printWarn << "Poisson optimizer crashed for seq " << ref_seq_bin << " and length " << insert_length << " with: " << e.what() << std::endl;
		tmp_calc.func_calls_ += 20000;
		tmp_calc.converged_ = false;
	}

	// Fit gc spline with fixed surrounding bias assuming poisson
	tmp_calc.optimizer_spline_.set_max_objective(BiasCalculationVectors::LogLikeGcSpline, reinterpret_cast<void*>(&tmp_calc));
	tmp_calc.GetGCSpline();

	// Fit negative binomial to downsampled sites with fixed biases to estimate dispersion parameters
	tmp_calc.optimizer_nbinom_.set_max_objective(BiasCalculationVectors::LogLikelihoodNbinom, reinterpret_cast<void*>(&tmp_calc));

	for(auto k=0; k < tmp_calc.gc_spline_pars_.size(); ++k){
		tmp_calc.fit_pars_.push_back( tmp_calc.gc_spline_pars_.at(k) );
	}
	tmp_calc.fit_pars_.push_back( 1.0 ); // alpha
	tmp_calc.fit_pars_.push_back( 1.0 ); // beta

	tmp_calc.func_calls_ = 0;

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), lower_bound);
	if(BiasCalculationVectors::kSurMult){
		tmp_calc.DeactivateZeroCounts(tmp_calc.bounds_, 0.0, 0);
		for(uintNumFits k=tmp_calc.sur_count_.size()+2; k < tmp_calc.fit_pars_.size()-2; ++k){
			tmp_calc.bounds_.at(k) = BiasCalculationVectors::kLowerBound;
		}
	}
	else{
		tmp_calc.bounds_.at(tmp_calc.sur_count_.size()+1) = BiasCalculationVectors::kBaseValue;
	}
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 1) = BiasCalculationVectors::kBaseValue;
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 2) = BiasCalculationVectors::kBaseValue;
	tmp_calc.optimizer_nbinom_.set_lower_bounds(tmp_calc.bounds_);

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), upper_bound);
	if(BiasCalculationVectors::kSurMult){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else{
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::kLowerBound, 1);
	}
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 1) = BiasCalculationVectors::kUpperBound;
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 2) = BiasCalculationVectors::kUpperBound;
	tmp_calc.optimizer_nbinom_.set_upper_bounds(tmp_calc.bounds_);

	tmp_calc.converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == tmp_calc.optimizer_nbinom_.optimize(tmp_calc.fit_pars_, tmp_calc.loglike_) ){
			//printWarn << "Optimizer reached maximum number of likelihood calculations for seq " << ref_seq_bin << " and length " << insert_length << std::endl;
			tmp_calc.func_calls_ += 10000;
			tmp_calc.converged_ = false;
		}
	}
	catch(const std::exception& e){
		//printWarn << "Nbinom optimizer crashed for seq " << ref_seq_bin << " and length " << insert_length << " with: " << e.what() << std::endl;
		tmp_calc.func_calls_ += 20000;
		tmp_calc.converged_ = false;
	}

	tmp_calc.PostProcessFit();
}

void FragmentDistributionStats::AcquireBiases(const BiasCalculationVectors &calc, mutex &print_mutex){
	if(calc.converged_){
		uintNumFits cur_res = current_bias_result_++;

		// Store biases in vectors from where a single value will be chosen later
		for(uintPercent gc=101; gc--; ){
			tmp_gc_bias_.at(gc).at(cur_res).first = calc.gc_bias_.at(gc);
			tmp_gc_bias_.at(gc).at(cur_res).second = calc.gc_weights_.at(gc);
			//tmp_gc_bias_.at(gc).at(cur_res).second = 1.0;
		}
		for(uintNumFits sur=calc.sur_bias_.size(); sur--; ){
			tmp_sur_bias_.at(sur).at(cur_res) = calc.sur_bias_.at(sur);
		}
		for(uintNumFits par=calc.dispersion_.size(); par--; ){
			tmp_dispersion_parameters_.at(par).at(cur_res) = calc.dispersion_.at(par);
		}
	}

	uintNumFits pars_fitted = ++params_fitted_;

	if( tmp_gc_bias_.at(0).size() > 20 && !(pars_fitted%(tmp_gc_bias_.at(0).size()/20)) ){
		lock_guard<mutex> lock(print_mutex);
		printInfo << "Finished " << static_cast<uintPercentPrint>(Percent(pars_fitted, static_cast<uintNumFits>(tmp_gc_bias_.at(0).size()))) << "% of the bias fits." << std::endl;
	}
}

bool FragmentDistributionStats::StoreBias(){
	// GC bias
	// Calculate weighted median
	for(uintPercent gc=0; gc<tmp_gc_bias_.size(); ++gc){
		tmp_gc_bias_.at(gc).resize(current_bias_result_);

		sort(tmp_gc_bias_.at(gc).begin(), tmp_gc_bias_.at(gc).end());
		double sum_of_weights(0.0);
		for(auto &bias_values : tmp_gc_bias_.at(gc)){
			sum_of_weights += bias_values.second;
		}
		sum_of_weights /= 2.0;

		uintNumFits num_fit=0;
		double sum = tmp_gc_bias_.at(gc).at(0).second;
		while(sum < sum_of_weights){
			sum += tmp_gc_bias_.at(gc).at(++num_fit).second;
		}

		gc_fragment_content_bias_[gc] = tmp_gc_bias_.at(gc).at( num_fit ).first;

		tmp_gc_bias_.at(gc).clear();
		tmp_gc_bias_.at(gc).shrink_to_fit();
	}

	// Calculate Max
	double max_val(*max_element(gc_fragment_content_bias_.begin(), gc_fragment_content_bias_.end()));
	// Normalize so that maximum bias is 1.0
	for(uintPercent gc=gc_fragment_content_bias_.size(); gc--;){
		gc_fragment_content_bias_[gc] /= max_val;
	}

	// Surrounding bias
	// Calculate Median
	std::array<double, 4*Surrounding::Length()> sur_median;
	for(uintNumFits sur=0; sur<tmp_sur_bias_.size(); ++sur){
		tmp_sur_bias_.at(sur).resize(current_bias_result_);

		sort(tmp_sur_bias_.at(sur).begin(), tmp_sur_bias_.at(sur).end());
		sur_median.at(sur) = tmp_sur_bias_.at(sur).at( tmp_sur_bias_.at(sur).size()/2 );

		tmp_sur_bias_.at(sur).clear();
		tmp_sur_bias_.at(sur).shrink_to_fit();
	}

	fragment_surroundings_bias_.CombinePositions(sur_median);

	// Dispersion parameters
	// Calculate Median
	for(uintNumFits par=0; par<tmp_dispersion_parameters_.size(); ++par){
		tmp_dispersion_parameters_.at(par).resize(current_bias_result_);

		sort(tmp_dispersion_parameters_.at(par).begin(), tmp_dispersion_parameters_.at(par).end());
		dispersion_parameters_.at(par) = tmp_dispersion_parameters_.at(par).at( tmp_dispersion_parameters_.at(par).size()/2 );

		tmp_dispersion_parameters_.at(par).clear();
		tmp_dispersion_parameters_.at(par).shrink_to_fit();
	}

	return true;
}

double FragmentDistributionStats::MedianFragmentCoverage(const Reference &ref) const{
	// Calculate weighted median coverage of fragments, so normal coverage divided by 2*average_read_length
	uintFragCount total_filtered_sites(0);
	vector<pair<double, uintFragCount>> frag_coverage;
	frag_coverage.reserve(ref.NumberSequences());
	for(uintRefSeqId ref_seq = 0; ref_seq < ref.NumberSequences(); ++ref_seq){
		if(corrected_abundance_.at(ref_seq)){
			total_filtered_sites += filtered_sites_.at(ref_seq);
			frag_coverage.emplace_back(static_cast<double>(corrected_abundance_.at(ref_seq))/ref.SequenceLength(ref_seq) , filtered_sites_.at(ref_seq));
		}
	}
	sort(frag_coverage.begin(), frag_coverage.end());

	if(0 == frag_coverage.size()){
		return 0.0;
	}
	else{
		total_filtered_sites /= 2;
		uintRefSeqId i(0);
		uintFragCount sum_filtered_sites = frag_coverage.at(0).second;
		while(sum_filtered_sites < total_filtered_sites){
			sum_filtered_sites += frag_coverage.at(++i).second;
		}
		return frag_coverage.at(i).first;
	}
}

bool FragmentDistributionStats::CalculateInsertLengthAndRefSeqBias(const Reference &reference, uintNumThreads num_threads){
	// Select the insert_length_sample_positions
	InsertLengthSpline insert_length_spline;
	if( !insert_length_spline.GetSamplePositions(insert_lengths_) ){
		return false;
	}

	// Sum up the other biases for the insert length, reference sequences pairs
	atomic<uintNumFits> current_param(0);
	atomic<uintNumFits> finished_params(0);
	vector<BiasCalculationParams> params;
	FillParams(params, reference, insert_length_spline.sample_positions_);

	std::vector<std::vector<double>> bias_sum;
	bias_sum.resize(abundance_.size());
	for(auto &sum : bias_sum){
		sum.resize(insert_lengths_.to(), 0.0);
	}

	// Run the predetermined parameters in defined number of threads
	mutex print_mutex;
	thread threads[num_threads];
	for(auto i = num_threads; i--; ){
		threads[i] = thread(BiasSumThread, std::cref(*this), std::cref(reference), std::cref(params), std::ref(current_param), std::ref(finished_params), std::ref(bias_sum), std::ref(print_mutex) );
	}
	for(auto i = num_threads; i--; ){
		threads[i].join();
	}

	// Get median_frag_coverage before modifying corrected_abundance_
	double median_frag_coverage = MedianFragmentCoverage(reference);

	// Initialize biases for fit
	uintFragCount sum_frag_len(0);
	insert_lengths_bias_.Clear();
	insert_lengths_bias_[insert_lengths_.to()-1]; // Resize Vect
	for( auto frag_len : insert_length_spline.sample_positions_ ){
		insert_lengths_bias_.at(frag_len) = 1.0;
		sum_frag_len += insert_lengths_.at(frag_len);
	}

	uintFragCount sum_ref_seq(0);
	ref_seq_bias_.clear();
	ref_seq_bias_.resize(abundance_.size(), 0.0);
	for(auto ref_seq=abundance_.size(); ref_seq--; ){
		if(0 == fragment_lengths_used_for_correction_.at(ref_seq) || Divide(uintFragCount(filtered_sites_.at(ref_seq)), uintFragCount(fragment_lengths_used_for_correction_.at(ref_seq))) < kMinSites){
			corrected_abundance_.at(ref_seq) = 0.0; // Do not fit reference sequences, where the corrected abundance is based on a low number of sites, instead later take the median ref_seq_bias
		}

		if( corrected_abundance_.at(ref_seq) ){
			ref_seq_bias_.at(ref_seq) = 1.0;
			sum_ref_seq += corrected_abundance_.at(ref_seq);
		}
	}

	// Expectation maximization
	double normalization( static_cast<double>(sum_frag_len)/sum_ref_seq ); // Due to adapter reads insert length and ref seq do have different sums, so we account for that
	double precision(numeric_limits<double>::max()), sum;
	uintNumFits iteration(0);
	while(precision > kPrecisionAimRefSeqFragLengthFit && iteration < kMaxIterationsRefSeqFragLengthFit){
		++iteration;
		precision = 1.0;

		// Adjust ref_seq_bias_
		for(auto ref_seq=abundance_.size(); ref_seq--; ){
			if(corrected_abundance_.at(ref_seq)){
				sum = 0.0;
				for( auto frag_len : insert_length_spline.sample_positions_ ){
					sum += insert_lengths_bias_.at(frag_len) * bias_sum.at(ref_seq).at(frag_len);
				}
				sum *= ref_seq_bias_.at(ref_seq);

				if(sum > corrected_abundance_.at(ref_seq)){
					SetToMax(precision, sum/corrected_abundance_.at(ref_seq));
				}
				else{
					SetToMax(precision, corrected_abundance_.at(ref_seq)/sum);
				}

				ref_seq_bias_.at(ref_seq) *= corrected_abundance_.at(ref_seq)/sum;
			}
		}

		// Adjust insert_lengths_bias_
		for( auto frag_len : insert_length_spline.sample_positions_ ){
			sum = 0.0;
			for(auto ref_seq=abundance_.size(); ref_seq--; ){
				sum += ref_seq_bias_.at(ref_seq) * bias_sum.at(ref_seq).at(frag_len);
			}
			sum *= insert_lengths_bias_.at(frag_len) * normalization;

			if(sum > insert_lengths_.at(frag_len)){
				SetToMax(precision, sum/insert_lengths_.at(frag_len));
			}
			else{
				SetToMax(precision, insert_lengths_.at(frag_len)/sum);
			}

			insert_lengths_bias_.at(frag_len) *= insert_lengths_.at(frag_len)/sum;
		}
	}

	// Adjust ref_seq_bias_ for the excluded part of the sequences
	for(uintRefSeqId ref_seq = 0; ref_seq < reference.NumberSequences(); ++ref_seq){
		if(0.0 != ref_seq_bias_.at(ref_seq)){
			ref_seq_bias_.at(ref_seq) /= static_cast<double>(reference.SequenceLength(ref_seq)) / (reference.SequenceLength(ref_seq) - reference.SumExcludedBases(ref_seq));
		}
	}

	// Replace the ref_seq_bias and corrected_abundance for non-fitted reference sequences with weighted median ref_seq_bias and the corresponding corrected abundance
	uintFragCount total_filtered_sites(0);
	vector<pair<double, uintFragCount>> ref_bias;
	ref_bias.reserve(reference.NumberSequences());
	for(uintRefSeqId ref_seq = 0; ref_seq < reference.NumberSequences(); ++ref_seq){
		if(0.0 != ref_seq_bias_.at(ref_seq)){
			total_filtered_sites += filtered_sites_.at(ref_seq);
			ref_bias.emplace_back(ref_seq_bias_.at(ref_seq), filtered_sites_.at(ref_seq));
		}
	}
	sort(ref_bias.begin(), ref_bias.end());

	total_filtered_sites /= 2;
	uintRefSeqId i(0);
	uintFragCount sum_filtered_sites = ref_bias.at(0).second;
	while(sum_filtered_sites < total_filtered_sites){
		sum_filtered_sites += ref_bias.at(++i).second;
	}
	double median_ref_bias = ref_bias.at(i).first;

	for(uintRefSeqId ref_seq = 0; ref_seq < reference.NumberSequences(); ++ref_seq){
		if(0 == fragment_lengths_used_for_correction_.at(ref_seq) || Divide(uintFragCount(filtered_sites_.at(ref_seq)), uintFragCount(fragment_lengths_used_for_correction_.at(ref_seq))) < kMinSites){
			ref_seq_bias_.at(ref_seq) = median_ref_bias;

			if( reference.ReferenceSequenceExcluded(ref_seq) ){
				// We did not calculate any biases for this one, so use median coverage
				corrected_abundance_.at(ref_seq) = reference.SequenceLength(ref_seq)*median_frag_coverage;
			}
			else{
				sum = 0.0;
				for(auto frag_len=insert_lengths_.to(); frag_len-- > insert_lengths_.from(); ){
					sum += insert_lengths_bias_.at(frag_len) * bias_sum.at(ref_seq).at(frag_len);
				}
				corrected_abundance_.at(ref_seq) = sum * reference.SequenceLength(ref_seq) / (reference.SequenceLength(ref_seq) - reference.SumExcludedBases(ref_seq)) * ref_seq_bias_.at(ref_seq);
			}
		}
	}

	// Fill missing insert length using a natural spline
	double max_len_bias(0.0);
	for( auto frag_len : insert_length_spline.sample_positions_ ){
		SetToMax(max_len_bias, insert_lengths_bias_.at(frag_len));
	}
	for( auto frag_len : insert_length_spline.sample_positions_ ){
		insert_lengths_bias_.at(frag_len) /= max_len_bias;
	}

	if( !insert_length_spline.FitInsertLengthWithSpline(insert_lengths_bias_, insert_lengths_) ){
		return false;
	}

	// Find maximum (to normalize it to 1.0 instead of normalizing the sum to 1.0 as it is a bias not a probability)
	max_len_bias = 0.0;
	for(auto frag_len = insert_lengths_bias_.to(); frag_len-- > insert_lengths_.from(); ){
		SetToMax(max_len_bias, insert_lengths_bias_.at(frag_len));
	}

	double max_seq_bias(0.0);
	for(auto ref_seq = ref_seq_bias_.size(); ref_seq--; ){
		SetToMax(max_seq_bias, ref_seq_bias_.at(ref_seq));
	}

	// Normalize so that maximum bias is 1.0
	for(auto frag_len = insert_lengths_bias_.to(); frag_len-- > insert_lengths_.from(); ){
		insert_lengths_bias_.at(frag_len) /= max_len_bias;
	}

	for(auto ref_seq = ref_seq_bias_.size(); ref_seq--; ){
		ref_seq_bias_.at(ref_seq) /= max_seq_bias;
	}

	printInfo << "Fitted reference sequence and insert length biases with a precision of " << precision << " needing " << iteration << " iteration" << (iteration == 1?"":"s") << std::endl;

	return true;
}

void FragmentDistributionStats::ReplaceUncertainCorrectedAbundanceWithMedian(const Reference &ref){
	double median_coverage = MedianFragmentCoverage(ref);

	// Replace corrected_abundance for reference sequences where it is based on a low number of sites
	for(uintRefSeqId ref_seq = 0; ref_seq < ref.NumberSequences(); ++ref_seq){
		if(0 == fragment_lengths_used_for_correction_.at(ref_seq) || Divide(uintFragCount(filtered_sites_.at(ref_seq)), uintFragCount(fragment_lengths_used_for_correction_.at(ref_seq))) < kMinSites){
			corrected_abundance_.at(ref_seq) = ref.SequenceLength(ref_seq)*median_coverage;
		}
	}
}

void FragmentDistributionStats::AddNewBiasCalculations(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, mutex &print_mutex, const Reference &reference){
	uintRefSeqBin ref_seq_bin = num_handled_reference_sequence_bins_;
	while( ref_seq_bin < still_needed_ref_bin ){
		// Check if a spot in the queue is free
		uint32_t queue_spot = 0;
		while( queue_spot < kMaxBinsQueuedForBiasCalc && claimed_bias_bins_.at(queue_spot).test_and_set() ){
			++queue_spot;
		}

		if(queue_spot < kMaxBinsQueuedForBiasCalc){
			// Additional reference sequences can be handled
			if( num_handled_reference_sequence_bins_.compare_exchange_strong(ref_seq_bin, ref_seq_bin+1) ){
				// Handle ref_seq_bin
				CheckLowQExclusions(ref_seq_bin, thread.tmp_frag_count_, reference);
				SortFragmentSites(ref_seq_bin, thread.num_sites_per_insert_length_);
				UpdateBiasCalculationParams(ref_seq_bin, queue_spot, thread.bias_calc_tmp_params_, print_mutex);
				ref_seq_bin = num_handled_reference_sequence_bins_; // Set ref_seq_bin after everything has been done to the new value (In case compare_exchange_strong fails this is done automatically)
			}
			else{
				claimed_bias_bins_.at(queue_spot).clear(); // Let go of claim as reference sequence bin was already handled by another thread
			}
		}
		else{
			break; // End loop as currently no queue spots are available
		}
	}
}

void FragmentDistributionStats::ExecuteBiasCalculations( const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, mutex &print_mutex ){
	auto queue_bin_size = bias_calc_params_.size()/kMaxBinsQueuedForBiasCalc;
	for(uint16_t queue_bin=0; queue_bin < kMaxBinsQueuedForBiasCalc; ++queue_bin ){
		auto cur_par = queue_bin*queue_bin_size + current_bias_param_.at(queue_bin)++;

		for(; cur_par < queue_bin*queue_bin_size + queue_bin_size; cur_par = queue_bin*queue_bin_size + current_bias_param_.at(queue_bin)++){
			--params_left_for_calculation_;

			// Execute calculations
			if( bias_calc_params_.at(cur_par).fragment_length_ ){
				if( bias_calc_params_.at(cur_par).bias_calculation_ ){
					CalculateBiasByBin(thread_values, reference, duplications, bias_calc_params_.at(cur_par).ref_seq_bin_, bias_calc_params_.at(cur_par).fragment_length_);
					AcquireBiases(thread_values, print_mutex);
				}
				else{
					CountDuplicates(duplications, bias_calc_params_.at(cur_par), reference);
				}
			}

			// Check if all calculations in queue bin are finished
			if( ++finished_bias_calcs_.at(queue_bin) == queue_bin_size ){
				// Free memory that is not needed anymore
				fragment_sites_by_ref_seq_bin_by_insert_length_.at( bias_calc_params_.at(cur_par).ref_seq_bin_ ).clear();
				fragment_sites_by_ref_seq_bin_by_insert_length_.at( bias_calc_params_.at(cur_par).ref_seq_bin_ ).shrink_to_fit();

				// All finished, allow to reuse queue bin
				claimed_bias_bins_.at(queue_bin).clear();
			}
		}
	}
}

void FragmentDistributionStats::HandleReferenceSequencesUntil(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, mutex &print_mutex){
	AddNewBiasCalculations(still_needed_ref_bin, thread, print_mutex, reference);
	ExecuteBiasCalculations( reference, duplications, thread.bias_calc_vects_, print_mutex );
}

void FragmentDistributionStats::BiasSumThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<uintNumFits> &current_param,
		atomic<uintNumFits> &finished_params,
		vector<vector<double>> &bias_sum,
		std::mutex &print_mutex ){
	decltype(params.size()) cur_par(current_param++);

	for(; cur_par < params.size(); cur_par = current_param++){
		bias_sum.at( params.at(cur_par).ref_seq_id ).at( params.at(cur_par).fragment_length ) = reference.SumBias(params.at(cur_par).ref_seq_id, params.at(cur_par).fragment_length, 1.0, self.gc_fragment_content_bias_, self.fragment_surroundings_bias_);

		if( params.size() > 20 && !(++finished_params%(params.size()/20)) ){
			lock_guard<mutex> lock(print_mutex);
			printInfo << "Finished " << static_cast<uintPercentPrint>(Percent(finished_params, static_cast<uintNumFits>(params.size()))) << "% of the bias sums." << std::endl;
		}
	}
}

reseq::uintRefSeqId FragmentDistributionStats::SplitCoverageGroups(vector<uintRefSeqId> &coverage_groups) const{
	vector<pair<double,uintRefSeqId>> sorted_ref_seq_bias;
	sorted_ref_seq_bias.reserve(ref_seq_bias_.size());
	for(auto ref_seq = ref_seq_bias_.size(); ref_seq--; ){
		sorted_ref_seq_bias.emplace_back(ref_seq_bias_.at(ref_seq), ref_seq);
	}
	sort(sorted_ref_seq_bias.begin(), sorted_ref_seq_bias.end());

	// Split all ref seqs in groups where the bias is at max a factor two apart
	coverage_groups.clear();
	coverage_groups.resize(ref_seq_bias_.size());
	double group_start = sorted_ref_seq_bias.front().first;
	uintRefSeqId group(0);
	for( auto &bias : sorted_ref_seq_bias ){
		if(bias.first > 2*group_start){
			group_start = bias.first;
			++group;
		}

		coverage_groups.at(bias.second) = group;
	}

	return group+1;
}

void FragmentDistributionStats::BiasNormalizationThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<uintNumFits> &current_param,
		vector<double> &norm,
		mutex &result_mutex,
		const vector<uintRefSeqId> &coverage_groups,
		vector<vector<double>> &max_bias){
	decltype(params.size()) cur_par(current_param++);

	vector<double> tmp_norm;
	tmp_norm.resize(norm.size(), 0.0);
	vector<vector<double>> tmp_max_bias;
	tmp_max_bias.resize(max_bias.size());
	for( auto cov_group = max_bias.size(); cov_group--; ){
		tmp_max_bias.at(cov_group).resize( max_bias.at(cov_group).size(), 0.0 );
	}

	for(; cur_par < params.size(); cur_par = current_param++){
		tmp_norm.at( params.at(cur_par).fragment_length ) += reference.SumBias(tmp_max_bias.at( coverage_groups.at(params.at(cur_par).ref_seq_id) ).at(params.at(cur_par).fragment_length), params.at(cur_par).ref_seq_id, params.at(cur_par).fragment_length, self.ref_seq_bias_.at(params.at(cur_par).ref_seq_id)*self.insert_lengths_bias_.at(params.at(cur_par).fragment_length), self.gc_fragment_content_bias_, self.fragment_surroundings_bias_);
	}

	result_mutex.lock();
	for( auto cov_group = max_bias.size(); cov_group--; ){
		for( auto frag_len = max_bias.at(cov_group).size(); frag_len--; ){
			SetToMax(max_bias.at(cov_group).at(frag_len), tmp_max_bias.at(cov_group).at(frag_len));
		}
	}
	for( auto frag_len = norm.size(); frag_len--; ){
		norm.at(frag_len) += tmp_norm.at(frag_len);
	}
	result_mutex.unlock();
}

double FragmentDistributionStats::CalculateNonZeroThreshold(double bias_normalization, double max_bias) const{
	// Threshold that random number from fragment generation has to reach or fragment count will be zero and does not have to be calculated
	double max_dispersion( Dispersion(bias_normalization*max_bias) );

	return pow(max_dispersion/(max_dispersion+bias_normalization*max_bias), max_dispersion); // Lowest possible F(0) for negative binomial
}

void FragmentDistributionStats::SetUniformBias(){
	ref_seq_bias_.clear();
	ref_seq_bias_.resize(abundance_.size(), 1.0);

	gc_fragment_content_bias_.Clear();
	gc_fragment_content_bias_.Set(0,101,1.0);

	fragment_surroundings_bias_.SetUniform();

	// For the insert_lengths_bias_ we follow the insert length distribution, because uniform is most likely not what people want here
	insert_lengths_bias_.Clear();
	if(insert_lengths_.size()){
		auto max = insert_lengths_.Max();
		for(auto len = insert_lengths_.to(); len-- > insert_lengths_.from(); ){
			insert_lengths_bias_[len] = static_cast<double>(insert_lengths_.at(len)) / max;
		}
		insert_lengths_bias_.Shrink();
	}

	// Approximating poisson here
	dispersion_parameters_.at(0) = 0.0;
	dispersion_parameters_.at(1) = 1e-100;
}

reseq::uintRefSeqBin FragmentDistributionStats::CreateRefBins( const Reference &ref, uintSeqLen max_ref_seq_bin_size ){
	max_ref_seq_bin_length_ = max_ref_seq_bin_size;

	ref_seq_start_bin_.resize(ref.NumberSequences()+1); // Add one to give number of bins at the end
	uintRefSeqBin start_id = 0;
	ref_seq_start_bin_.at(0) = 0;
	for(uintRefSeqId ref_seq=1; ref_seq < ref_seq_start_bin_.size(); ++ref_seq){
		start_id += DivideAndCeil(ref.SequenceLength(ref_seq-1)-ref.SumExcludedBases(ref_seq-1), max_ref_seq_bin_length_);
		ref_seq_start_bin_.at(ref_seq) = start_id;
	}

	ref_seq_bin_def_.resize(ref_seq_start_bin_.back());
	uintRefSeqBin cur_start_bin = 0;
	uintSeqLen cur_exclusion_region;
	for(uintRefSeqId ref_seq=0; ref_seq < ref_seq_start_bin_.size()-1; ++ref_seq){
		auto num_bins = DivideAndCeil(ref.SequenceLength(ref_seq)-ref.SumExcludedBases(ref_seq), max_ref_seq_bin_length_);
		if(num_bins){
			ref_seq_bin_def_.at(cur_start_bin).first = ref_seq;
			ref_seq_bin_def_.at(cur_start_bin).second = 0;
			cur_exclusion_region = 0;
		}

		for( uintRefSeqBin cur_bin=1; cur_bin < num_bins; ++cur_bin){
			auto start_pos = ref_seq_bin_def_.at(cur_start_bin+cur_bin-1).second + RefSeqSplitLength(ref_seq, ref);
			ref.CorrectPositionAddingExcludedRegions(start_pos, cur_exclusion_region, ref_seq);

			ref_seq_bin_def_.at(cur_start_bin+cur_bin).first = ref_seq;
			ref_seq_bin_def_.at(cur_start_bin+cur_bin).second = start_pos;
		}
		cur_start_bin += num_bins;
	}

	return ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1);
}

void FragmentDistributionStats::Prepare( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const vector<uintFragCount> &reads_per_frag_len_bin, const vector<uintFragCount> &lowq_reads_per_frag_len_bin ){
	tmp_abundance_.resize(ref.NumberSequences());
	tmp_insert_lengths_.resize(maximum_insert_length+1);
	tmp_gc_fragment_content_.resize(101);

	for(auto strand=2; strand--; ){
		for(auto nuc=4; nuc--; ){
			tmp_outskirt_content_.at(strand).at(nuc).resize(2*kOutskirtRange);
		}
	}

	PrepareBiasCalculation( ref, maximum_insert_length, max_ref_seq_bin_size, reads_per_frag_len_bin, lowq_reads_per_frag_len_bin );
}

void FragmentDistributionStats::FillInOutskirtContent( const Reference &reference, const BamAlignmentRecord &record_start, uintSeqLen fragment_start_pos, uintSeqLen fragment_end_pos ){
	bool reversed_fragment = hasFlagLast(record_start); // record is forward, so if it is sequenced last, fragment is reversed
	auto ref_pos = fragment_start_pos;
	uintSeqLen outskirt_pos;
	auto &ref_seq = reference.ReferenceSequence(record_start.rID);

	// Fill in fragment start
	Surrounding surroundings;
	reference.ForwardSurroundingWithN(surroundings, record_start.rID, ref_pos);
	tmp_fragment_surroundings_.Count(surroundings);

	if(reversed_fragment){
		// Reverse strand: Fill in postfragment-content
		if( ref_pos < kOutskirtRange ){
			outskirt_pos = kOutskirtRange + ref_pos;
			ref_pos = 0;
		}
		else{
			outskirt_pos = 2*kOutskirtRange;
			ref_pos -= kOutskirtRange;
		}

		for(; outskirt_pos-- > kOutskirtRange; ){
			if( !IsN(ref_seq, ref_pos) ){
				++tmp_outskirt_content_.at(true).at(Complement::Dna5(at(ref_seq, ref_pos))).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
	else{
		// Forward strand: Fill in prefragment-content
		if( ref_pos < kOutskirtRange ){
			outskirt_pos = kOutskirtRange - ref_pos;
			ref_pos = 0;
		}
		else{
			outskirt_pos = 0;
			ref_pos -= kOutskirtRange;
		}

		for(; outskirt_pos < kOutskirtRange; ++outskirt_pos){
			if( !IsN(ref_seq, ref_pos) ){
				++tmp_outskirt_content_.at(false).at(at(ref_seq, ref_pos)).at(outskirt_pos);
			}
			++ref_pos;
		}
	}

	ref_pos = fragment_end_pos;

	// Fill in fragment end
	reference.ReverseSurroundingWithN(surroundings, record_start.rNextId, ref_pos-1);
	tmp_fragment_surroundings_.Count(surroundings);

	if(reversed_fragment){
		// Reverse strand: Fill in prefragment-content
		outskirt_pos = kOutskirtRange;
		for(; outskirt_pos-- && ref_pos < reference.SequenceLength(record_start.rNextId);){
			if( !IsN(ref_seq, ref_pos) ){
				++tmp_outskirt_content_.at(true).at(Complement::Dna5(at(ref_seq, ref_pos))).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
	else{
		// Forward strand:Fill in postfragment-content
		for(; outskirt_pos < 2*kOutskirtRange && ref_pos < reference.SequenceLength(record_start.rNextId); ++outskirt_pos){
			if( !IsN(ref_seq, ref_pos) ){
				++tmp_outskirt_content_.at(false).at(at(ref_seq, ref_pos)).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
}

void FragmentDistributionStats::HandleReferenceSequencesUntil(uintRefSeqId still_needed_reference_sequence, uintSeqLen still_needed_position, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, mutex &print_mutex){
	// Check if new biases can be added for calculation
	auto still_needed_ref_bin = ref_seq_start_bin_.at(still_needed_reference_sequence);
	while(still_needed_ref_bin+1 < ref_seq_start_bin_.at(still_needed_reference_sequence+1) && still_needed_position > ref_seq_bin_def_[still_needed_ref_bin+1].second){
		// Not the last bin of the reference and next bin starts before still_needed_position
		++still_needed_ref_bin;
	}

	HandleReferenceSequencesUntil(still_needed_ref_bin, thread, reference, duplications, print_mutex);
}

void FragmentDistributionStats::FinishThreads(ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, mutex &print_mutex){
	while(true){
		HandleReferenceSequencesUntil(ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1), thread, reference, duplications, print_mutex);

		if(0 == params_left_for_calculation_){
			break;
		}
		else{
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}
}

void FragmentDistributionStats::Finalize(){
	// Acquire data from temporary vectors
	Acquire(abundance_, tmp_abundance_);
	insert_lengths_.Acquire(tmp_insert_lengths_);
	gc_fragment_content_.Acquire(tmp_gc_fragment_content_);

	fragment_surroundings_.Acquire(tmp_fragment_surroundings_);

	for(auto strand=2; strand--; ){
		for(auto nuc=4; nuc--; ){
			outskirt_content_.at(strand).at(nuc).Acquire(tmp_outskirt_content_.at(strand).at(nuc));
		}
	}
}

bool FragmentDistributionStats::FinalizeBiasCalculation(const Reference &reference, uintNumThreads num_threads, FragmentDuplicationStats &duplications){
	if(calculate_bias_){
		if( 0 == current_bias_result_ ){
			printWarn << "No bias fit converged. Continuing with uniform coverage." << std::endl;
			SetUniformBias();
			return true;
		}
		else{
			printInfo << current_bias_result_ << " of " << tmp_gc_bias_.at(0).size() << " bias fits [" << static_cast<uintPercentPrint>(Percent(current_bias_result_, static_cast<uintNumFits>(tmp_gc_bias_.at(0).size()))) << "%] converged." << std::endl;

			if( !StoreBias() ){
				return false;
			}

			if( !CalculateInsertLengthAndRefSeqBias(reference, num_threads) ){
				return false;
			}
		}
	}
	else{
		printInfo << "No bias calculation done. Storing uniform coverage." << std::endl;
		SetUniformBias();
		ReplaceUncertainCorrectedAbundanceWithMedian(reference);
	}

	// Free not needed memory at the end
	fragment_sites_by_ref_seq_bin_.clear();
	fragment_sites_by_ref_seq_bin_.shrink_to_fit();
	fragment_sites_by_ref_seq_bin_cur_id_.clear();
	fragment_sites_by_ref_seq_bin_cur_id_.shrink_to_fit();
	fragment_sites_by_ref_seq_bin_by_insert_length_.clear();
	fragment_sites_by_ref_seq_bin_by_insert_length_.shrink_to_fit();

	lowq_site_start_by_ref_seq_bin_.clear();
	lowq_site_start_by_ref_seq_bin_.shrink_to_fit();
	lowq_site_start_by_ref_seq_bin_cur_id_.clear();
	lowq_site_start_by_ref_seq_bin_cur_id_.shrink_to_fit();
	lowq_site_end_by_ref_seq_bin_.clear();
	lowq_site_end_by_ref_seq_bin_.shrink_to_fit();
	lowq_site_end_by_ref_seq_bin_cur_id_.clear();
	lowq_site_end_by_ref_seq_bin_cur_id_.shrink_to_fit();
	lowq_site_start_exclusion_.clear();
	lowq_site_start_exclusion_.shrink_to_fit();
	lowq_site_end_exclusion_.clear();
	lowq_site_end_exclusion_.shrink_to_fit();

	filtered_sites_.clear();
	filtered_sites_.shrink_to_fit();
	fragment_lengths_used_for_correction_.clear();
	fragment_lengths_used_for_correction_.shrink_to_fit();

	ref_seq_in_nxx_.clear();
	ref_seq_in_nxx_.shrink_to_fit();
	ref_seq_start_bin_.clear();
	ref_seq_start_bin_.shrink_to_fit();
	ref_seq_bin_def_.clear();
	ref_seq_bin_def_.shrink_to_fit();

	duplications.FinalizeDuplicationVector();

	if(calculate_bias_){
		auto ref_size = reference.TotalSize();
		printInfo << "Finished bias calculation excluding " << excluded_lowq_positions_ << " of " << ref_size << " [" << static_cast<uintPercentPrint>(Percent(excluded_lowq_positions_,ref_size)) << "%] reference positions within " << excluded_lowq_regions_ << " regions due to low quality mappings" << std::endl;
	}
	else{
		printInfo << "Finished collecting duplicates" << std::endl;
	}

	return true;
}

double FragmentDistributionStats::CorrectedCoverage(const Reference &ref, uintReadLen average_read_len){
	// Calculate corrected coverage
	double corrected_abundance(0.0);
	for(uintRefSeqId ref_seq = 0; ref_seq < ref.NumberSequences(); ++ref_seq){
		corrected_abundance += corrected_abundance_.at(ref_seq);
	}

	corrected_abundance_.clear();
	corrected_abundance_.shrink_to_fit();

	return corrected_abundance * average_read_len * 2 / ref.TotalSize();
}

bool FragmentDistributionStats::UpdateRefSeqBias(RefSeqBiasSimulation model, const std::string &bias_file, const Reference &ref, mt19937_64 &rgen){
	switch(model){
	case kKeep:
		if(ref.NumberSequences() == ref_seq_bias_.size()){
			break;
		}
		else{
			printWarn << "Cannot keep reference biases, because they don't match the reference. Will remove all reference biases. Reference has " << ref.NumberSequences() << " and there are " << ref_seq_bias_.size() << " biases stored." << std::endl;
		}
	case kNo:
		ref_seq_bias_.clear();
		ref_seq_bias_.resize(ref.NumberSequences(), 1.0);
		break;
	case kDraw:
	{
		// Temporarily store the old one to draw from
		std::vector<double> old_bias_;
		old_bias_.resize( ref_seq_bias_.size() );
		for(auto seq=ref_seq_bias_.size(); seq--; ){
			old_bias_.at(seq) = ref_seq_bias_.at(seq);
		}

		// Clear out the old bias
		ref_seq_bias_.clear();
		ref_seq_bias_.resize(ref.NumberSequences());

		// Draw new biases from old biases with replacement
		uniform_int_distribution<uintRefSeqId> rdist(0,old_bias_.size());
		for(auto seq=ref_seq_bias_.size(); seq--; ){
			ref_seq_bias_.at(seq) = old_bias_.at( rdist(rgen) );
		}

		break;
	}
	case kFile:
	{
		// Clear out the old bias
		ref_seq_bias_.clear();
		ref_seq_bias_.resize(ref.NumberSequences(), 0.0);

		ifstream fbias(bias_file);
		if( fbias.is_open() ){
			uintErrorCount errors = 0;

			// Check which reference sequences are needed and store their ref_seq_bias_ position
			unordered_map<string, uintRefSeqId> ref_seq_ids;
			string tmp_id;
			for(uintRefSeqId seq=0; seq < ref.NumberSequences(); ++seq){
				tmp_id = toCString(ref.ReferenceId(seq));
				ref_seq_ids.emplace(tmp_id.substr(0,tmp_id.find(' ')), seq);
			}

			// Read in information from file (One line corresponds to one reference sequence [identifier bias])
			double bias(0.0);
			string line;
			uintRefSeqId nline(0);
			bool empty_line = false; // Allow only a single empty line at the end of the file
			std::vector<bool> bias_found;
			bias_found.resize(ref_seq_bias_.size(), false);
			while( getline(fbias,line) ){
				if(empty_line){
					if( errors < kMaxErrorsShownPerFile){
						printErr << "Error reading in reference sequence biases. Line " << nline << " is empty in file " << bias_file << std::endl;
					}
					IncreaseErrorCounter(errors);
				}
				else{
					++nline;

					if(0 == line.size()){
						empty_line = true;
					}
					else{
						auto bias_sep = line.find_last_of(" \t");
						if( bias_sep == string::npos ){
							if( errors < kMaxErrorsShownPerFile){
								printErr << "Error reading in reference sequence biases. Line " << nline << " does not have a field separator in file " << bias_file << std::endl;
							}
							IncreaseErrorCounter(errors);
						}
						else{
							try{
								bias = stod( line.substr(bias_sep+1) );
							}
							catch(const exception &e){
								if( errors < kMaxErrorsShownPerFile){
									printErr << "Could not convert bias starting at postion " << bias_sep+1 << " in line " << nline << " in file " << bias_file << " into double\nException thrown: " << e.what() << std::endl;
								}
								IncreaseErrorCounter(errors);
								bias = 0.0;
							}

							if(0.0 > bias){
								if( errors < kMaxErrorsShownPerFile){
									printErr << "Bias in line " << nline << " in file " << bias_file << " is negative: " << bias << std::endl;
								}
								IncreaseErrorCounter(errors);
							}

							auto id_len = line.find(' ');
							if( id_len == string::npos ){
								id_len = bias_sep; // No additional sequence description included, so just take the identifier
							}

							uint16_t id_start = 0;
							if('>' == line.at(0)){
								++id_start;
								--id_len;
							}

							// Find reference sequence id and insert bias at that position
							auto it = ref_seq_ids.find( line.substr(id_start, id_len) );
							if (it != ref_seq_ids.end()){
								ref_seq_bias_.at( it->second ) = bias;
								bias_found.at( it->second ) = true;
							}
						}
					}
				}
			}
			fbias.close();

			for(uintRefSeqId seq=0; seq < ref.NumberSequences(); ++seq){
				if(!bias_found.at(seq)){
					if( errors < kMaxErrorsShownPerFile){
						printErr << "Could not find bias for reference sequence " << ref.ReferenceId(seq) << " in file " << bias_file << std::endl;
					}
					IncreaseErrorCounter(errors);
				}
			}

			if(errors){
				return false;
			}

			printInfo << "Read in " << ref_seq_bias_.size() << " biases from " << nline-(empty_line?1:0) << " biases in file " << bias_file << std::endl;
		}
		else{
			printErr << "Unable to open reference bias file " << bias_file << std::endl;
			return false;
		}

		break;
	}
	default:
		printErr << "Unknown option chosen for reference sequence bias" << std::endl;
		return false;
	}

	return true;
}

double FragmentDistributionStats::CalculateBiasNormalization(vector<uintRefSeqId> &coverage_groups, vector<vector<double>> &non_zero_thresholds, const Reference &reference, uintNumThreads num_threads, uintFragCount total_reads) const{
	// Select the insert_length_sample_positions
	InsertLengthSpline insert_length_spline;
	if( !insert_length_spline.GetSamplePositions(insert_lengths_) ){
		return 0.0;
	}

	atomic<uintNumFits> current_param(0);
	vector<BiasCalculationParams> params;
	FillParamsSimulation(params, reference, insert_length_spline.sample_positions_);

	auto num_groups = SplitCoverageGroups(coverage_groups);
	non_zero_thresholds.clear();
	non_zero_thresholds.resize(num_groups);
	for( auto &thresholds : non_zero_thresholds ){
		thresholds.resize(insert_lengths_.to(), 0.0); // First use for maximum bias: non_zero_thresholds[coverageGroup][insertLength] = maxBias;
	}

	mutex result_mutex;
	vector<double> normalization_by_frag_len;
	normalization_by_frag_len.resize(insert_lengths_.to(), 0.0);

	// Run the predetermined parameters in defined number of threads
	thread threads[num_threads];
	for(auto i = num_threads; i--; ){
		threads[i] = thread(BiasNormalizationThread, std::cref(*this), std::cref(reference), std::cref(params), std::ref(current_param), std::ref(normalization_by_frag_len), std::ref(result_mutex), std::cref(coverage_groups), std::ref(non_zero_thresholds) );
	}
	for(auto i = num_threads; i--; ){
		threads[i].join();
	}

	// Estimate not counted fragment length using a natural spline
	if( !insert_length_spline.InterpolateNormalizationWithSpline(normalization_by_frag_len, insert_lengths_bias_) ){
		return 0.0;
	}

	for( auto &group : non_zero_thresholds ){
		// Get maximum ratio
		double max_ratio = 0.0; // max_bias/frag_bias
		for(uintSeqLen s=0; s < insert_length_spline.sample_positions_.size(); ++s){
			auto frag_len = insert_length_spline.sample_positions_.at(s);
			double ratio = group.at(frag_len) / insert_lengths_bias_.at(frag_len);
			SetToMax(max_ratio, ratio);
		}

		// Use maximum ratio to set not counted values
		for(uintSeqLen s=1; s < insert_length_spline.sample_positions_.size(); ++s){
			for(auto frag_len = insert_length_spline.sample_positions_.at(s-1)+1; frag_len < insert_length_spline.sample_positions_.at(s); ++frag_len){
				group.at(frag_len) = max_ratio * insert_lengths_bias_.at(frag_len);
			}
		}

		for(auto frag_len = insert_length_spline.sample_positions_.at( insert_length_spline.sample_positions_.size()-1 )+1; frag_len < group.size(); ++frag_len){
			group.at(frag_len) = max_ratio * insert_lengths_bias_.at(frag_len);
		}
	}

	// Calculate normalization and threshold
	double normalization(0.0);
	for(auto norm : normalization_by_frag_len){
		normalization += norm;
	}
	double full_normalization = total_reads / (normalization * 2) / reference.NumAlleles(); // Bias is equal for both strands so it is only calculated for one, but the number of fragments at a site is calculated separately per strand: This means we have to half the calculated number of reads

	for( auto &thresholds : non_zero_thresholds ){
		for( auto &thresh : thresholds){
			if(0.0 == thresh){
				thresh = 1.0; // If the bias is zero we cannot get fragment counts drawn for this
			}
			else{
				thresh = CalculateNonZeroThreshold(full_normalization, thresh);
			}
		}
	}

	return full_normalization;
}

reseq::uintDupCount FragmentDistributionStats::NegativeBinomial(double p, double r, double probability_chosen){
	// (1-p)^r * mult_{i=1}^{k}[p * (r+i-1)/i]
	double probability_count( pow(1-p, r) ), probability_sum(probability_count);
	uintDupCount count(0);

	while(probability_sum < probability_chosen){
		probability_count *= p * ((r-1) / ++count + 1);
		probability_sum += probability_count;
	}

	return count;
}

reseq::uintDupCount FragmentDistributionStats::GetFragmentCounts(double bias_normalization, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintPercent gc, const Surrounding &fragment_start, const Surrounding &fragment_end, double probability_chosen) const{
	double bias( Reference::Bias( ref_seq_bias_.at(ref_seq_id), insert_lengths_bias_[fragment_length], gc_fragment_content_bias_[gc], fragment_surroundings_bias_.Bias(fragment_start), fragment_surroundings_bias_.Bias(fragment_end) ) );
	if(0.0 < bias){
		double mean( bias * bias_normalization );
		double dispersion( Dispersion(mean) );

		return NegativeBinomial(mean/(mean+dispersion), dispersion, probability_chosen);
	}
	else{
		return 0;
	}
}

void FragmentDistributionStats::PreparePlotting(){
	array<double, 4*Surrounding::Length()> separated_bias;
	fragment_surroundings_bias_.SeparatePositions(separated_bias);

	for( auto &sur_vect : fragment_surrounding_bias_by_base_ ){
		sur_vect.resize(Surrounding::Length(), 0.0);
	}

	// Split by nucleotide
	for( uintSurPos i=0; i < separated_bias.size(); ++i){
		fragment_surrounding_bias_by_base_.at(i%4).at(i/4) = separated_bias.at(i);
	}
}

bool FragmentDistributionStats::WriteRefSeqBias(const std::string &bias_file, const Reference &reference){
	if(reference.NumberSequences() != ref_seq_bias_.size()){
		printErr << "Reference and reseq file do not match. The reference has " << reference.NumberSequences() << " sequences and the reseq file has biases for " << ref_seq_bias_.size() << " sequences" << std::endl;
		return false;
	}

	if(bias_file.empty()){
		for( uintRefSeqId ref_seq = 0; ref_seq < ref_seq_bias_.size(); ++ref_seq ){
			std::cout << reference.ReferenceIdFirstPart(ref_seq) << '\t' << ref_seq_bias_.at(ref_seq) << '\n';
		}
	}
	else{
		ofstream fbias(bias_file);
		if( !fbias.is_open() ){
			printErr << "Unable to open reference bias file " << bias_file << std::endl;
			return false;
		}
		else{
			for( uintRefSeqId ref_seq = 0; ref_seq < ref_seq_bias_.size(); ++ref_seq ){
				fbias << reference.ReferenceIdFirstPart(ref_seq) << '\t' << ref_seq_bias_.at(ref_seq) << '\n';
			}

			fbias.close();
		}
	}

	return true;
}
