#include "FragmentDistributionStats.h"
using reseq::FragmentDistributionStats;
using reseq::BiasCalculationVectors;

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

void BiasCalculationVectors::PrepareSplines(){
	// https://en.wikipedia.org/wiki/Spline_(mathematics) "Computation of Natural Cubic Splines" slightly adapted to get the linear combination of the spline parameters a
	const uintPercent df = kGCSplineDf;
	auto &x = gc_knots_;

	array<double, df> l; // j=[0, n]
	array<array<double, df>, df> c, z; // c[j][a_index]: j=[0, n]
	array<double, df-1> mu, h; // i=[0, n-1]
	auto &b = lin_comb_gc_splines_.at(0); // b[i][a_index]: i=[0, n-1]
	auto &d = lin_comb_gc_splines_.at(2); // d[i][a_index]: i=[0, n-1]
	array<array<double, df>, df-1> beta; //  beta[k][a_index]: k=[1, n-1]

	for(uintPercent i=0; i < h.size(); ++i){
		h.at(i) = x.at(i+1) - x.at(i);
	}

	for(uintPercent k=1; k < beta.size(); ++k){
		beta.at(k).fill(0.0);
		beta.at(k).at(k+1) = 3/h.at(k);
		beta.at(k).at(k) = -3/h.at(k) - 3/h.at(k-1);
		beta.at(k).at(k-1) = 3/h.at(k-1);
	}

	l.at(0) = 0.0;
	mu.at(0) = 0.0;
	z.at(0).fill(0.0);

	for(uintPercent k=1; k < beta.size(); ++k){
		l.at(k) = 2*(x.at(k+1)-x.at(k-1)) - h.at(k-1)*mu.at(k-1);
		mu.at(k) = h.at(k)/l.at(k);
		for(uintPercent ai=0; ai<df; ++ai){
			z.at(k).at(ai) = (beta.at(k).at(ai) - h.at(k-1) * z.at(k-1).at(ai))/l.at(k);
		}
	}

	l.at(l.size()-1) = 1.0;
	z.at(z.size()-1).fill(0.0);
	c.at(c.size()-1).fill(0.0);

	for(uintPercent i=h.size(); i--;){
		for(uintPercent ai=0; ai<df; ++ai){
			c.at(i).at(ai) = z.at(i).at(ai) - mu.at(i)*c.at(i+1).at(ai);
			b.at(i).at(ai) = -h.at(i)*(c.at(i+1).at(ai) + 2*c.at(i).at(ai))/3;
			d.at(i).at(ai) = (c.at(i+1).at(ai) - c.at(i).at(ai))/3/h.at(i);
		}

		b.at(i).at(i+1) += 1/h.at(i);
		b.at(i).at(i) -= 1/h.at(i);
	}

	// Store values from c in lin_comb_gc_splines_ (as c needs to be one larger for the calculation it is not a reference like b and d)
	for(uintPercent i=h.size(); i--;){
		for(uintPercent ai=0; ai<df; ++ai){
			lin_comb_gc_splines_.at(1).at(i).at(ai) = c.at(i).at(ai);
		}
	}
}

void BiasCalculationVectors::GetSplineCoefficients(double &a, double &b, double &c, double &d, uintPercent k, const vector<double> &spline_pars){
	a = spline_pars.at(k+1);
	b = 0.0;
	c = 0.0;
	d = 0.0;

	for(uintPercent ai=1; ai < gc_spline_pars_.size(); ++ai){
		b += spline_pars.at(ai) * lin_comb_gc_splines_.at(0).at(k).at(ai-1);
		c += spline_pars.at(ai) * lin_comb_gc_splines_.at(1).at(k).at(ai-1);
		d += spline_pars.at(ai) * lin_comb_gc_splines_.at(2).at(k).at(ai-1);
	}
}

void BiasCalculationVectors::CalculateSpline(const vector<double> &spline_pars){
	// Constant bias values for gc below any observed sites: Extrapolate the linear function at the lowest measured GC by half a GC unit and use this value
	double norm = spline_pars.at(0);
	double a, b, c, d;
	GetSplineCoefficients(a, b, c, d, 0, spline_pars);
	double lower_const_no_logit = a - 0.5*b;
	double lower_const = norm*InvLogit2(lower_const_no_logit) + kBaseValue;
	for(uintPercent gc=0; gc < gc_knots_.at(0); ++gc){
		gc_bias_no_logit_.at(gc) = lower_const_no_logit;
		gc_bias_.at(gc) = lower_const;
	}

	// Fill observed range of gc with spline values
	gc_bias_no_logit_.at( gc_knots_.at(0) ) = a;
	gc_bias_.at( gc_knots_.at(0) ) = norm*InvLogit2(a) + kBaseValue;
	uintPercent k=0;
	for(uintPercent gc=gc_knots_.at(0)+1; gc < gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		if(gc == gc_knots_.at(k+1)){
			// Reached new spline part
			GetSplineCoefficients(a, b, c, d, ++k, spline_pars);
			gc_bias_no_logit_.at(gc) = a;
			gc_bias_.at(gc) = norm*InvLogit2(a) + kBaseValue;
		}
		else{
			// Continue on current spline part
			uintPercent cur_gc = gc-gc_knots_.at(k);
			gc_bias_no_logit_.at(gc) = a + b*cur_gc + c*cur_gc*cur_gc + d*cur_gc*cur_gc*cur_gc;
			gc_bias_.at(gc) = norm*InvLogit2( gc_bias_no_logit_.at(gc) ) + kBaseValue;
		}
	}
	++k;

	gc_bias_no_logit_.at( gc_knots_.at(k) ) = spline_pars.at(k);
	gc_bias_.at( gc_knots_.at(k) ) = norm*InvLogit2(spline_pars.at(k)) + kBaseValue;

	// Constant bias values for gc above any observed sites: Extrapolate the linear function at the highest measured GC by half a GC unit and use this value
	uintPercent cur_gc = gc_knots_.at(k)-gc_knots_.at(k-1);
	double upper_const_no_logit = gc_spline_pars_.at(k) + 0.5*(b + c*cur_gc); // d=0 since second derivative must be 0 (natural spline)
	double upper_const = norm*InvLogit2( upper_const_no_logit ) + kBaseValue;
	for(uintPercent gc=gc_knots_.at(k)+1; gc < gc_bias_.size(); ++gc){
		gc_bias_no_logit_.at(gc) = upper_const_no_logit;
		gc_bias_.at(gc) = upper_const;
	}
}

void BiasCalculationVectors::CalculateSplineGrad(const vector<double> &spline_pars, vector<double> &grad, uintNumFits grad_offset){
	for( uintPercent gc=0; gc < gc_bias_grad_.size(); ++gc ){
		grad.at(grad_offset) += gc_bias_grad_.at(gc) / spline_pars.at(0);
		gc_bias_grad_.at(gc) *= InvLogit2(-gc_bias_no_logit_.at(gc))/2;
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

	double max_val = 0.0;
	for( uintPercent k=gc_spline_pars_.size()-1; k--; ){
		SetToMax(max_val, gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second);
	}
	max_val /= 1.6;
	gc_spline_pars_.at(0) = max_val;

	max_val *= 2;
	for( uintPercent k=gc_spline_pars_.size()-1; k--; ){
		if( gc_count_.at( gc_knots_.at(k) ) ){
			double val = gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second;

			gc_spline_pars_.at(k+1) = log(val/(max_val-val)); // Set gc bias to set the gradient to zero for this gc

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
	std::vector<double> grad;
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

	double grad_term2 = 2*log(1/(bias/r+1)); // Be careful to handle r->inf without numerical instabilities, should result in r*grad_term2 -> -2*bias
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

inline void FragmentDistributionStats::IncreaseErrorCounter(uintErrorCount &errors){
	if( errors++ == kMaxErrorsShownPerFile ){
		printErr << "Maximum number of errors reached. Additional errors are not shown for this file." << std::endl;
	}
}

void FragmentDistributionStats::PrepareBiasCalculation( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const vector<uintFragCount> &reads_per_frag_len_bin ){
	// Bin reference sequence
	auto num_bins = CreateRefBins(ref, max_ref_seq_bin_size);

	vector<uintFragCount> reads_per_ref_seq_bin;
	reads_per_ref_seq_bin.resize(num_bins);

	uintRefLenCalc ref_seq_start_bin(0);
	uintRefSeqBin ref_seq_bin(0);
	for(uintRefSeqId ref_seq = 0; ref_seq < ref.NumberSequences(); ++ref_seq){
		if( !ref.ReferenceSequenceExcluded(ref_seq) ){
			for(uintSeqLen cur_bin = 0; cur_bin <= ref.SequenceLength(ref_seq)/maximum_insert_length; ++cur_bin){
				ref_seq_bin = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length, ref_seq_bin);
				reads_per_ref_seq_bin.at(ref_seq_bin) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);

				// Add it also to all bins the fragment might overlap, as we don' know the exact start position of the fragment
				if(cur_bin < ref.SequenceLength(ref_seq)/maximum_insert_length){
					auto ref_seq_bin2 = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length - maximum_insert_length, ref_seq_bin);
					if(ref_seq_bin != ref_seq_bin2){
						reads_per_ref_seq_bin.at(ref_seq_bin2) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
					}
				}

				if(cur_bin){
					auto ref_seq_bin2 = GetRefSeqBin(ref_seq, cur_bin*maximum_insert_length + maximum_insert_length, ref_seq_bin);
					if(ref_seq_bin != ref_seq_bin2){
						reads_per_ref_seq_bin.at(ref_seq_bin2) += reads_per_frag_len_bin.at(ref_seq_start_bin+cur_bin);
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

	// Prepare bias fitting
	bias_calc_params_.resize( maximum_insert_length * kMaxBinsQueuedForBiasCalc );

	ref.RefSeqsInNxx(ref_seq_in_nxx_, BiasCalculationVectors::kNXXRefSeqs);

	// Temporary bias fit result vectors
	params_fitted_ = 0;
	auto max_fits = ref.NumRefSeqBinsInNxx(ref_seq_in_nxx_, max_ref_seq_bin_length_) * BiasCalculationVectors::kNumFitsInsertLength;
	for(uintPercent gc=tmp_gc_bias_.size(); gc--;){
		tmp_gc_bias_.at(gc).resize(max_fits, nan("0"));
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
			// Sort them by counts and only take top kNumFitsInsertLength
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
				if(n < BiasCalculationVectors::kNumFitsInsertLength){
					bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, calculate_bias_);
				}
				else{
					bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, false);
				}
			}

			if(tmp_params.size() < BiasCalculationVectors::kNumFitsInsertLength){
				params_fitted_ += BiasCalculationVectors::kNumFitsInsertLength-tmp_params.size(); // Do not have to be fitted, so are already done

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

void FragmentDistributionStats::FillParams(vector<BiasCalculationParams> &params, const Reference &reference) const{
	params.reserve( insert_lengths_.size() * abundance_.size() );

	for( auto frag_length=max(static_cast<uintSeqLen>(1),static_cast<uintSeqLen>(insert_lengths_.from())); frag_length < insert_lengths_.to(); ++frag_length){
		if(insert_lengths_.at(frag_length)){
			for( uintRefSeqId ref_id=abundance_.size(); ref_id--; ){
				if(abundance_.at(ref_id) && !reference.ReferenceSequenceExcluded(ref_id)){
					// Add parameters for later calculation
					params.push_back({ref_id, frag_length});
				}
			}
		}
	}
	params.shrink_to_fit();
}

void FragmentDistributionStats::FillParamsSimulation(vector<BiasCalculationParams> &params) const{
	params.reserve( insert_lengths_.size() * abundance_.size() );

	for( auto frag_length=max(static_cast<uintSeqLen>(1),static_cast<uintSeqLen>(insert_lengths_.from())); frag_length < insert_lengths_.to(); ++frag_length){
		if(insert_lengths_.at(frag_length)){
			for( uintRefSeqId ref_id=abundance_.size(); ref_id--; ){
				if(abundance_.at(ref_id)){
					// Add parameters for later calculation
					params.push_back({ref_id, frag_length});
				}
			}
		}
	}
	params.shrink_to_fit();
}

void FragmentDistributionStats::CountDuplicates(FragmentDuplicationStats &duplications, const BiasCalculationParamsSplitSeqs &params, const Reference &reference){
	duplications.AddDuplicates( fragment_sites_by_ref_seq_bin_by_insert_length_.at(params.ref_seq_bin_).at(params.fragment_length_), GetRefSeqId(params.ref_seq_bin_), params.fragment_length_, ref_seq_bin_def_.at(params.ref_seq_bin_).second, reference  );
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

	reference.GetFragmentSites(tmp_calc.sites_, ref_seq_id, insert_length, ref_start, ref_end);
	AddFragmentsToSites(tmp_calc.sites_, fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length) );
	duplications.AddSites(tmp_calc.sites_, insert_length);
	tmp_calc.GetCounts();
	tmp_calc.RemoveUnnecessarySites();
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

	tmp_calc.optimizer_spline_.set_max_objective(BiasCalculationVectors::LogLikeGcSpline, reinterpret_cast<void*>(&tmp_calc));
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
			tmp_gc_bias_.at(gc).at(cur_res) = calc.gc_bias_.at(gc);
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
	if( 0 == current_bias_result_ ){
		printWarn << "No bias fit converged. Continuing with uniform coverage." << std::endl;
		SetUniformBias();
		return true;
	}
	else{
		printInfo << current_bias_result_ << " of " << tmp_gc_bias_.at(0).size() << " bias fits [" << static_cast<uintPercentPrint>(Percent(current_bias_result_, static_cast<uintNumFits>(tmp_gc_bias_.at(0).size()))) << "%] converged." << std::endl;
	}

	// GC bias
	// Calculate Median
	for(uintPercent gc=0; gc<tmp_gc_bias_.size(); ++gc){
		tmp_gc_bias_.at(gc).resize(current_bias_result_);

		sort(tmp_gc_bias_.at(gc).begin(), tmp_gc_bias_.at(gc).end());
		gc_fragment_content_bias_[gc] = tmp_gc_bias_.at(gc).at( tmp_gc_bias_.at(gc).size()/2 );

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

void FragmentDistributionStats::AddNewBiasCalculations(uintRefSeqBin still_needed_ref_bin, ThreadData &thread, mutex &print_mutex){
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
	AddNewBiasCalculations(still_needed_ref_bin, thread, print_mutex);
	ExecuteBiasCalculations( reference, duplications, thread.bias_calc_vects_, print_mutex );
}

void FragmentDistributionStats::BiasSumThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<uintNumFits> &current_param,
		vector<double> &insert_length_sum,
		vector<double> &ref_seq_sum,
		vector<vector<VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc,
		std::mutex &result_mutex ){
	decltype(params.size()) cur_par(current_param++);

	double sum;
	vector<double> len, ref;
	len.resize(insert_length_sum.size(), 0.0);
	ref.resize(ref_seq_sum.size(), 0.0);

	for(; cur_par < params.size(); cur_par = current_param++){
		sum = reference.SumBias(site_count_by_insert_length_gc.at(params.at(cur_par).fragment_length), params.at(cur_par).ref_seq_id, params.at(cur_par).fragment_length, 1.0, self.gc_fragment_content_bias_, self.fragment_surroundings_bias_);
		len.at( params.at(cur_par).fragment_length ) += sum;
		ref.at( params.at(cur_par).ref_seq_id ) += sum;

		if( params.size() > 20 && !(cur_par%(params.size()/20)) ){
			lock_guard<mutex> lock(result_mutex);
			printInfo << "Finished " << static_cast<uintPercentPrint>(Percent(cur_par, params.size())) << "% of the bias sums." << std::endl;
		}
	}

	result_mutex.lock();
	for(uintSeqLen l=len.size(); l--; ){
		insert_length_sum.at(l) += len.at(l);
	}
	for(uintRefSeqId r=ref.size(); r--; ){
		ref_seq_sum.at(r) += ref.at(r);
	}
	result_mutex.unlock();
}

void FragmentDistributionStats::BiasNormalizationThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<uintNumFits> &current_param,
		double &norm,
		mutex &result_mutex,
		double &max_bias){
	decltype(params.size()) cur_par(current_param++);

	double tmp_norm(0.0), tmp_max_bias(0.0);

	for(; cur_par < params.size(); cur_par = current_param++){
		tmp_norm += reference.SumBias(tmp_max_bias, params.at(cur_par).ref_seq_id, params.at(cur_par).fragment_length, self.ref_seq_bias_.at(params.at(cur_par).ref_seq_id)*self.insert_lengths_bias_.at(params.at(cur_par).fragment_length), self.gc_fragment_content_bias_, self.fragment_surroundings_bias_);
	}

	result_mutex.lock();
	if(tmp_max_bias > max_bias){
		max_bias = tmp_max_bias;
	}
	norm += tmp_norm;
	result_mutex.unlock();
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

void FragmentDistributionStats::Prepare( const Reference &ref, uintSeqLen maximum_insert_length, uintSeqLen max_ref_seq_bin_size, const vector<uintFragCount> &reads_per_frag_len_bin ){
	tmp_abundance_.resize(ref.NumberSequences());
	tmp_insert_lengths_.resize(maximum_insert_length+1);
	tmp_gc_fragment_content_.resize(101);

	for(auto strand=2; strand--; ){
		for(auto nuc=4; nuc--; ){
			tmp_outskirt_content_.at(strand).at(nuc).resize(2*kOutskirtRange);
		}
	}

	PrepareBiasCalculation( ref, maximum_insert_length, max_ref_seq_bin_size, reads_per_frag_len_bin );
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
		for(;  outskirt_pos-- && ref_pos < reference.SequenceLength(record_start.rNextId);){
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

	ref_seq_in_nxx_.clear();
	ref_seq_in_nxx_.shrink_to_fit();
	ref_seq_start_bin_.clear();
	ref_seq_start_bin_.shrink_to_fit();
	ref_seq_bin_def_.clear();
	ref_seq_bin_def_.shrink_to_fit();
}

void FragmentDistributionStats::CalculateInsertLengthAndRefSeqBias(const Reference &reference, uintNumThreads num_threads, vector<vector<VectorAtomic<uintFragCount>>> &site_count_by_insert_length_gc){
	// Sum up the other biases for the insert length, reference sequences pairs
	atomic<uintNumFits> current_param(0);
	vector<BiasCalculationParams> params;
	FillParams(params, reference);

	mutex result_mutex;

	std::vector<double> bias_sum_ref_seq, bias_sum_insert_length;
	bias_sum_ref_seq.resize(abundance_.size(), 0.0);
	bias_sum_insert_length.resize(insert_lengths_.to(), 0.0);

	site_count_by_insert_length_gc.resize(insert_lengths_.to());
	for(auto &site_count_by_gc : site_count_by_insert_length_gc){
		site_count_by_gc.resize(101);
	}

	// Run the predetermined parameters in defined number of threads
	thread threads[num_threads];
	for(auto i = num_threads; i--; ){
		threads[i] = thread(BiasSumThread, std::cref(*this), std::cref(reference), std::cref(params), std::ref(current_param), std::ref(bias_sum_insert_length), std::ref(bias_sum_ref_seq), std::ref(site_count_by_insert_length_gc), std::ref(result_mutex) );
	}
	for(auto i = num_threads; i--; ){
		threads[i].join();
	}

	// Divide the counts by the other biases to get ref_seq and insert_length bias
	// Bias is equal for both strands so it is only calculated for one, but the number of fragments at a site is calculated separately per strand: This means we have to half the calculated number of reads
	insert_lengths_bias_.Clear();
	for(auto i=insert_lengths_.to(); --i; ){
		if(0.0 != bias_sum_insert_length.at(i)){
			insert_lengths_bias_[i] = insert_lengths_.at(i)/bias_sum_insert_length.at(i)/2;
		}
		else{
			insert_lengths_bias_[i] = 0.0;
		}
	}
	ref_seq_bias_.clear();
	ref_seq_bias_.resize(abundance_.size(), 0.0);
	for(auto i=abundance_.size(); i--; ){
		if(0.0 != bias_sum_ref_seq.at(i)){
			ref_seq_bias_.at(i) = abundance_.at(i)/bias_sum_ref_seq.at(i)/2;
		}
	}

	// Find maximum (to normalize it to 1.0 instead of normalizing the sum to 1.0 as it is normal for probabilities)
	double max_len_bias(0.0);
	for(auto i=insert_lengths_bias_.to(); --i; ){
		SetToMax(max_len_bias, insert_lengths_bias_.at(i));
	}

	double max_seq_bias(0.0);
	for(auto i=ref_seq_bias_.size(); i--; ){
		SetToMax(max_seq_bias, ref_seq_bias_.at(i));
	}

	// Normalize so that maximum bias is 1.0
	for(auto i=insert_lengths_bias_.to(); --i; ){
		insert_lengths_bias_.at(i) /= max_len_bias;
	}

	for(auto i=ref_seq_bias_.size(); i--; ){
		ref_seq_bias_.at(i) /= max_seq_bias;
	}
}

bool FragmentDistributionStats::FinalizeBiasCalculation(const Reference &reference, uintNumThreads num_threads, FragmentDuplicationStats &duplications){
	vector<vector<VectorAtomic<uintFragCount>>> site_count_by_insert_length_gc;
	if(calculate_bias_){
		if( !StoreBias() ){
			return false;
		}

		CalculateInsertLengthAndRefSeqBias(reference, num_threads, site_count_by_insert_length_gc);
	}
	else{
		printInfo << "No bias calculation done. Storing uniform coverage." << std::endl;
		SetUniformBias();
	}

	// Free not needed memory at the end
	fragment_sites_by_ref_seq_bin_.clear();
	fragment_sites_by_ref_seq_bin_.shrink_to_fit();
	fragment_sites_by_ref_seq_bin_by_insert_length_.clear();
	fragment_sites_by_ref_seq_bin_by_insert_length_.shrink_to_fit();
	duplications.FinalizeDuplicationVector(site_count_by_insert_length_gc);

	if(calculate_bias_){
		printInfo << "Finished bias calculation" << std::endl;
	}
	else{
		printInfo << "Finished collecting duplicates" << std::endl;
	}

	return true;
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
									printErr <<  "Could not convert bias starting at postion " << bias_sep+1 << " in line " << nline << " in file " << bias_file << " into double\nException thrown: " << e.what() << std::endl;
								}
								IncreaseErrorCounter(errors);
								bias = 0.0;
							}

							if(0.0 > bias){
								if( errors < kMaxErrorsShownPerFile){
									printErr <<  "Bias in line " << nline << " in file " << bias_file << " is negative: " << bias << std::endl;
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
						printErr <<  "Could not find bias for reference sequence " << ref.ReferenceId(seq) << " in file " << bias_file << std::endl;
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

double FragmentDistributionStats::CalculateBiasNormalization(double &max_bias, const Reference &reference, uintNumThreads num_threads, uintFragCount total_reads) const{
	atomic<uintNumFits> current_param(0);
	vector<BiasCalculationParams> params;
	FillParamsSimulation(params);

	mutex result_mutex;
	double normalization(0.0);

	// Run the predetermined parameters in defined number of threads
	max_bias = 0.0;
	thread threads[num_threads];
	for(auto i = num_threads; i--; ){
		threads[i] = thread(BiasNormalizationThread, std::cref(*this), std::cref(reference), std::cref(params), std::ref(current_param), std::ref(normalization), std::ref(result_mutex), std::ref(max_bias) );
	}
	for(auto i = num_threads; i--; ){
		threads[i].join();
	}

	return total_reads / (normalization * 2); // Bias is equal for both strands so it is only calculated for one, but the number of fragments at a site is calculated separately per strand: This means we have to half the calculated number of reads
}

double FragmentDistributionStats::CalculateNonZeroThreshold(double bias_normalization, double max_bias) const{
	// Threshold that random number from fragment generation has to reach or fragment count will be zero and does not have to be calculated
	double max_dispersion( Dispersion(bias_normalization*max_bias) );

	return pow(max_dispersion/(max_dispersion+bias_normalization*max_bias), max_dispersion); // Lowest possible F(0) for negative binomial
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

reseq::uintDupCount FragmentDistributionStats::GetFragmentCounts(double bias_normalization, uintRefSeqId ref_seq_id, uintSeqLen fragment_length, uintPercent gc, const Surrounding &fragment_start, const Surrounding &fragment_end, double probability_chosen, double non_zero_threshold) const{
	if(probability_chosen < non_zero_threshold){
		return 0;
	}
	else{
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
