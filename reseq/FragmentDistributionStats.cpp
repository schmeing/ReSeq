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
#include <math.h>
using std::abs;
using std::exp;
using std::isnan;
using std::log;
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
using seqan::FunctorComplement;

//include "utilities.h"
using reseq::utilities::Divide;
using reseq::utilities::IntPow;
using reseq::utilities::Percent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;
using reseq::utilities::Sign;
using reseq::utilities::InvLogit2;
using reseq::utilities::VectorAtomic;

constexpr double BiasCalculationVectors::lower_bound_;
constexpr double BiasCalculationVectors::upper_bound_;
constexpr double BiasCalculationVectors::base_value_;

void BiasCalculationVectors::AddCountsFromSite(const FragmentSite &site, array<uint64_t, 101> &gc_count, array<uint64_t, 4*Reference::num_surrounding_blocks_*Reference::surrounding_range_> &sur_count){
	if(0 < site.count_forward_ + site.count_reverse_){
		gc_count.at(site.gc_) += site.count_forward_ + site.count_reverse_;

		for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
			auto sur_start = site.start_surrounding_.at(block);
			auto sur_end = site.end_surrounding_.at(block);

			for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
				auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);

				sur_count.at(base_pos + sur_start%4) += site.count_forward_ + site.count_reverse_;
				sur_count.at(base_pos + sur_end%4) += site.count_forward_ + site.count_reverse_;

				sur_start /= 4;
				sur_end /= 4;
			}
		}
	}
}

void BiasCalculationVectors::GetCounts(){
	gc_count_.fill(0);
	sur_count_.fill(0);
	duplication_count_.fill(0);
	for(auto &site : sites_){
		if(0.0 != site.bias_){
			AddCountsFromSite(site, gc_count_, sur_count_);
		}

		++duplication_count_.at(min(site.count_forward_, max_duplications_+1));
		++duplication_count_.at(min(site.count_reverse_, max_duplications_+1));
	}

	total_counts_ = SumVect(gc_count_);
}

void BiasCalculationVectors::RemoveUnnecessarySites(){
	gc_sites_.fill(0);
	sur_sites_.fill(0);

	uint32_t needed_sites(0);
	bool zero_site;
	for(uint32_t cur_site = 0; cur_site < sites_.size(); ++cur_site){
		if(0.0 != sites_.at(cur_site).bias_){ // Filter out sides where at least one surrounding couldn't be determined
			zero_site = false; // sites with zero for gc are ok as we use splines there

			for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
				auto sur_start = sites_.at(cur_site).start_surrounding_.at(block);
				auto sur_end = sites_.at(cur_site).end_surrounding_.at(block);

				for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
					auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);
					if( !sur_count_.at(base_pos + sur_start%4) ){
						zero_site = true;
					}
					if( !sur_count_.at(base_pos + sur_end%4) && sur_start%4 != sur_end%4 ){
						zero_site = true;
					}

					sur_start /= 4;
					sur_end /= 4;
				}
			}

			if(!zero_site){
				sites_.at(needed_sites++) = sites_.at(cur_site);

				gc_sites_.at(sites_.at(cur_site).gc_) += 2;

				for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
					auto sur_start = sites_.at(cur_site).start_surrounding_.at(block);
					auto sur_end = sites_.at(cur_site).end_surrounding_.at(block);

					for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
						auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);

						sur_sites_.at(base_pos + sur_start%4) += 2;
						sur_sites_.at(base_pos + sur_end%4) += 2;

						sur_start /= 4;
						sur_end /= 4;
					}
				}
			}
		}
	}
	sites_.resize(needed_sites);

	total_sites_ = 2*needed_sites; // Forward + Reverse
}

void BiasCalculationVectors::AddBaseValue(uint32_t count){
	if(1 < count){
		double k_fac = 2;
		for(uint16_t k = 3; k <= count; ++k){
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

void BiasCalculationVectors::DeactivateZeroCounts(vector<double> &par, double deactive_value, uint16_t sur_shift){
	for( uint16_t sur=sur_count_.size(); sur--; ){
		if(!sur_count_.at(sur)){
			par.at(sur+sur_shift) = deactive_value;
		}
	}
}

void BiasCalculationVectors::NormGC(){
	// Sort gc by number of sites
	array<pair<uint64_t, uint16_t>, 101> sites;
	for(uint16_t gc = gc_sites_.size(); gc--; ){
		sites.at(gc) = {gc_sites_.at(gc), gc};
	}
	sort(sites.begin(), sites.end());

	// Sum the biases up to half the sites starting with the gc with most sites
	uint64_t sum_sites = 0;
	double sum_bias = 0.0;
	uint16_t num_bias = 0;
	for(uint16_t gci = gc_sites_.size(); gci-- && sum_sites < Divide(total_sites_*percent_gc_sites_for_normalization_,100); ){
		sum_sites += sites.at(gci).first;
		sum_bias += gc_bias_.at( sites.at(gci).second );
		++num_bias;
	}

	// Normalize using the biases of the gc values having 50% of the sites to take only well determined values
	for(uint16_t gc = gc_sites_.size(); gc--; ){
		gc_bias_.at(gc) *= num_bias/sum_bias;
	}
}

void BiasCalculationVectors::NormSurroundings(const vector<double> &x){
	if(sur_sum_){
		for(auto sur_pos = 0; sur_pos < Reference::num_surrounding_blocks_*Reference::surrounding_range_; ++sur_pos){
			auto from = 1+sur_pos*4;
			auto to = 1+(sur_pos+1)*4;

			double sur_sum(0.0);
			uint16_t valid_sur = 0;
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
	else if(sur_mult_){
		for(auto sur_pos = 0; sur_pos < Reference::num_surrounding_blocks_*Reference::surrounding_range_; ++sur_pos){
			double sur_sum = 0.0;
			uint16_t valid_sur = 0;
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
}

void BiasCalculationVectors::UnnormSurroundingGradients(vector<double> &grad, const vector<double> &x){
	if(sur_mult_){
		for(auto sur_pos = 0; sur_pos < Reference::num_surrounding_blocks_*Reference::surrounding_range_; ++sur_pos){
			auto from = sur_pos*4;
			auto to = (sur_pos+1)*4;

			double grad_sum(0.0), bias_sum(0.0);
			uint16_t valid_sur = 0;
			for(auto sur = from; sur < to; ++sur){
				if(sur_count_.at(sur)){
					grad_sum += sur_bias_.at(sur) * grad.at(sur);
					bias_sum += x.at(sur);
					++valid_sur;
				}
			}

			for(auto sur = from; sur < to; ++sur){
				grad.at(sur) = (valid_sur*grad.at(sur)-grad_sum)/bias_sum;
			}
		}
	}
	else if(sur_sum_){
		for(auto sur_pos = 0; sur_pos < Reference::num_surrounding_blocks_*Reference::surrounding_range_; ++sur_pos){
			auto from = sur_pos*4;
			auto to = (sur_pos+1)*4;

			double grad_mean(0.0);
			uint16_t valid_sur = 0;
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
		for(uint16_t sur = 0; sur < 4; ++sur){
			grad.at(0) += sur_grad_.at(sur);
		}
	}
}

double BiasCalculationVectors::SurBiasAtSite(const pair<double, double>& bias){
	if(sur_sum_){
		return static_cast<double>(total_counts_) / total_sites_ * InvLogit2(bias.first) * InvLogit2(bias.second) + base_value_;
	}
	else if(sur_mult_){
		return static_cast<double>(total_counts_) / total_sites_ * bias.first * bias.second;
	}
}

pair<double, double> BiasCalculationVectors::SurBiasAtSiteSplit(const FragmentSite& site){
	pair<double, double> bias;
	if(sur_mult_){
		bias = {1.0, 1.0};
	}
	else if(sur_sum_){
		bias = {0.0, 0.0};
	}

	for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
		auto sur_start = site.start_surrounding_.at(block);
		auto sur_end = site.end_surrounding_.at(block);

		for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
			auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);

			if(sur_sum_){
				bias.first += sur_bias_.at(base_pos + sur_start%4);
				bias.second += sur_bias_.at(base_pos + sur_end%4);
			}
			else if(sur_mult_){
				bias.first *= sur_bias_.at(base_pos + sur_start%4);
				bias.second *= sur_bias_.at(base_pos + sur_end%4);
			}

			sur_start /= 4;
			sur_end /= 4;
		}
	}

	return bias;
}

double BiasCalculationVectors::SurBiasAtSite(const FragmentSite& site){
	auto bias = SurBiasAtSiteSplit(site);
	return SurBiasAtSite(bias);
}

void BiasCalculationVectors::DefineStartingKnots(){
	// Define knots for gc
	uint64_t sum(0);
	uint16_t k(0);
	for(uint16_t gc=0; gc < gc_sites_.size(); ++gc){
		if( gc_sites_.at(gc) && gc_count_.at(gc) ){
			// Get an equal spacing of sites (First site where the cumulated sum is at least the k-th quartile)
			sum += gc_sites_.at(gc);
			gc_knots_.at(k) = gc;
			if(sum >= Divide(k*total_sites_,gc_spline_df_-1)){
				++k;
			}
		}
	}

	// In case there are not enough gc bins with counts
	if( k < gc_knots_.size() ){
		// Fill rest of array with invalid gc, to not jump over an unfilled bin that accidentally happens to have the right gc
		for(auto k_fill=k; k_fill<gc_knots_.size(); ++k_fill){
			gc_knots_.at(k_fill) = 102;
		}

		// Start from the first knot and assign unused gc bins with zero counts
		uint16_t next_k = 1;
		uint16_t gc = gc_knots_.at(0) + 1;
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
	const uint16_t df = gc_spline_df_;
	auto &x = gc_knots_;

	array<double, df> l; // j=[0, n]
	array<array<double, df>, df> c, z; // c[j][a_index]: j=[0, n]
	array<double, df-1> mu, h; // i=[0, n-1]
	auto &b = lin_comb_gc_splines_.at(0); // b[i][a_index]: i=[0, n-1]
	auto &d = lin_comb_gc_splines_.at(2); // d[i][a_index]: i=[0, n-1]
	array<array<double, df>, df-1> beta; //  beta[k][a_index]: k=[1, n-1]

	for(uint16_t i=0; i < h.size(); ++i){
		h.at(i) = x.at(i+1) - x.at(i);
	}

	for(uint16_t k=1; k < beta.size(); ++k){
		beta.at(k).fill(0.0);
		beta.at(k).at(k+1) = 3/h.at(k);
		beta.at(k).at(k) = -3/h.at(k) - 3/h.at(k-1);
		beta.at(k).at(k-1) = 3/h.at(k-1);
	}

	l.at(0) = 0.0;
	mu.at(0) = 0.0;
	z.at(0).fill(0.0);

	for(uint16_t k=1; k < beta.size(); ++k){
		l.at(k) = 2*(x.at(k+1)-x.at(k-1)) - h.at(k-1)*mu.at(k-1);
		mu.at(k) = h.at(k)/l.at(k);
		for(uint16_t ai=0; ai<df; ++ai){
			z.at(k).at(ai) = (beta.at(k).at(ai) - h.at(k-1) * z.at(k-1).at(ai))/l.at(k);
		}
	}

	l.at(l.size()-1) = 1.0;
	z.at(z.size()-1).fill(0.0);
	c.at(c.size()-1).fill(0.0);

	for(uint16_t i=h.size(); i--;){
		for(uint16_t ai=0; ai<df; ++ai){
			c.at(i).at(ai) = z.at(i).at(ai) - mu.at(i)*c.at(i+1).at(ai);
			b.at(i).at(ai) = -h.at(i)*(c.at(i+1).at(ai) + 2*c.at(i).at(ai))/3;
			d.at(i).at(ai) = (c.at(i+1).at(ai) - c.at(i).at(ai))/3/h.at(i);
		}

		b.at(i).at(i+1) += 1/h.at(i);
		b.at(i).at(i) -= 1/h.at(i);
	}

	// Store values from c in lin_comb_gc_splines_ (as c needs to be one larger for the calculation it is not a reference like b and d)
	for(uint16_t i=h.size(); i--;){
		for(uint16_t ai=0; ai<df; ++ai){
			lin_comb_gc_splines_.at(1).at(i).at(ai) = c.at(i).at(ai);
		}
	}
}

void BiasCalculationVectors::GetSplineCoefficients(double &a, double &b, double &c, double &d, uint16_t k, const vector<double> &spline_pars){
	a = spline_pars.at(k+1);
	b = 0.0;
	c = 0.0;
	d = 0.0;

	for(uint16_t ai=1; ai < gc_spline_pars_.size(); ++ai){
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
	double lower_const = norm*InvLogit2(lower_const_no_logit) + base_value_;
	for(uint16_t gc=0; gc < gc_knots_.at(0); ++gc){
		gc_bias_no_logit_.at(gc) = lower_const_no_logit;
		gc_bias_.at(gc) = lower_const;
	}

	// Fill observed range of gc with spline values
	gc_bias_no_logit_.at( gc_knots_.at(0) ) = a;
	gc_bias_.at( gc_knots_.at(0) ) = norm*InvLogit2(a) + base_value_;
	uint16_t k=0;
	for(uint16_t gc=gc_knots_.at(0)+1; gc < gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		if(gc == gc_knots_.at(k+1)){
			// Reached new spline part
			GetSplineCoefficients(a, b, c, d, ++k, spline_pars);
			gc_bias_no_logit_.at(gc) = a;
			gc_bias_.at(gc) = norm*InvLogit2(a) + base_value_;
		}
		else{
			// Continue on current spline part
			uint32_t cur_gc = gc-gc_knots_.at(k);
			gc_bias_no_logit_.at(gc) = a + b*cur_gc + c*cur_gc*cur_gc + d*cur_gc*cur_gc*cur_gc;
			gc_bias_.at(gc) = norm*InvLogit2( gc_bias_no_logit_.at(gc) ) + base_value_;
		}
	}
	++k;

	gc_bias_no_logit_.at( gc_knots_.at(k) ) = spline_pars.at(k);
	gc_bias_.at( gc_knots_.at(k) ) = norm*InvLogit2(spline_pars.at(k)) + base_value_;

	// Constant bias values for gc above any observed sites: Extrapolate the linear function at the highest measured GC by half a GC unit and use this value
	uint32_t cur_gc = gc_knots_.at(k)-gc_knots_.at(k-1);
	double upper_const_no_logit = gc_spline_pars_.at(k) + 0.5*(b+c*cur_gc);
	double upper_const = norm*InvLogit2( upper_const_no_logit ) + base_value_;
	for(uint16_t gc=gc_knots_.at(k)+1; gc < gc_bias_.size(); ++gc){
		gc_bias_no_logit_.at(gc) = upper_const_no_logit;
		gc_bias_.at(gc) = upper_const;
	}
}

void BiasCalculationVectors::CalculateSplineGrad(const vector<double> &spline_pars, vector<double> &grad, uint16_t grad_offset){
	for( uint16_t gc=0; gc < gc_bias_grad_.size(); ++gc ){
		grad.at(grad_offset) += gc_bias_grad_.at(gc) / spline_pars.at(0);
		gc_bias_grad_.at(gc) *= InvLogit2(-gc_bias_no_logit_.at(gc))/2;
	}

	// Constant bias values for gc below any observed sites
	for( uint16_t gc=0; gc < gc_knots_.at(0); ++gc ){
		grad.at(grad_offset+1) += gc_bias_grad_.at(gc);

		for(uint16_t ai=0; ai < gc_spline_pars_.size()-1; ++ai){
			grad.at(grad_offset+ai+1) -= 0.5*lin_comb_gc_splines_.at(0).at(0).at(ai) * gc_bias_grad_.at(gc);
		}
	}

	// Observed range of gc
	int16_t k=-1;
	for( uint16_t gc=gc_knots_.at(0); gc <= gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		if(gc == gc_knots_.at(k+1)){
			// Reached new spline part
			++k;
		}
		else{
			// Continue on current spline part
			uint32_t cur_gc = gc-gc_knots_.at(k);
			for(uint16_t ai=0; ai < gc_spline_pars_.size()-1; ++ai){
				grad.at(grad_offset+ai+1) += (lin_comb_gc_splines_.at(0).at(k).at(ai)*cur_gc + lin_comb_gc_splines_.at(1).at(k).at(ai)*cur_gc*cur_gc + lin_comb_gc_splines_.at(2).at(k).at(ai)*cur_gc*cur_gc*cur_gc) * gc_bias_grad_.at(gc);
			}
		}
		grad.at(grad_offset+k+1) += gc_bias_grad_.at(gc);
	}

	// Constant bias values for gc above any observed sites
	uint32_t cur_gc = gc_knots_.at(k)-gc_knots_.at(k-1);
	for( uint16_t gc=gc_knots_.at( gc_knots_.size()-1 )+1; gc < gc_count_.size(); ++gc ){
		grad.at(grad_offset+k+1) += gc_bias_grad_.at(gc);
		for(uint16_t ai=0; ai < gc_spline_pars_.size()-1; ++ai){
			grad.at(grad_offset+ai+1) += 0.5*(lin_comb_gc_splines_.at(0).at(k-1).at(ai) + lin_comb_gc_splines_.at(1).at(k-1).at(ai)*cur_gc) * gc_bias_grad_.at(gc);
		}
	}
}

double BiasCalculationVectors::LogLikeGcSpline(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_spline_;

	calc.CalculateSpline(x);

	double loglike(calc.loglike_pois_base_);
	for( uint16_t gc=0; gc < calc.gc_count_.size(); ++gc){
		loglike += calc.gc_bias_sum_.at(gc).first + calc.gc_count_.at(gc)*log(calc.gc_bias_.at(gc)) - calc.gc_bias_.at(gc)*calc.gc_bias_sum_.at(gc).second;
	}

	if( grad.size() ){
		for( auto &g : grad ){
			g = 0.0;
		}

		for( uint16_t gc=0; gc < calc.gc_count_.size(); ++gc){
			//double grad_gc = gc_count_.at(gc)/gc_bias_.at(gc) - gc_bias_sum_.at(gc).second;
			calc.gc_bias_grad_.at(gc) = (calc.gc_count_.at(gc) - calc.gc_bias_sum_.at(gc).second*calc.gc_bias_.at(gc));
		}

		calc.CalculateSplineGrad(x, grad, 0);

		for(uint16_t gr=grad.size(); gr--; ){
			grad.at(gr) *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;
		}
	}

	loglike *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;

	return loglike;
}

double BiasCalculationVectors::OptimizeSpline(){
	PrepareSplines();

	double max_val = 0.0;
	for( uint16_t k=gc_spline_pars_.size()-1; k--; ){
		SetToMax(max_val, gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second);
	}
	max_val /= 1.6;
	gc_spline_pars_.at(0) = max_val;

	max_val *= 2;
	for( uint16_t k=gc_spline_pars_.size()-1; k--; ){
		if( gc_count_.at( gc_knots_.at(k) ) ){
			double val = gc_count_.at( gc_knots_.at(k) )/gc_bias_sum_.at( gc_knots_.at(k) ).second;

			gc_spline_pars_.at(k+1) = log(val/(max_val-val)); // Set gc bias to set the gradient to zero for this gc

			if(upper_bound_ < gc_spline_pars_.at(k+1)){
				gc_spline_pars_.at(k+1) = upper_bound_;
			}
			else if(lower_bound_ > gc_spline_pars_.at(k+1)){
				gc_spline_pars_.at(k+1) = lower_bound_;
			}
		}
		else{
			gc_spline_pars_.at(k+1) = lower_bound_;
		}
	}

	func_calls_spline_ = 0;
	std::vector<double> grad;
	double loglike;

	converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == optimizer_spline_.optimize(gc_spline_pars_, loglike) ){
			func_calls_spline_ += 10000;
			converged_ = false;
		}
	}
	catch(const std::exception& e){
		func_calls_spline_ += 20000;
		converged_ = false;
	}

	return loglike;
}

double BiasCalculationVectors::OptimizeSpline(uint16_t knot, int16_t shift){
	// knot is assumed to be an inner knot
	if( shift < 0 && gc_knots_.at(knot-1) >= gc_knots_.at(knot)+shift || shift > 0 && gc_knots_.at(knot+1) <= gc_knots_.at(knot)+shift || 0 == gc_count_.at( gc_knots_.at(knot) + shift ) ){
		return -1 * numeric_limits<double>::infinity();
	}
	else{
		gc_knots_.at(knot) += shift;
		double loglike = OptimizeSpline();
		gc_knots_.at(knot) -= shift;
		return loglike;
	}
}

uint16_t BiasCalculationVectors::OptimizeSplineIdToKnot(int16_t id, const uint16_t range){
	return (id+2*range-1)/(2*range);
}

int16_t BiasCalculationVectors::OptimizeSplineIdToShift(int16_t id, const uint16_t range){
	return (id-1)%(2*range) - range + ((id-1)%(2*range) < range?0:1);
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

	for( uint16_t gc=gc_knots_.at(0); gc <= gc_knots_.at( gc_knots_.size()-1 ); ++gc){
		gc_bias_sum_.at(gc).second *= 2; // Multiplication with 2 here to save one multiplication per site
	}

	if(parameter_info_file || dispersion_info_file){
		// Get "perfect" GC bias by setting gradient to zero for each gc bin
		for( uint16_t gc=0; gc < gc_count_.size(); ++gc){
			if(gc_count_.at(gc)){
				gc_bias_.at(gc) = gc_count_.at(gc)/gc_bias_sum_.at(gc).second;
			}
			else{
				gc_bias_.at(gc) = 0;
			}
		}

		if(dispersion_info_file){
			WriteOutInformation("poisson");
		}

		NormGC();
		for(uint16_t gc=gc_bias_.size(); gc--; ){
			gc_bias_pois_.at(gc) = gc_bias_.at(gc);
		}
	}

	// Optimize gc spline
	// Knots are adjusted by comparing the likelihoods to the ones with all inner knots in turn moved by 1 to range up or down and taking the best likelihood in a greedy fashion
	// Iterate until likelihood cannot be improved by moving an inner knot by 1 to range gc
	const uint16_t range(max_knot_shift_);
	array<double, (gc_spline_df_-2)*2*range+1 > loglikes;
	uint16_t bestlike(1);
	loglikes.at(0) = OptimizeSpline();
	loglikes.at(1) = OptimizeSpline(1, -static_cast<int16_t>(range));
	while(bestlike != 0){ // loglikes[0] is current knot position
		// Calculate loglikelihoods with inner knots in turn moved by one up or down
		for(uint16_t cur_like=loglikes.size(); --cur_like; ){ // Leave out 0
			if(cur_like != bestlike){ // Otherwise we still have the value from last iteration
				loglikes.at(cur_like) = OptimizeSpline( OptimizeSplineIdToKnot(cur_like, range), OptimizeSplineIdToShift(cur_like, range) );
			}
		}

		bestlike = distance( loglikes.begin(), max_element(loglikes.begin(), loglikes.end()) );
		if(bestlike != 0){
			// Save the two likelihoods that will be identical in the next round in the proper spots
			uint16_t flipped_bestlike = bestlike - 2*(static_cast<int16_t>((bestlike-1)%(2*range)) - range) - 1;
			loglikes.at(flipped_bestlike) = loglikes.at(0);
			loglikes.at(0) = loglikes.at(bestlike);
			// Apply the knot movement
			gc_knots_.at( OptimizeSplineIdToKnot(bestlike, range) ) += OptimizeSplineIdToShift(bestlike, range);
			// Store at which new position the original likelihood is, so that we don't recalculate it
			bestlike = flipped_bestlike;
		}
	}

	// Get gc bias for best spline (We require a higher precision here as the surrounding gradient assumes that the gc gradients are exactly zero, so we don't want to deviate much)
	loglike_spline_ = OptimizeSpline();

	if(parameter_info_file || dispersion_info_file){
		CalculateSpline(gc_spline_pars_);

		if(dispersion_info_file){
			WriteOutInformation("gcspline");
		}

		NormGC();
		for(uint16_t gc=gc_bias_.size(); gc--; ){
			gc_bias_spline_.at(gc) = gc_bias_.at(gc);
		}
	}
}

double BiasCalculationVectors::LogLikelihoodPoisson(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_pois_;

	calc.NormSurroundings(x);

	if( grad.size() ){
		if(BiasCalculationVectors::sur_mult_){
			for(uint16_t sur=calc.sur_count_.size(); sur--;){
				grad.at(sur) = calc.sur_count_.at(sur)/calc.sur_bias_.at(sur);
			}
		}
		else if(BiasCalculationVectors::sur_sum_){
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
			auto cur_grad_start = InvLogit2(-bias_split.first)/2;
			auto cur_grad_end = InvLogit2(-bias_split.second)/2;

			uint64_t sur_start, sur_end;
			for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
				sur_start = site.start_surrounding_.at(block);
				sur_end = site.end_surrounding_.at(block);

				for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
					auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);
					auto i_start = base_pos + sur_start%4;
					auto i_end = base_pos + sur_end%4;

					if(BiasCalculationVectors::sur_mult_){
						calc.grad_gc_bias_sum_.at(site.gc_).at(i_start).first += bias/calc.sur_bias_.at(i_start);
						calc.grad_gc_bias_sum_.at(site.gc_).at(i_end).first += bias/calc.sur_bias_.at(i_end);
					}
					else if(BiasCalculationVectors::sur_sum_){
						calc.grad_gc_bias_sum_.at(site.gc_).at(i_start).first += bias*cur_grad_start;
						calc.grad_gc_bias_sum_.at(site.gc_).at(i_end).first += bias*cur_grad_end;
					}

					if(site.count_forward_+site.count_reverse_){
						if(BiasCalculationVectors::sur_sum_){
							calc.sur_grad_.at(i_start) += (site.count_forward_+site.count_reverse_)*cur_grad_start;
							calc.sur_grad_.at(i_end) += (site.count_forward_+site.count_reverse_)*cur_grad_end;
						}
					}

					sur_start /= 4;
					sur_end /= 4;
				}
			}
		}
	}

	for( uint16_t gc=calc.gc_count_.size(); gc--; ){
		calc.gc_bias_sum_.at(gc).second *= 2; // Multiplication with 2 here to save one multiplication per site
	}

	// Calculate gc bias by setting gradient to zero for every gc bin
	double loglike = calc.loglike_pois_base_;
	for( uint16_t gc=calc.gc_count_.size(); gc--; ){
		if(calc.gc_count_.at(gc)){
			calc.gc_bias_.at(gc) = calc.gc_count_.at(gc)/calc.gc_bias_sum_.at(gc).second;
			loglike += calc.gc_bias_sum_.at(gc).first + calc.gc_count_.at(gc)*log(calc.gc_bias_.at(gc)) - calc.gc_bias_.at(gc)*calc.gc_bias_sum_.at(gc).second;
		}
		else{
			calc.gc_bias_.at(gc) = 0;
		}
	}

	// Calculate gradients for surroundings based on the gc bias
	if( grad.size() ){
		for( uint16_t gc=calc.gc_count_.size(); gc--; ){
			for(auto sur=calc.sur_bias_.size(); sur--;){
				calc.grad_gc_bias_sum_.at(gc).at(sur).first *= 2; // Multiplication with 2 here to save two multiplication per site
				if(BiasCalculationVectors::sur_mult_){
					grad.at(sur) -= calc.gc_bias_.at(gc) * calc.grad_gc_bias_sum_.at(gc).at(sur).first;
				}
				else if(BiasCalculationVectors::sur_sum_){
					calc.sur_grad_.at(sur) -= calc.gc_bias_.at(gc) * calc.grad_gc_bias_sum_.at(gc).at(sur).first;
				}
			}
		}

		calc.UnnormSurroundingGradients(grad, x);

		for(uint16_t gr=grad.size(); gr--; ){
			grad.at(gr) *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;
		}
	}

	loglike *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;

	return loglike;
}

double BiasCalculationVectors::GetDispersion(double bias, double a, double b){
	double r = bias/(a+b*bias);
	if(r > bias*1e14){
		// r is so much larger than the bias that r + bias is having precision issues, so shift it back into precision (it anyways is approximately poisson at this point)
		r = bias*1e14;
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

	for(uint16_t i=1; i <= site.count_forward_; ++i){
		grad_term2 += 1/(r+(i-1)); // Be careful to handle r->0 (brackets around i-1, so that r+1-1=r and not 0)
		loglike += log( (r+(i-1))/i );
	}
	for(uint16_t i=1; i <= site.count_reverse_; ++i){
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
		auto cur_grad_start = InvLogit2(-bias_split.first)/2;
		auto cur_grad_end = InvLogit2(-bias_split.second)/2;

		uint64_t sur_start, sur_end;
		for(uint16_t block = Reference::num_surrounding_blocks_; block--;){
			sur_start = site.start_surrounding_.at(block);
			sur_end = site.end_surrounding_.at(block);

			for(uint16_t sur_pos = Reference::surrounding_range_; sur_pos--; ){
				auto base_pos = 4*(block*Reference::surrounding_range_ + sur_pos);
				auto i_start = base_pos + sur_start%4;
				auto i_end = base_pos + sur_end%4;

				if(sur_mult_){
					grad.at(i_start) += grad_term1; // Division by bias later
					grad.at(i_end) += grad_term1;
				}
				else if(sur_sum_){
					sur_grad_.at(i_start) += grad_term1*cur_grad_start;
					sur_grad_.at(i_end) += grad_term1*cur_grad_end;
				}

				sur_start /= 4;
				sur_end /= 4;
			}
		}
	}
}

double BiasCalculationVectors::LogLikelihoodNbinom(const vector<double> &x, vector<double> &grad, void* f_data){
	BiasCalculationVectors &calc(*reinterpret_cast<BiasCalculationVectors *>(f_data));

	++calc.func_calls_nbinom_;

	calc.NormSurroundings(x);

	uint16_t spline_par_offset(calc.sur_count_.size());
	if(BiasCalculationVectors::sur_sum_){
		spline_par_offset += 1;
	}

	for(uint16_t k=0; k<calc.gc_spline_pars_.size(); ++k){
		calc.gc_spline_pars_.at(k) = x.at( spline_par_offset + k );
	}
	calc.CalculateSpline(calc.gc_spline_pars_);

	double loglike(0.0), a( x.at(x.size()-2) ), b( x.at(x.size()-1) ), ga(0.0), gb(0.0);
	if( grad.size() ){
		for(uint16_t sur=spline_par_offset; sur--;){
			grad.at(sur) = 0.0;
		}
		for(uint16_t k=spline_par_offset; k<grad.size(); ++k){
			grad.at(k) = 0.0;
		}

		if(BiasCalculationVectors::sur_sum_){
			calc.sur_grad_.fill(0.0);
		}
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
		if(BiasCalculationVectors::sur_mult_){
			for(auto sur=calc.sur_count_.size(); sur--;){
				grad.at(sur) /= calc.sur_bias_.at(sur); // Division by bias here instead of in the site loop
			}
		}
		calc.UnnormSurroundingGradients(grad, x);

		grad.at( grad.size()-2 ) = ga;
		grad.at( grad.size()-1 ) = gb;

		for(uint16_t gr=grad.size(); gr--; ){
			grad.at(gr) *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;
		}
	}

	loglike *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;

	return loglike;
}

void BiasCalculationVectors::PostProcessFit(){
	// Normalize sur biases to 4
	NormSurroundings(fit_pars_);

	// Get correct gc bias values
	for(uint16_t k=0; k<gc_spline_pars_.size(); ++k){
		gc_spline_pars_.at(k) = fit_pars_.at( sur_count_.size() + (sur_sum_?1:0) + k );
	}
	CalculateSpline(gc_spline_pars_);

	if(dispersion_info_file){
		WriteOutInformation("nbinom");
	}

	NormGC();

	dispersion_.resize(2);
	dispersion_.at(0) = fit_pars_.at( fit_pars_.size()-2 );
	dispersion_.at(1) = fit_pars_.at( fit_pars_.size()-1 );
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
	for(uint16_t i=1; i < calc.duplication_count_part_.size(); ++i){
		gr += calc.duplication_count_part_.at(i)/(r+i-1);
		loglike += calc.duplication_count_part_.at(i) * log( (r+i-1)/i );
	}

	if( grad.size() ){
		gr *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;
		grad.at(0) = gr;
	}

	loglike *= BiasCalculationVectors::loglike_mult_/calc.total_sites_;

	return loglike;
}

void BiasCalculationVectors::WriteOutInformation(const char *context){
	// Set predicted site bias
	for(auto &site : sites_){
		site.bias_ = SurBiasAtSite(site) * gc_bias_.at(site.gc_);
	}

	// Sort sites from lowest to highest bias
	sort(sites_.begin(), sites_.end(), FragmentSite::LowerBias);

	// Prepare optimizer
	dispersion_.resize(1);
	nlopt::opt optimizer(nlopt::LD_LBFGS, dispersion_.size());
	optimizer.set_xtol_rel(precision_aim_);
	optimizer.set_ftol_abs(stop_criterion_);
	optimizer.set_maxeval(max_likelihood_calculations_);

	vector<double> bounds;
	bounds.resize(1);
	bounds.at(0) = 1e-200;
	optimizer.set_lower_bounds(bounds);
	bounds.at(0) = 1e200;
	optimizer.set_upper_bounds(bounds);

	optimizer.set_max_objective(BiasCalculationVectors::LogLikelihoodConstDispersion, reinterpret_cast<void*>(this));

	// Prepare fits
	double loglike, sample_mean, prediction_mean;
	uint32_t max_duplication;
	vector<double> grad;
	ofstream myfile;
	myfile.open(dispersion_info_file, ofstream::out | ofstream::app);

	// Constant number of sites per bin
	const uint32_t num_sites = 200000;
	start_site_dispersion_fit_ = ((total_sites_)%num_sites)/2; // Ignore equally many sites in the beginning and end
	end_site_dispersion_fit_ = start_site_dispersion_fit_ + num_sites/2;

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
			++duplication_count_part_.at(min(sites_.at(s).count_forward_, max_duplications_+1));
			++duplication_count_part_.at(min(sites_.at(s).count_reverse_, max_duplications_+1));
		}
		sample_mean /= 2*(end_site_dispersion_fit_ - start_site_dispersion_fit_);
		prediction_mean /= (end_site_dispersion_fit_ - start_site_dispersion_fit_);

		myfile << context << ", " << prediction_mean << ", " << sample_mean << ", " << 2*(end_site_dispersion_fit_ - start_site_dispersion_fit_) << ", " << max_duplication;
		for(uint16_t dups=1; dups<4; ++ dups){
			myfile << ", " << duplication_count_part_.at(dups);

		}

		// Predicted means
		use_sample_mean_ = nan("0");
		dispersion_.at(0) = 1.0;

		func_calls_const_disp_ = 0;
		try{
			if( nlopt::MAXEVAL_REACHED == optimizer.optimize(dispersion_, loglike) ){
				//printWarn << "Optimizer reached maximum number of chi^2 calculations for seq " << ref_seq_id << " and length " << insert_length << std::endl;
				func_calls_const_disp_ += 10000;
			}
		}
		catch(const std::exception& e){
			//printWarn << "Optimizer crashed for seq " << ref_seq_id << " and length " << insert_length << " with: " << e.what() << std::endl;
			func_calls_const_disp_ += 20000;
		}
		myfile << ", " << func_calls_const_disp_ << ", " << dispersion_.at(0);

		// Sample mean
		use_sample_mean_ = sample_mean;
		dispersion_.at(0) = 1.0;

		func_calls_const_disp_ = 0;
		try{
			if( nlopt::MAXEVAL_REACHED == optimizer.optimize(dispersion_, loglike) ){
				//printWarn << "Optimizer reached maximum number of chi^2 calculations for seq " << ref_seq_id << " and length " << insert_length << std::endl;
				func_calls_const_disp_ += 10000;
			}
		}
		catch(const std::exception& e){
			//printWarn << "Optimizer crashed for seq " << ref_seq_id << " and length " << insert_length << " with: " << e.what() << std::endl;
			func_calls_const_disp_ += 20000;
		}

		myfile << ", " << func_calls_const_disp_ << ", " << dispersion_.at(0) << std::endl;

		// Prepare next iteration
		start_site_dispersion_fit_ += num_sites/2;
		end_site_dispersion_fit_ += num_sites/2;
	}

	myfile.close();
}

void FragmentDistributionStats::PrepareBiasCalculation( const Reference &ref, uint32_t maximum_insert_length, const vector<uint64_t> &reads_per_ref_seq_bin ){
	// Read in vectors to calculate bias from
	auto num_bins = ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1);
	fragment_sites_by_ref_seq_bin_.resize(num_bins);
	for( uint32_t ref_bin=0; ref_bin < fragment_sites_by_ref_seq_bin_.size(); ++ref_bin){
		fragment_sites_by_ref_seq_bin_.at(ref_bin).resize(reads_per_ref_seq_bin.at(ref_bin)/2); // Sites are defined by a read pair, so /2
	}
	fragment_sites_by_ref_seq_bin_cur_id_.resize(num_bins);
	fragment_sites_by_ref_seq_bin_by_insert_length_.resize(num_bins);

	// Prepare bias fitting
	bias_calc_params_.resize( maximum_insert_length * max_bins_queued_for_bias_calc_ );

	if( BiasCalculationVectors::dispersion_info_file ){
		ref.RefSeqsInNxx(ref_seq_in_nxx_, 0.0); // Get only the largest reference sequence
	}
	else{
		ref.RefSeqsInNxx(ref_seq_in_nxx_, BiasCalculationVectors::nxx_ref_seqs_);
	}

	// Temporary bias fit result vectors
	auto max_fits = ref.NumRefSeqBinsInNxx(ref_seq_in_nxx_, max_ref_seq_bin_length_) * BiasCalculationVectors::num_fits_insert_length_;
	for(uint16_t gc=tmp_gc_bias_.size(); gc--;){
		tmp_gc_bias_.at(gc).resize(max_fits, nan("0"));
	}
	for(uint16_t sur=tmp_sur_bias_.size(); sur--;){
		tmp_sur_bias_.at(sur).resize(max_fits, nan("0"));
	}
	for(uint16_t par=tmp_dispersion_parameters_.size(); par--;){
		tmp_dispersion_parameters_.at(par).resize(max_fits, nan("0"));
	}

	// Initialize bias calculation queue information
	params_left_for_calculation_ = ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1) * maximum_insert_length;
	for( uint16_t bin = 0; bin < current_bias_param_.size(); ++bin){
		current_bias_param_.at(bin) = maximum_insert_length+1; // So no thread tries to calculate on the unfilled bins
		finished_bias_calcs_.at(bin) = 0;
		claimed_bias_bins_.at(bin).clear();
	}

	// Creating files for additional output
	if( BiasCalculationVectors::parameter_info_file ){
		ofstream myfile;
		myfile.open(BiasCalculationVectors::parameter_info_file);
		myfile << "Fit, RefSeq, InsertLength, Counts, Sites, FunctionCallsSampled, FunctionCalls, LogLikelihood, DispersionA, DispersionB";
		for(auto gc = 0; gc < 101; ++gc){
			myfile << ", " << "GCbias" << gc << ", " << "GCcount" << gc << ", " << "GCsites" << gc;
		}
		for(auto k = 0; k < BiasCalculationVectors::gc_spline_df_; ++k){
			myfile << ", " << "GCknot" << k;
		}
		for(auto sur = 0; sur < 4*Reference::num_surrounding_blocks_*Reference::surrounding_range_; ++sur){
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

	if( BiasCalculationVectors::dispersion_info_file ){
		ofstream myfile;
		myfile.open(BiasCalculationVectors::dispersion_info_file);
		myfile << "mean_fit, prediction_mean, sample_mean, num_sites, max_duplication, counts1, counts2, counts3, calls_prediction, dispersion_prediction, calls_sample, dispersion_sample" << std::endl;
		myfile.close();
	}
}

void FragmentDistributionStats::SortFragmentSites(uint32_t ref_seq_bin, vector<uint32_t> &num_sites_per_insert_length){
	// Sort sites into insert_length bins
	// Start by getting size of each bin
	num_sites_per_insert_length.clear();
	num_sites_per_insert_length.resize(tmp_insert_lengths_.size(), 0);
	for( uint64_t site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
		++num_sites_per_insert_length.at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).second);
	}

	uint32_t max_used_length = tmp_insert_lengths_.size();
	while(--max_used_length && 0 == num_sites_per_insert_length.at(max_used_length) ); // Ignore empty bins at the end
	if(max_used_length){
		// Reserve needed space
		fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).resize(max_used_length+1);
		for(uint32_t insert_length = 1; insert_length < max_used_length+1; ++insert_length ){
			fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).reserve(num_sites_per_insert_length.at(insert_length));
		}

		// Sort into bins
		for( uint64_t site_id=0; site_id < fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_bin); ++site_id ){
			fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).second).push_back(fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).at(site_id).first);
		}
	}

	// Free space that is not needed anymore
	fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).clear();
	fragment_sites_by_ref_seq_bin_.at(ref_seq_bin).shrink_to_fit();
}

void FragmentDistributionStats::UpdateBiasCalculationParams(uint32_t ref_seq_bin, uint32_t queue_spot, vector<pair<uint64_t, pair<uint32_t, uint32_t>>> &tmp_params, mutex &print_mutex){
	auto qbin_size=bias_calc_params_.size()/max_bins_queued_for_bias_calc_;
	if(fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size()){
		// Select ref_seq, insert_length bins to use for bias calculation
		if(ref_seq_in_nxx_.at( GetRefSeqId(ref_seq_bin) )){
			// Sort them by counts and only take top num_fits_insert_length_
			tmp_params.clear();
			for(uint32_t insert_length=1; insert_length < fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); ++insert_length){
				if( fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size() ){
					tmp_params.push_back({fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size(), {ref_seq_bin, insert_length}});
				}
				else{
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Clear(ref_seq_bin); // If we won't set anything to it later we have to clear it from it's previous content
				}
			}

			sort(tmp_params.begin(), tmp_params.end(), std::greater<pair<uint64_t, pair<uint32_t, uint32_t>>>());

			if( BiasCalculationVectors::dispersion_info_file ){
				bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(0).second.second-1).Set(tmp_params.at(0).second.first, tmp_params.at(0).second.second, calculate_bias_);
				for( uint32_t n=1; n < tmp_params.size(); ++n ){
					bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, false);
				}
			}
			else{
				for( uint32_t n=0; n < tmp_params.size(); ++n ){
					if(n < BiasCalculationVectors::num_fits_insert_length_){
						bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, calculate_bias_);
					}
					else{
						bias_calc_params_.at(queue_spot*qbin_size + tmp_params.at(n).second.second-1).Set(tmp_params.at(n).second.first, tmp_params.at(n).second.second, false);
					}
				}

				if(tmp_params.size() < BiasCalculationVectors::num_fits_insert_length_){
					params_fitted_ += BiasCalculationVectors::num_fits_insert_length_-tmp_params.size(); // Do not have to be fitted, so are already done
					if(!BiasCalculationVectors::parameter_info_file){
						lock_guard<mutex> lock(print_mutex);
						printInfo << "Finished " << static_cast<uint16_t>(Percent(params_fitted_, tmp_gc_bias_.at(0).size())) << "% of the bias fits.\n";
					}
				}
			}
		}
		else{
			for(uint32_t insert_length=1; insert_length < fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); ++insert_length){
				if( fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length).size() ){
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Set(ref_seq_bin, insert_length, false);
				}
				else{
					bias_calc_params_.at(queue_spot*qbin_size + insert_length-1).Clear(ref_seq_bin); // If we won't set anything to it later we have to clear it from it's previous content
				}
			}
		}

		// Clear all the other insert_length, so no old values from other ref bins are in there
		for(uint32_t insert_length=fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).size(); insert_length <= qbin_size; ++insert_length){
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

	for( auto frag_length=max(static_cast<uint64_t>(1),insert_lengths_.from()); frag_length < insert_lengths_.to(); ++frag_length){
		if(insert_lengths_.at(frag_length)){
			for( auto ref_id=abundance_.size(); ref_id--; ){
				if(abundance_.at(ref_id) && reference.SequenceLength(ref_id) >= insert_lengths_.to()+2*reference.MinDistToRefSeqEnds()){
					// Add parameters for later calculation
					params.push_back({ref_id, frag_length});
				}
			}
		}
	}
	params.shrink_to_fit();
}

void FragmentDistributionStats::CountDuplicates(FragmentDuplicationStats &duplications, const BiasCalculationParamsSplitSeqs &params, const Reference &reference){
	duplications.AddDuplicates( fragment_sites_by_ref_seq_bin_by_insert_length_.at(params.ref_seq_bin_).at(params.fragment_length_), GetRefSeqId(params.ref_seq_bin_), params.fragment_length_, reference  );
}

void FragmentDistributionStats::AddFragmentsToSites(vector<FragmentSite> &sites, const vector<uint32_t> &fragment_positions, uint32_t min_dist_to_ref_seq_ends){
	for(auto pos : fragment_positions){
		if(pos%2){
			++sites.at(pos/2-min_dist_to_ref_seq_ends).count_forward_;
		}
		else{
			++sites.at(pos/2-min_dist_to_ref_seq_ends).count_reverse_;
		}
	}
}

void FragmentDistributionStats::CalculateBiasByBin(BiasCalculationVectors &tmp_calc, const Reference &reference, FragmentDuplicationStats &duplications, uint32_t ref_seq_bin, uint32_t insert_length){
	// Prepare fits
	auto ref_seq_id = GetRefSeqId(ref_seq_bin);
	uint32_t ref_start = (ref_seq_bin-ref_seq_start_bin_.at(ref_seq_id))*RefSeqSplitLength(ref_seq_id, reference);
	uint32_t ref_end = reference.SequenceLength(ref_seq_id);
	if(ref_seq_bin+1 < ref_seq_start_bin_.at(ref_seq_id+1)){ // Not the last bin of the reference
		ref_end = (ref_seq_bin+1-ref_seq_start_bin_.at(ref_seq_id))*RefSeqSplitLength(ref_seq_id, reference);
	}

	reference.GetFragmentSites(tmp_calc.sites_, ref_seq_id, insert_length, ref_start, ref_end);
	AddFragmentsToSites(tmp_calc.sites_, fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_bin).at(insert_length), max(static_cast<uint32_t>(reference.MinDistToRefSeqEnds()), ref_start) );
	duplications.AddSites(tmp_calc.sites_, insert_length);
	tmp_calc.GetCounts();
	tmp_calc.RemoveUnnecessarySites();
	tmp_calc.GetLogLikeBase();

	// Fit poisson with independent gc bins
	tmp_calc.fit_pars_.clear();
	if(BiasCalculationVectors::sur_mult_){
		tmp_calc.fit_pars_.resize(tmp_calc.sur_bias_.size(), 1.0);
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else if(BiasCalculationVectors::sur_sum_){
		tmp_calc.fit_pars_.resize(tmp_calc.sur_bias_.size()+1, 0.0);
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::upper_bound_, 1);
		tmp_calc.fit_pars_.at(0) = 1.0;
	}

	tmp_calc.func_calls_pois_ = 0;

	tmp_calc.optimizer_spline_.set_max_objective(BiasCalculationVectors::LogLikeGcSpline, reinterpret_cast<void*>(&tmp_calc));
	tmp_calc.optimizer_poisson_.set_max_objective(BiasCalculationVectors::LogLikelihoodPoisson, reinterpret_cast<void*>(&tmp_calc));

	double lower_bound;
	if(BiasCalculationVectors::sur_mult_){
		lower_bound = BiasCalculationVectors::base_value_;
	}
	else if(BiasCalculationVectors::sur_sum_){
		lower_bound = BiasCalculationVectors::lower_bound_;
	}
	const double upper_bound = BiasCalculationVectors::upper_bound_;
	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), lower_bound);
	if(BiasCalculationVectors::sur_mult_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.bounds_, 0.0, 0);
	}
	tmp_calc.optimizer_poisson_.set_lower_bounds(tmp_calc.bounds_);

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), upper_bound);
	if(BiasCalculationVectors::sur_mult_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else if(BiasCalculationVectors::sur_sum_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::upper_bound_, 1);
	}
	tmp_calc.optimizer_poisson_.set_upper_bounds(tmp_calc.bounds_);

	tmp_calc.converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == tmp_calc.optimizer_poisson_.optimize(tmp_calc.fit_pars_, tmp_calc.loglike_pois_) ){
			//printWarn << "Optimizer reached maximum number of likelihood calculations for seq " << ref_seq_bin << " and length " << insert_length << std::endl;
			tmp_calc.func_calls_pois_ += 10000;
			tmp_calc.converged_ = false;
		}

	}
	catch(const std::exception& e){
		//printWarn << "Poisson optimizer crashed for seq " << ref_seq_bin << " and length " << insert_length << " with: " << e.what() << std::endl;
		tmp_calc.func_calls_pois_ += 20000;
		tmp_calc.converged_ = false;
	}

	// Fit gc spline with fixed surrounding bias assuming poisson
	tmp_calc.GetGCSpline();

	for(uint16_t sur=tmp_calc.sur_bias_.size(); sur--; ){
		tmp_calc.sur_bias_pois_.at(sur) = tmp_calc.sur_bias_.at(sur);
	}

	// Fit negative binomial to downsampled sites with fixed biases to estimate dispersion parameters
	tmp_calc.optimizer_nbinom_.set_max_objective(BiasCalculationVectors::LogLikelihoodNbinom, reinterpret_cast<void*>(&tmp_calc));

	for(auto k=0; k < tmp_calc.gc_spline_pars_.size(); ++k){
		tmp_calc.fit_pars_.push_back( tmp_calc.gc_spline_pars_.at(k) );
	}
	tmp_calc.fit_pars_.push_back( 1.0 ); // alpha
	tmp_calc.fit_pars_.push_back( 1.0 ); // beta

	tmp_calc.func_calls_nbinom_ = 0;

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), lower_bound);
	if(BiasCalculationVectors::sur_mult_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.bounds_, 0.0, 0);
		for(uint16_t k=tmp_calc.sur_count_.size()+2; k < tmp_calc.fit_pars_.size()-2; ++k){
			tmp_calc.bounds_.at(k) = BiasCalculationVectors::lower_bound_;
		}
	}
	else if(BiasCalculationVectors::sur_sum_){
		tmp_calc.bounds_.at(tmp_calc.sur_count_.size()+1) = BiasCalculationVectors::base_value_;
	}
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 1) = BiasCalculationVectors::base_value_;
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 2) = BiasCalculationVectors::base_value_;
	tmp_calc.optimizer_nbinom_.set_lower_bounds(tmp_calc.bounds_);

	tmp_calc.bounds_.clear();
	tmp_calc.bounds_.resize(tmp_calc.fit_pars_.size(), upper_bound);
	if(BiasCalculationVectors::sur_mult_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, 0.0, 0);
	}
	else if(BiasCalculationVectors::sur_sum_){
		tmp_calc.DeactivateZeroCounts(tmp_calc.fit_pars_, BiasCalculationVectors::lower_bound_, 1);
	}
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 1) = BiasCalculationVectors::upper_bound_;
	tmp_calc.bounds_.at(tmp_calc.bounds_.size() - 2) = BiasCalculationVectors::upper_bound_;
	tmp_calc.optimizer_nbinom_.set_upper_bounds(tmp_calc.bounds_);

	tmp_calc.converged_ = true;
	try{
		if( nlopt::MAXEVAL_REACHED == tmp_calc.optimizer_nbinom_.optimize(tmp_calc.fit_pars_, tmp_calc.loglike_nbinom_) ){
			//printWarn << "Optimizer reached maximum number of likelihood calculations for seq " << ref_seq_bin << " and length " << insert_length << std::endl;
			tmp_calc.func_calls_nbinom_ += 10000;
			tmp_calc.converged_ = false;
		}
	}
	catch(const std::exception& e){
		//printWarn << "Nbinom optimizer crashed for seq " << ref_seq_bin << " and length " << insert_length << " with: " << e.what() << std::endl;
		tmp_calc.func_calls_nbinom_ += 20000;
		tmp_calc.converged_ = false;
	}

	tmp_calc.PostProcessFit();
}

void FragmentDistributionStats::AcquireBiases(const BiasCalculationVectors &calc, const BiasCalculationParamsSplitSeqs &params, mutex &print_mutex){
	if(calc.converged_){
		uint32_t cur_res = current_bias_result_++;

		// Store biases in vectors from where a single value will be chosen later
		for(uint16_t gc=101; gc--; ){
			tmp_gc_bias_.at(gc).at(cur_res) = calc.gc_bias_.at(gc);
		}
		for(uint16_t sur=calc.sur_bias_.size(); sur--; ){
			tmp_sur_bias_.at(sur).at(cur_res) = calc.sur_bias_.at(sur);
		}
		for(uint16_t par=calc.dispersion_.size(); par--; ){
			tmp_dispersion_parameters_.at(par).at(cur_res) = calc.dispersion_.at(par);
		}
	}

	uint32_t pars_fitted = ++params_fitted_;
	if(!BiasCalculationVectors::parameter_info_file){
		if( tmp_gc_bias_.at(0).size() > 20 && !(pars_fitted%(tmp_gc_bias_.at(0).size()/20)) ){
			lock_guard<mutex> lock(print_mutex);
			printInfo << "Finished " << static_cast<uint16_t>(Percent(pars_fitted, tmp_gc_bias_.at(0).size())) << "% of the bias fits.\n";
		}
	}

	// Print information
	if(BiasCalculationVectors::parameter_info_file){
		lock_guard<mutex> lock(print_mutex);

		auto ref_seq_id = GetRefSeqId(params.ref_seq_bin_);

		std::cout << ref_seq_id << ' ' << params.fragment_length_ << ' ' << calc.total_counts_ << ' ' << calc.total_sites_ << ' ' << calc.func_calls_pois_ << ' ' << calc.loglike_pois_ << ' ' << calc.func_calls_spline_ << ' ' << calc.loglike_spline_ << ' ' << calc.func_calls_nbinom_sampled_ << ' ' << calc.func_calls_nbinom_ << ' ' << calc.loglike_nbinom_ << std::endl;

		ofstream myfile;
		myfile.open(BiasCalculationVectors::parameter_info_file, ofstream::out | ofstream::app);
		myfile << "poisson, " << ref_seq_id << ", " << params.fragment_length_ << ", " << calc.total_counts_ << ", " << calc.total_sites_ << ", 0, " << calc.func_calls_pois_ << ", " << calc.loglike_pois_ << ", " << 0 << ", " << 0;
		for(auto gc = 0; gc < calc.gc_bias_pois_.size(); ++gc){
			myfile << ", " << calc.gc_bias_pois_.at(gc) << ", " << calc.gc_count_.at(gc) << ", " << calc.gc_sites_.at(gc);
		}
		for(auto k = 0; k < calc.gc_knots_.size(); ++k){
			myfile << ", " << calc.gc_knots_.at(k);
		}
		for(auto sur = 0; sur < calc.sur_bias_pois_.size(); ++sur){
			myfile << ", " << calc.sur_bias_pois_.at(sur) << ", " << calc.sur_count_.at(sur) << ", " << calc.sur_sites_.at(sur);
		}
		myfile << std::endl;

		myfile << "gcspline, " << ref_seq_id << ", " << params.fragment_length_ << ", " << calc.total_counts_ << ", " << calc.total_sites_ << ", 0, " << calc.func_calls_spline_ << ", " << calc.loglike_spline_ << ", " << 0 << ", " << 0;
		for(auto gc = 0; gc < 101; ++gc){
			myfile << ", " << calc.gc_bias_spline_.at(gc) << ", " << calc.gc_count_.at(gc) << ", " << calc.gc_sites_.at(gc);
		}
		for(auto k = 0; k < calc.gc_knots_.size(); ++k){
			myfile << ", " << calc.gc_knots_.at(k);
		}
		for(auto sur = 0; sur < calc.sur_bias_pois_.size(); ++sur){
			myfile << ", " << calc.sur_bias_pois_.at(sur) << ", " << calc.sur_count_.at(sur) << ", " << calc.sur_sites_.at(sur);
		}
		myfile << std::endl;

		myfile << "nbinom, " << ref_seq_id << ", " << params.fragment_length_ << ", " << calc.total_counts_ << ", " << calc.total_sites_ << ", " <<  calc.func_calls_nbinom_sampled_ << ", " << calc.func_calls_nbinom_ << ", " << calc.loglike_nbinom_ << ", " << calc.dispersion_.at(0) << ", " << calc.dispersion_.at(1);
		for(auto gc = 0; gc < 101; ++gc){
			myfile << ", " << calc.gc_bias_.at(gc) << ", " << calc.gc_count_.at(gc) << ", " << calc.gc_sites_.at(gc);
		}
		for(auto k = 0; k < calc.gc_knots_.size(); ++k){
			myfile << ", " << calc.gc_knots_.at(k);
		}
		for(auto sur = 0; sur < calc.sur_bias_.size(); ++sur){
			myfile << ", " << calc.sur_bias_.at(sur) << ", " << calc.sur_count_.at(sur) << ", " << calc.sur_sites_.at(sur);
		}
		myfile << std::endl;
		myfile.close();
	}
}

bool FragmentDistributionStats::StoreBias(){
	printInfo << "Acquired " << static_cast<uint16_t>(Percent(current_bias_result_, tmp_gc_bias_.at(0).size())) << "% of the bias fits.\n";

	// Remove the nan at the end for not converged fits, this will affect all gc and sur biases in the same way so do it only once
	uint32_t used_size = tmp_gc_bias_.at(0).size();
	while(used_size-- && isnan( tmp_gc_bias_.at(0).at( used_size ) ));
	++used_size; // Go from index to size

	if( 0 == used_size ){
		printErr << "No bias fit converged. Cannot continue" << std::endl;
		return false;
	}

	// GC bias
	// Calculate Median
	for(uint16_t gc=0; gc<tmp_gc_bias_.size(); ++gc){
		tmp_gc_bias_.at(gc).resize(used_size);

		sort(tmp_gc_bias_.at(gc).begin(), tmp_gc_bias_.at(gc).end());
		gc_fragment_content_bias_[gc] = tmp_gc_bias_.at(gc).at( tmp_gc_bias_.at(gc).size()/2 );

		tmp_gc_bias_.at(gc).clear();
		tmp_gc_bias_.at(gc).shrink_to_fit();
	}

	// Calculate Max
	double max_val(*max_element(gc_fragment_content_bias_.begin(), gc_fragment_content_bias_.end()));
	// Normalize so that maximum bias is 1.0
	for(uint16_t gc=gc_fragment_content_bias_.size(); gc--;){
		gc_fragment_content_bias_[gc] /= max_val;
	}

	// Surrounding bias
	// Calculate Median
	std::array<double, 4*Reference::num_surrounding_blocks_*Reference::surrounding_range_> sur_median;
	for(uint16_t sur=0; sur<tmp_sur_bias_.size(); ++sur){
		tmp_sur_bias_.at(sur).resize(used_size);

		sort(tmp_sur_bias_.at(sur).begin(), tmp_sur_bias_.at(sur).end());
		sur_median[sur] = tmp_sur_bias_.at(sur).at( tmp_sur_bias_.at(sur).size()/2 );

		tmp_sur_bias_.at(sur).clear();
		tmp_sur_bias_.at(sur).shrink_to_fit();
	}

	CombineSurroundingPositions(fragment_surroundings_bias_, sur_median);

	// Dispersion parameters
	// Calculate Median
	for(uint16_t par=0; par<tmp_dispersion_parameters_.size(); ++par){
		tmp_dispersion_parameters_.at(par).resize(used_size);

		sort(tmp_dispersion_parameters_.at(par).begin(), tmp_dispersion_parameters_.at(par).end());
		dispersion_parameters_[par] = tmp_dispersion_parameters_.at(par).at( tmp_dispersion_parameters_.at(par).size()/2 );

		tmp_dispersion_parameters_.at(par).clear();
		tmp_dispersion_parameters_.at(par).shrink_to_fit();
	}

	return true;
}

void FragmentDistributionStats::AddNewBiasCalculations(uint32_t still_needed_ref_bin, ThreadData &thread, mutex &print_mutex){
	uint32_t ref_seq_bin = num_handled_reference_sequence_bins_;
	while( ref_seq_bin < still_needed_ref_bin ){
		// Check if a spot in the queue is free
		uint32_t queue_spot = 0;
		while( queue_spot < max_bins_queued_for_bias_calc_ && claimed_bias_bins_.at(queue_spot).test_and_set() ){
			++queue_spot;
		}

		if(queue_spot < max_bins_queued_for_bias_calc_){
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
	auto queue_bin_size = bias_calc_params_.size()/max_bins_queued_for_bias_calc_;
	for(uint16_t queue_bin=0; queue_bin < max_bins_queued_for_bias_calc_; ++queue_bin ){
		auto cur_par = queue_bin*queue_bin_size + current_bias_param_.at(queue_bin)++;

		for(; cur_par < queue_bin*queue_bin_size + queue_bin_size; cur_par = queue_bin*queue_bin_size + current_bias_param_.at(queue_bin)++){
			--params_left_for_calculation_;

			// Execute calculations
			if( bias_calc_params_.at(cur_par).fragment_length_ ){
				if( bias_calc_params_.at(cur_par).bias_calculation_ ){
					CalculateBiasByBin(thread_values, reference, duplications, bias_calc_params_.at(cur_par).ref_seq_bin_, bias_calc_params_.at(cur_par).fragment_length_);
					AcquireBiases(thread_values, bias_calc_params_.at(cur_par), print_mutex);
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

void FragmentDistributionStats::HandleReferenceSequencesUntil(uint32_t still_needed_ref_bin, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, mutex &print_mutex){
	AddNewBiasCalculations(still_needed_ref_bin, thread, print_mutex);
	ExecuteBiasCalculations( reference, duplications, thread.bias_calc_vects_, print_mutex );
}

void FragmentDistributionStats::BiasSumThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<vector<BiasCalculationParams>::size_type> &current_param,
		vector<double> &insert_length_sum,
		vector<double> &ref_seq_sum,
		vector<vector<VectorAtomic<uint64_t>>> &site_count_by_insert_length_gc,
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
			printInfo << "Finished " << static_cast<uint16_t>(Percent(cur_par, params.size())) << "% of the bias sums.\n";
		}
	}

	result_mutex.lock();
	for(uint32_t l=len.size(); l--; ){
		insert_length_sum.at(l) += len.at(l);
	}
	for(uint32_t r=ref.size(); r--; ){
		ref_seq_sum.at(r) += ref.at(r);
	}
	result_mutex.unlock();
}

void FragmentDistributionStats::BiasNormalizationThread(
		const FragmentDistributionStats &self,
		const Reference &reference,
		const vector<BiasCalculationParams> &params,
		atomic<vector<BiasCalculationParams>::size_type> &current_param,
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

uint32_t FragmentDistributionStats::CreateRefBins( const Reference &ref, uint64_t max_ref_seq_bin_size ){
	max_ref_seq_bin_length_ = max_ref_seq_bin_size;

	ref_seq_start_bin_.resize(ref.NumberSequences()+1); // Add one to give number of bins at the end
	uint32_t start_id = 0;
	ref_seq_start_bin_.at(0) = 0;
	for(uint32_t ref_seq=1; ref_seq < ref_seq_start_bin_.size(); ++ref_seq){
		start_id += ref.SequenceLength(ref_seq-1)/max_ref_seq_bin_length_ + 1;
		ref_seq_start_bin_.at(ref_seq) = start_id;
	}

	return ref_seq_start_bin_.at(ref_seq_start_bin_.size()-1);
}

void FragmentDistributionStats::Prepare( const Reference &ref, uint32_t maximum_insert_length, const vector<uint64_t> &reads_per_ref_seq_bin ){
	surrounding_range_ = ref.SurroundingRange();

	tmp_abundance_.resize(ref.NumberSequences());
	tmp_insert_lengths_.resize(maximum_insert_length+1);
	tmp_gc_fragment_content_.resize(101);
	for( auto &tmp_sur : tmp_fragment_surroundings_){
		tmp_sur.resize(Reference::SurroundingSize());
	}

	for(auto strand=2; strand--; ){
		for(auto nuc=4; nuc--; ){
			tmp_outskirt_content_.at(strand).at(nuc).resize(2*outskirt_range_);
		}
	}

	PrepareBiasCalculation( ref, maximum_insert_length, reads_per_ref_seq_bin );
}

void FragmentDistributionStats::FillInOutskirtContent( const Reference &reference, const BamAlignmentRecord &record_start, uint32_t fragment_end_pos ){
	bool reversed_fragment = hasFlagLast(record_start); // record is forward, so if it is sequenced last, fragment is reversed
	auto ref_pos = record_start.beginPos;
	uint16_t outskirt_pos;
	FunctorComplement<Dna5> complement;
	auto &ref_seq = reference.ReferenceSequence(record_start.rID);

	// Fill in fragment start
	array<int32_t, Reference::num_surrounding_blocks_> surroundings;
	reference.ForwardSurroundingWithN(surroundings, record_start.rID, ref_pos);
	for( auto block = surroundings.size(); block--; ){
		if(0 <= surroundings.at(block)){
			++(tmp_fragment_surroundings_.at(block).at( surroundings.at(block) ));
		}
	}

	if(reversed_fragment){
		// Reverse strand: Fill in postfragment-content
		if( ref_pos < outskirt_range_ ){
			outskirt_pos = outskirt_range_ + ref_pos;
			ref_pos = 0;
		}
		else{
			outskirt_pos = 2*outskirt_range_;
			ref_pos -= outskirt_range_;
		}

		for(; outskirt_pos-- > outskirt_range_; ){
			if( 4 >static_cast<uint16_t>(ref_seq[ref_pos]) ){
				++tmp_outskirt_content_.at(true).at(complement(ref_seq[ref_pos])).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
	else{
		// Forward strand: Fill in prefragment-content
		if( ref_pos < outskirt_range_ ){
			outskirt_pos = outskirt_range_ - ref_pos;
			ref_pos = 0;
		}
		else{
			outskirt_pos = 0;
			ref_pos -= outskirt_range_;
		}

		for(; outskirt_pos < outskirt_range_; ++outskirt_pos){
			if( 4 >static_cast<uint16_t>(ref_seq[ref_pos]) ){
				++tmp_outskirt_content_.at(false).at(ref_seq[ref_pos]).at(outskirt_pos);
			}
			++ref_pos;
		}
	}

	ref_pos = fragment_end_pos;

	// Fill in fragment end
	reference.ReverseSurroundingWithN(surroundings, record_start.rNextId, ref_pos-1);
	for( auto block = surroundings.size(); block--; ){
		if(0 <= surroundings.at(block)){
			++(tmp_fragment_surroundings_.at(block).at( surroundings.at(block) ));
		}
	}

	if(reversed_fragment){
		// Reverse strand: Fill in prefragment-content
		outskirt_pos = outskirt_range_;
		for(;  outskirt_pos-- && ref_pos < reference.SequenceLength(record_start.rNextId);){
			if( 4 >static_cast<uint16_t>(ref_seq[ref_pos]) ){
				++tmp_outskirt_content_.at(true).at(complement(ref_seq[ref_pos])).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
	else{
		// Forward strand:Fill in postfragment-content
		for(; outskirt_pos < 2*outskirt_range_ && ref_pos < reference.SequenceLength(record_start.rNextId); ++outskirt_pos){
			if( 4 >static_cast<uint16_t>(ref_seq[ref_pos]) ){
				++tmp_outskirt_content_.at(false).at(ref_seq[ref_pos]).at(outskirt_pos);
			}
			++ref_pos;
		}
	}
}

void FragmentDistributionStats::HandleReferenceSequencesUntil(uint32_t still_needed_reference_sequence, uint32_t still_needed_position, ThreadData &thread, const Reference &reference, FragmentDuplicationStats &duplications, mutex &print_mutex){
	// Check if new biases can be added for calculation
	auto still_needed_ref_bin = ref_seq_start_bin_.at(still_needed_reference_sequence);
	still_needed_ref_bin += still_needed_position / RefSeqSplitLength(still_needed_reference_sequence, reference);
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

void FragmentDistributionStats::Shrink(){
	insert_lengths_.Shrink();
	gc_fragment_content_.Shrink();
	for( int direction=2; direction--; ){
		for( int base=4; base--; ){
			outskirt_content_.at(direction).at(base).Shrink();
		}
	}
}

void FragmentDistributionStats::Finalize(){
	// Acquire data from temporary vectors
	Acquire(abundance_, tmp_abundance_);
	insert_lengths_.Acquire(tmp_insert_lengths_);
	gc_fragment_content_.Acquire(tmp_gc_fragment_content_);
	for(auto i=fragment_surroundings_.size(); i--; ){
		Acquire(fragment_surroundings_.at(i), tmp_fragment_surroundings_.at(i));
	}

	for(auto strand=2; strand--; ){
		for(auto nuc=4; nuc--; ){
			outskirt_content_.at(strand).at(nuc).Acquire(tmp_outskirt_content_.at(strand).at(nuc));
		}
	}
}

void FragmentDistributionStats::CalculateInsertLengthAndRefSeqBias(const Reference &reference, uint16_t num_threads, vector<vector<VectorAtomic<uint64_t>>> &site_count_by_insert_length_gc){
	// Sum up the other biases for the insert length, reference sequences pairs
	atomic<vector<BiasCalculationParams>::size_type> current_param(0);
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

bool FragmentDistributionStats::FinalizeBiasCalculation(const Reference &reference, uint16_t num_threads, FragmentDuplicationStats &duplications){

	vector<vector<VectorAtomic<uint64_t>>> site_count_by_insert_length_gc;
	if(calculate_bias_){
		if( !StoreBias() ){
			return false;
		}

		CalculateInsertLengthAndRefSeqBias(reference, num_threads, site_count_by_insert_length_gc);
	}

	// Free not needed memory at the end
	fragment_sites_by_ref_seq_bin_.clear();
	fragment_sites_by_ref_seq_bin_.shrink_to_fit();
	fragment_sites_by_ref_seq_bin_by_insert_length_.clear();
	fragment_sites_by_ref_seq_bin_by_insert_length_.shrink_to_fit();
	duplications.FinalizeDuplicationVector(site_count_by_insert_length_gc);

	printInfo << "Finished bias calculation" << std::endl;

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
		uniform_int_distribution<uint32_t> rdist(0,old_bias_.size());
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
			bool error = false;

			// Check which reference sequences are needed and store their ref_seq_bias_ position
			unordered_map<string, uint32_t> ref_seq_ids;
			string tmp_id;
			for(uint32_t seq=0; seq < ref.NumberSequences(); ++seq){
				tmp_id = toCString(ref.ReferenceId(seq));
				ref_seq_ids.emplace(tmp_id.substr(0,tmp_id.find(' ')), seq);
			}

			// Read in information from file (One line corresponds to one reference sequence [identifier bias])
			double bias;
			string line;
			uint32_t nline;
			bool empty_line = false; // Allow only a single empty line at the end of the file
			std::vector<bool> bias_found;
			bias_found.resize(ref_seq_bias_.size(), false);
			while( getline(fbias,line) ){
				if(empty_line){
					printErr << "Error reading in reference sequence biases. Line " << nline << " is empty in file " << bias_file << std::endl;
					error = true;
				}
				else{
					++nline;

					if(0 == line.size()){
						empty_line = true;
					}
					else{
						auto bias_sep = line.find_last_of(" \t");
						if( bias_sep == string::npos ){
							printErr << "Error reading in reference sequence biases. Line " << nline << " does not have a field separator in file " << bias_file << std::endl;
							error = true;
						}
						else{
							try{
								bias = stod( line.substr(bias_sep+1) );
							}
							catch(const exception &e){
								printErr <<  "Could not convert bias starting at postion " << bias_sep+1 << " in line " << nline << " in file " << bias_file << " into double\nException thrown: " << e.what() << std::endl;
								error = true;
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

			for(uint32_t seq=0; seq < ref.NumberSequences(); ++seq){
				if(!bias_found.at(seq)){
					printErr <<  "Could not find bias for reference sequence " << ref.ReferenceId(seq) << " in file " << bias_file << std::endl;
					error = true;
				}
			}

			if(error){
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

double FragmentDistributionStats::CalculateBiasNormalization(double &max_bias, const Reference &reference, uint16_t num_threads, uint64_t total_reads) const{
	atomic<vector<BiasCalculationParams>::size_type> current_param(0);
	vector<BiasCalculationParams> params;
	FillParams(params, reference);

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

uint32_t FragmentDistributionStats::NegativeBinomial(double p, double r, double probability_chosen){
	// (1-p)^r * mult_{i=1}^{k}[p * (r+i-1)/i]
	double probability_count( pow(1-p, r) ), probability_sum(probability_count);
	uint32_t count(0);

	while(probability_sum < probability_chosen){
		probability_count *= p * ((r-1) / ++count + 1);
		probability_sum += probability_count;
	}

	return count;
}

uint32_t FragmentDistributionStats::GetFragmentCounts(const Reference &reference, double bias_normalization, uint32_t ref_seq_id, uint32_t fragment_length, uint16_t gc, const array<uint32_t, 3> &fragment_start, const array<uint32_t, 3> &fragment_end, double probability_chosen, double non_zero_threshold) const{
	if(probability_chosen < non_zero_threshold){
		return 0;
	}
	else{
		array<double, 3> start_bias, end_bias;
		for(auto block = start_bias.size(); block--; ){
			start_bias.at(block) = fragment_surroundings_bias_.at(block).at(fragment_start.at(block));
			end_bias.at(block) = fragment_surroundings_bias_.at(block).at(fragment_end.at(block));
		}

		double bias( reference.Bias( ref_seq_bias_.at(ref_seq_id), insert_lengths_bias_[fragment_length], gc_fragment_content_bias_[gc], start_bias, end_bias ) );

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

void FragmentDistributionStats::NegativeBinomial(vector<uint32_t> &counts, vector<pair<double,uint16_t>> &probabilities_chosen, uint16_t current, double p, double r) const{
	// (1-p)^r * mult_{i=1}^{k}[p * (r+i-1)/i]
	double probability_count( pow(1-p, r) ), probability_sum(probability_count);
	uint32_t count(0);

	while(current < counts.size()){
		while(probability_sum < probabilities_chosen.at(current).first){
			probability_count *= p * ((r-1) / ++count + 1);
			probability_sum += probability_count;
		}

		counts.at( probabilities_chosen.at(current++).second ) = count;
	}
}

void FragmentDistributionStats::GetFragmentCounts(vector<uint32_t> &counts, vector<pair<double,uint16_t>> &probabilities_chosen, const Reference &reference, double bias_normalization, uint32_t ref_seq_id, uint32_t fragment_length, uint16_t gc, const array<uint32_t, 3> &fragment_start, const array<uint32_t, 3> &fragment_end, double non_zero_threshold) const{
	sort(probabilities_chosen.begin(), probabilities_chosen.end());
	uint16_t current = 0;

	while(current < counts.size() && probabilities_chosen.at(current).first < non_zero_threshold){
		counts.at( probabilities_chosen.at(current++).second ) = 0;
	}

	if(current < counts.size()){
		array<double, 3> start_bias, end_bias;
		for(auto block = start_bias.size(); block--; ){
			start_bias.at(block) = fragment_surroundings_bias_.at(block).at(fragment_start.at(block));
			end_bias.at(block) = fragment_surroundings_bias_.at(block).at(fragment_end.at(block));
		}

		double bias( reference.Bias( ref_seq_bias_.at(ref_seq_id), insert_lengths_bias_[fragment_length], gc_fragment_content_bias_[gc], start_bias, end_bias ) );

		if(0.0 < bias){
			double mean( bias * bias_normalization );
			double dispersion( Dispersion(mean) );

			NegativeBinomial(counts, probabilities_chosen, current, mean/(mean+dispersion), dispersion);
		}
		else{
			for(auto &count : counts){
				count = 0;
			}
		}
	}
}

void FragmentDistributionStats::PreparePlotting(){
	vector<double> separated_bias;
	SeparateSurroundingPositions(separated_bias, fragment_surroundings_bias_);

	for( auto &sur_vect : fragment_surrounding_bias_by_base_ ){
		sur_vect.resize(separated_bias.size()/4, 0.0);
	}

	// Split by nucleotide
	for( uint16_t i=0; i < separated_bias.size(); ++i){
		fragment_surrounding_bias_by_base_.at(i%4).at(i/4) = separated_bias.at(i);
	}
}
