#include "FragmentDistributionStatsTest.h"
using reseq::FragmentDistributionStatsTest;

#include <algorithm>
using std::max_element;
using std::min;
using std::max;
//include <array>
using std::array;
#include <cmath>
using std::round;
//include <mutex>
using std::mutex;
#include <numeric>
using std::accumulate;
//include <random>
using std::mt19937_64;
using std::uniform_real_distribution;
#include <string>
using std::string;
#include <thread>
using std::thread;
#include <vector>
using std::vector;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::length;

#include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::InvLogit2;
using reseq::utilities::Percent;
using reseq::utilities::SafePercent;

namespace { // Unnamed namespace so it only exists for this file
reseq::uintNumThreads test_with_num_threads = 4;
}

void FragmentDistributionStatsTest::Register(uintNumThreads num_threads){
	// Guarantees that library is included
	test_with_num_threads = num_threads;
}

void FragmentDistributionStatsTest::CreateTestObject(const Reference *ref){
	ASSERT_TRUE( test_ = new FragmentDistributionStats ) << "Could not allocate memory for FragmentDistributionStats object\n";

	test_->CreateRefBins(*ref, 100000000);

	vector<uintFragCount> reads_per_ref_seq_;
	reads_per_ref_seq_.resize(ref->NumberSequences(), 0);
	test_->Prepare(*ref, 100, reads_per_ref_seq_);

	test_->Finalize();
	test_->tmp_insert_lengths_.resize(101); // Reverted in Finalize, so it has to be done here
}

void FragmentDistributionStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void FragmentDistributionStatsTest::TestOutskirtContent(
		const FragmentDistributionStats &test,
		uintTempSeq template_segment,
		uintSeqLen at_pos,
		uintNucCount cont_a,
		uintNucCount cont_c,
		uintNucCount cont_g,
		uintNucCount cont_t,
		const char * context ){
	EXPECT_EQ(cont_a, test.outskirt_content_.at(template_segment).at(0)[at_pos]) << "outskirt_content_[" << static_cast<uintTempSeqPrint>(template_segment) << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_c, test.outskirt_content_.at(template_segment).at(1)[at_pos]) << "outskirt_content_[" << static_cast<uintTempSeqPrint>(template_segment) << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_g, test.outskirt_content_.at(template_segment).at(2)[at_pos]) << "outskirt_content_[" << static_cast<uintTempSeqPrint>(template_segment) << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_t, test.outskirt_content_.at(template_segment).at(3)[at_pos]) << "outskirt_content_[" << static_cast<uintTempSeqPrint>(template_segment) << "] position " << at_pos << " wrong for " << context << '\n';
}

void FragmentDistributionStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void FragmentDistributionStatsTest::CreateCoverageDataHelper(uintPercent gc_perc, uintSeqLen start_pos, uintSeqLen fragment_length, mt19937_64 &rgen){
	uintRefSeqId ref_seq_id(0);

	double bias;
	if(27 > gc_perc){
		bias = 50;
	}
	else if(31 > gc_perc){
		bias = 50+(gc_perc-26)*10;
	}
	else if(35 > gc_perc){
		bias = 20+(35-gc_perc)*20;
	}
	else{
		bias = 20;
	}
	bias /= 100;

	// Start bias
	double sur_bias = 0.0;
	switch( static_cast<char>(at(species_reference_.ReferenceSequence(0), start_pos)) ){
	case 'A':
		sur_bias -= 0.4;
		break;
	case 'C':
		sur_bias -= 0.4;
		break;
	case 'G':
		sur_bias += 0.4;
		break;
	case 'T':
		sur_bias += 0.4;
		break;
	}
	bias *= InvLogit2(sur_bias);

	// End bias (must be reverse of start bias)
	sur_bias = 0.0;
	switch( static_cast<char>(at(species_reference_.ReferenceSequence(0), start_pos+fragment_length-1)) ){
	case 'A':
		sur_bias += 0.4;
		break;
	case 'C':
		sur_bias += 0.4;
		break;
	case 'G':
		sur_bias -= 0.4;
		break;
	case 'T':
		sur_bias -= 0.4;
		break;
	}
	bias *= InvLogit2(sur_bias);

	// Add new fragments based on used coverage model
	uniform_real_distribution<double> zero_to_one;

	double mean( bias * 1.5 );
	double dispersion( BiasCalculationVectors::GetDispersion(bias, 0.5, 1.0) );
	double probability = mean/(mean+dispersion);

	auto count = FragmentDistributionStats::NegativeBinomial(probability, dispersion, zero_to_one(rgen));
	for(uintSeqLen entries=count; entries--; ){
		test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).push_back({start_pos << 1, fragment_length}); // Forward
		++test_->fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_id);
	}

	count = FragmentDistributionStats::NegativeBinomial(probability, dispersion, zero_to_one(rgen));
	for(uintSeqLen entries=count; entries--; ){
		test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).push_back({(start_pos << 1)+1, fragment_length}); // Reverse
		++test_->fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_id);
	}
}

void FragmentDistributionStatsTest::CreateCoverageData(uintSeqLen fragment_length){
	mt19937_64 rgen;
	rgen.seed( 201907171113 );

	uintRefSeqId ref_seq_id(0);
	test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).reserve(4*species_reference_.SequenceLength(ref_seq_id));
	uintSeqLen dist_ref_seq_ends(50);

	uintPercent gc_perc;
	uintSeqLen n_count(0);
	uintSeqLen gc = species_reference_.GCContentAbsolut(n_count, ref_seq_id, dist_ref_seq_ends, dist_ref_seq_ends+fragment_length);
	gc_perc = SafePercent(gc, fragment_length-n_count);

	CreateCoverageDataHelper(gc_perc, dist_ref_seq_ends, fragment_length, rgen);

	const Dna5String &ref_seq(species_reference_.ReferenceSequence(ref_seq_id));
	for( uintSeqLen start_pos=dist_ref_seq_ends; start_pos < length(ref_seq)-fragment_length-dist_ref_seq_ends; ){
		species_reference_.UpdateGC( gc, n_count, ref_seq_id, start_pos, start_pos+fragment_length );
		gc_perc = SafePercent(gc, fragment_length-n_count);

		CreateCoverageDataHelper(gc_perc, ++start_pos, fragment_length, rgen);
	}
}

template<size_t N> void FragmentDistributionStatsTest::CheckDrawnCounts(double bias, std::array<double, N> thresholds, const Surrounding &start_sur, const Surrounding &end_sur){
	double delta = 0.000001;
	double bias_normalization = 0.2;
	double other_bias_negation = 2.0*2.0*2.0*2.0*5.0;

	test_->gc_fragment_content_bias_[43] = bias * other_bias_negation;

	// Test GetFragmentCounts
	for(uintDupCount i = 0; i < thresholds.size(); ++i){
		EXPECT_EQ(i, test_->GetFragmentCounts(bias_normalization, 0, 367, 43, start_sur, end_sur, thresholds.at(i)-delta, 0.0)) << "Lower gate: " << i << ' ' << thresholds.at(i) << " for bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;
		EXPECT_EQ(i+1, test_->GetFragmentCounts(bias_normalization, 0, 367, 43, start_sur, end_sur, thresholds.at(i)+delta, 0.0)) << "Upper gate: " << i << ' ' << thresholds.at(i) << " for bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;
	}
	EXPECT_EQ(0, test_->GetFragmentCounts(1e-10, 0, 367, 43, start_sur, end_sur, thresholds.back(), 1.0)) << "Bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;

	// Test CalculateNonZeroThreshold
	auto non_zero_threshold = test_->CalculateNonZeroThreshold(bias_normalization, test_->gc_fragment_content_bias_[43]/other_bias_negation/bias_normalization);
	EXPECT_NEAR(thresholds[0], non_zero_threshold, delta) << "Bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;
	EXPECT_EQ(0, test_->GetFragmentCounts(bias_normalization-delta, 0, 367, 43, start_sur, end_sur, non_zero_threshold, 0.0)) << "Bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;
	EXPECT_EQ(0, test_->GetFragmentCounts(1e-10, 0, 367, 43, start_sur, end_sur, non_zero_threshold, 0.0)) << "Bias " << bias << " dispersion pars " << test_->dispersion_parameters_.at(0) << ' ' << test_->dispersion_parameters_.at(1) << std::endl;
}

void FragmentDistributionStatsTest::TestSrr490124Equality(const FragmentDistributionStats &test, const char *context){
	EXPECT_EQ(1, test.abundance_.size() ) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.abundance_.at(0) ) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{if(1==NR%2){store=$4}else{print store, $4, $4-store+length($10)}}'
	EXPECT_EQ(59, test.insert_lengths_.size()) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.insert_lengths_[126]) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.insert_lengths_[184]) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{if(1==NR%2){store=$4}else{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" store "-" $4+length($10)-1)}}' | seqtk seq | awk '(0==NR%2){len=length($1); print (gsub("G","",$1)+gsub("C","",$1))/len}'
	EXPECT_EQ(2, test.gc_fragment_content_.size()) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_fragment_content_[52]) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_fragment_content_[53]) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.bam -q 10 -f 144 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq) <(samtools view ecoli-SRR490124-4pairs.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq -r) | awk '(0==NR%2){print ">0", substr($0,1,10); print substr($0,length($0)-9,10); print ">1", substr($0,11,10); print substr($0,length($0)-19,10);print ">2", substr($0,21,10); print substr($0,length($0)-29,10);}' | seqtk seq -r | awk '{if(1==NR%2){sur=substr($0, 2, 1); print "start", sur, $2}else{print "end", sur, $0}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}' | sort
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(39619)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(1009062)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(561168)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(996771)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(307041)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(856654)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(607462)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(983881)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(354618)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(765926)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(333453)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(405510)) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(1==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4-20 "-" $4-1)}(0==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4+length($10) "-" $4+length($10)+19)}' | seqtk seq -r | awk '(2==NR%4){store=$0}(0==NR%4){print $0 store}' | awk '{for(pos=1;pos<=length($0);pos+=1){print substr($0,pos,1), pos-1}}' | sort | uniq -c | sort -k2,2 -k3,3n
	EXPECT_EQ(0, test.outskirt_content_.at(0).at(0).size()) << "SRR490124-4pairs outskirt_content_[0][0].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_.at(0).at(1).size()) << "SRR490124-4pairs outskirt_content_[0][1].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_.at(0).at(2).size()) << "SRR490124-4pairs outskirt_content_[0][2].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_.at(0).at(3).size()) << "SRR490124-4pairs outskirt_content_[0][3].size() not correct for " << context << '\n';
	EXPECT_EQ(37, test.outskirt_content_.at(1).at(0).size()) << "SRR490124-4pairs outskirt_content_[1][0].size() not correct for " << context << '\n';
	EXPECT_EQ(36, test.outskirt_content_.at(1).at(1).size()) << "SRR490124-4pairs outskirt_content_[1][1].size() not correct for " << context << '\n';
	EXPECT_EQ(40, test.outskirt_content_.at(1).at(2).size()) << "SRR490124-4pairs outskirt_content_[1][2].size() not correct for " << context << '\n';
	EXPECT_EQ(37, test.outskirt_content_.at(1).at(3).size()) << "SRR490124-4pairs outskirt_content_[1][3].size() not correct for " << context << '\n';
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(1==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4-20 "-" $4-1)}(0==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4+length($10) "-" $4+length($10)+19)}' | seqtk seq -r | awk '(2==NR%4){store=$0}(0==NR%4){print $0 store}' | awk '{for(pos=1;pos<=length($0);pos+=1){print substr($0,pos,1), pos-1}}' | sort | uniq -c | sort -k3,3n -k2,2
	TestOutskirtContent(test, 1, 0, 0, 0, 1, 1, context);
	TestOutskirtContent(test, 1, 1, 0, 1, 1, 0, context);
	TestOutskirtContent(test, 1, 2, 1, 0, 1, 0, context);
	TestOutskirtContent(test, 1, 3, 1, 0, 0, 1, context);
	TestOutskirtContent(test, 1, 4, 0, 0, 0, 2, context);
	TestOutskirtContent(test, 1, 18, 0, 1, 1, 0, context);
	TestOutskirtContent(test, 1, 19, 0, 1, 1, 0, context);
	TestOutskirtContent(test, 1, 20, 1, 1, 0, 0, context);
	TestOutskirtContent(test, 1, 21, 0, 0, 1, 1, context);
	TestOutskirtContent(test, 1, 22, 0, 1, 0, 1, context);
	TestOutskirtContent(test, 1, 35, 0, 0, 1, 1, context);
	TestOutskirtContent(test, 1, 36, 0, 1, 0, 1, context);
	TestOutskirtContent(test, 1, 37, 1, 0, 1, 0, context);
	TestOutskirtContent(test, 1, 38, 2, 0, 0, 0, context);
	TestOutskirtContent(test, 1, 39, 0, 0, 2, 0, context);
}

void FragmentDistributionStatsTest::TestDuplicates(const FragmentDistributionStats &test){
	// samtools view -q 10 -f 3 ecoli-duplicates.bam | awk '(0 != substr($1,12,1) && ($4>$8 || $4==$8 && int($2%256/128))){len=0; for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e || "D"==e){len+=num}; num=0}}; print $4-$8+len}' | sort -n | uniq -c
	EXPECT_EQ(252, test.insert_lengths_.size()) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(3, test.insert_lengths_[50]) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(2, test.insert_lengths_[99]) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(2, test.insert_lengths_[150]) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(4, test.insert_lengths_[151]) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(2, test.insert_lengths_[201]) << "insert_lengths_ wrong in duplicates test\n";
	EXPECT_EQ(1, test.insert_lengths_[301]) << "insert_lengths_ wrong in duplicates test\n";
}

void FragmentDistributionStatsTest::TestCrossDuplicates(const FragmentDistributionStats &test){
	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(1230, test.abundance_.size() ) << "abundance_ wrong in cross duplicates test\n";
	EXPECT_EQ(0, SumVect(test.abundance_)) << "abundance_ wrong in cross duplicates test\n";

	EXPECT_EQ(0, test.insert_lengths_.size() ) << "insert_lengths_ wrong in cross duplicates test\n";
}

void FragmentDistributionStatsTest::TestCoverage(const FragmentDistributionStats &test){
	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(1230, test.abundance_.size() ) << "abundance_ wrong in coverage test\n";
	EXPECT_EQ(1, test.abundance_.at(2) ) << "abundance_ wrong in coverage test\n";
	EXPECT_EQ(1, test.abundance_.at(3) ) << "abundance_ wrong in coverage test\n";
	EXPECT_EQ(1, test.abundance_.at(8) ) << "abundance_ wrong in coverage test\n";
}
// Adapter with mapping quality 8 still used for adapters
void FragmentDistributionStatsTest::TestAdapters(const FragmentDistributionStats &test, const char *context, bool bwa){
	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(177, test.insert_lengths_.size()) << "insert_lengths_.size() wrong with adapters for " << context;
	EXPECT_EQ(2, test.insert_lengths_[0]) << "insert_lengths_[0] wrong with adapters for " << context;
	EXPECT_EQ(1, test.insert_lengths_[58]) << "insert_lengths_[58] wrong with adapters\n for " << context;
	EXPECT_EQ(1, test.insert_lengths_[80]) << "insert_lengths_[80] wrong with low mapping quality adapters for " << context;
	EXPECT_EQ(1, test.insert_lengths_[81]) << "insert_lengths_[81] wrong with adapters for " << context;
	EXPECT_EQ(1, test.insert_lengths_[95]) << "insert_lengths_[95] wrong with adapters for " << context;
	EXPECT_EQ(1, test.insert_lengths_[176]) << "insert_lengths_[176] wrong with adapters for " << context;
	EXPECT_EQ(7, SumVect(test.insert_lengths_)) << "Sum of insert_lengths_ wrong with adapters for " << context;

	// Check if adapter part of read has been ignored for these statistics
	if(bwa){
		// samtools view ecoli-SRR490124-adapter-bwa.bam -q 10 -f 16 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8 "-" $8-$9-1)}' | seqtk seq | awk '(0==NR%2){len=length($0); print int(gsub(/[CG]/,"",$0)*100/len+0.5)}' | sort -n
		EXPECT_EQ(20, test.gc_fragment_content_.size() ) << "for " << context;
		EXPECT_EQ(1, test.gc_fragment_content_[48] ) << "for " << context;
		EXPECT_EQ(1, test.gc_fragment_content_[66] ) << "for " << context;
	}
	else{
		// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 16 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8 "-" $4+length($10)-1)}' | seqtk seq | awk '(0==NR%2){len=length($0); print int(gsub(/[CG]/,"",$0)*100/len+0.5)}' | sort -n
		EXPECT_EQ(7, test.gc_fragment_content_.size() ) << "for " << context;
	}
	EXPECT_EQ(1, test.gc_fragment_content_[47] ) << "for " << context;
	EXPECT_EQ(1, test.gc_fragment_content_[53] ) << "for " << context;

	// cat <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 144 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq) <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq -r) | awk '(0==NR%2){print ">0", substr($0,1,10); print substr($0,length($0)-9,10); print ">1", substr($0,11,10); print substr($0,length($0)-19,10);print ">2", substr($0,21,10); print substr($0,length($0)-29,10);}' | seqtk seq -r | awk '{if(1==NR%2){sur=substr($0, 2, 1); print "start", sur, $2}else{print "end", sur, $0}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}' | sort
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(817646)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(940919)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(10510)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(477197)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(495577)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(718722)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(301332)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(0).at(824719)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(339634)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(1).at(595486)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(158717)) << "for " << context;
	EXPECT_EQ(1, test.fragment_surroundings_.counts_.at(2).at(289540)) << "for " << context;

	if(bwa){
		// samtools view ecoli-SRR490124-adapter-bwa.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-20 "-" $8-$9-1+20)}' | seqtk seq -r | awk '(0==NR%2){print substr($0,1,20) substr($0,length($0)-19,20)}' | awk '{for(pos=1;pos<=length($0);++pos){print substr($0,pos,1), pos-1}}' | sort -k1,1 -k2,2n | uniq -c
		EXPECT_EQ(39, test.outskirt_content_.at(1).at(1).size() ) << "for " << context;
		EXPECT_EQ(2, test.outskirt_content_.at(1).at(1)[31] ) << "for " << context;
		EXPECT_EQ(2, test.outskirt_content_.at(1).at(1)[34] ) << "for " << context;
	}
	else{
		// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-20 "-" $4+length($10)-1+20)}' | seqtk seq -r | awk '(0==NR%2){print substr($0,1,20) substr($0,length($0)-19,20)}' | awk '{for(pos=1;pos<=length($0);++pos){print 1, substr($0,pos,1), pos-1}}' | sort -k1,2 -k3,3n
		EXPECT_EQ(31, test.outskirt_content_.at(1).at(1).size() ) << "for " << context;
		EXPECT_EQ(1, test.outskirt_content_.at(1).at(1)[31] ) << "for " << context;
		EXPECT_EQ(1, test.outskirt_content_.at(1).at(1)[34] ) << "for " << context;
	}

	EXPECT_EQ(1, test.outskirt_content_.at(1).at(1)[4] ) << "for " << context;
	EXPECT_EQ(1, test.outskirt_content_.at(1).at(1)[7] ) << "for " << context;
}

void FragmentDistributionStatsTest::BiasCalculationThread(FragmentDistributionStats &test, const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, mutex &print_mutex){
	test.ExecuteBiasCalculations( reference, duplications, thread_values, print_mutex );
}

void FragmentDistributionStatsTest::TestBiasCalculationVectorsPreprocessing(){
	BiasCalculationVectors test;
	test.sites_.reserve(species_reference_.SequenceLength(0));

	species_reference_.GetFragmentSites(test.sites_, 0, 10, 0, species_reference_.SequenceLength(0));

	test.sites_.at(0).count_forward_ = 1;
	test.sites_.at(40).count_forward_ = 1;
	test.sites_.at(50).count_forward_ = 1;
	test.sites_.at(112).count_forward_ = 1;
	test.sites_.at(192).count_forward_ = 1; // Blocked by N's (should be ignored)
	test.sites_.at(331).count_reverse_ = 1;

	test.GetCounts();
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0==$1 || 40==$1 || 50==$1 || 112==$1 || 331==$1){print gsub(/[GC]/,"",$2)}' | sort | uniq -c
	EXPECT_EQ(2, test.gc_count_.at(0));
	EXPECT_EQ(1, test.gc_count_.at(30));
	EXPECT_EQ(1, test.gc_count_.at(40));
	EXPECT_EQ(1, test.gc_count_.at(50));
	EXPECT_EQ(5, SumVect(test.gc_count_));
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0==$1 || 40==$1 || 50==$1 || 112==$1 || 331==$1){pos=5; print substr($3,pos+1,1); print substr($4,pos+1,1)}' | sort | uniq -c
	EXPECT_EQ(1, test.sur_count_.at(5*4+0));
	EXPECT_EQ(3, test.sur_count_.at(5*4+1));
	EXPECT_EQ(2, test.sur_count_.at(5*4+2));
	EXPECT_EQ(4, test.sur_count_.at(5*4+3));
	EXPECT_EQ(5, test.sur_count_.at(13*4+0));
	EXPECT_EQ(0, test.sur_count_.at(13*4+1));
	EXPECT_EQ(2, test.sur_count_.at(13*4+2));
	EXPECT_EQ(3, test.sur_count_.at(13*4+3));

	test.RemoveUnnecessarySites();
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0==$1 || 40==$1 || 50==$1 || 112==$1 || 331==$1){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10)}'
	// removes also not needed gc, which is not the implementation anymore: cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0 == gsub("N","",$3)){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10), $1}' | awk '((0 == $1 || 3 == $1 || 4 == $1 || 5 == $1) && substr($2,3,1) != "C" && substr($5,3,1) != "C" && substr($2,7,1) != "C" && substr($5,7,1) != "C" && substr($2,7,1) != "G" && substr($5,7,1) != "G" && substr($3,4,1) != "C" && substr($6,4,1) != "C" && substr($3,7,1) != "G" && substr($6,7,1) != "G" && substr($4,2,1) != "C" && substr($7,2,1) != "C" && substr($4,4,1) != "C" && substr($7,4,1) != "C" && substr($4,4,1) != "G" && substr($7,4,1) != "G" && substr($4,8,1) != "G" && substr($7,8,1) != "G")' | wc -l
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0 == gsub("N","",$3)){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10), $1}' | awk '(substr($2,3,1) != "C" && substr($5,3,1) != "C" && substr($2,7,1) != "C" && substr($5,7,1) != "C" && substr($2,7,1) != "G" && substr($5,7,1) != "G" && substr($3,4,1) != "C" && substr($6,4,1) != "C" && substr($3,7,1) != "G" && substr($6,7,1) != "G" && substr($4,2,1) != "C" && substr($7,2,1) != "C" && substr($4,4,1) != "C" && substr($7,4,1) != "C" && substr($4,4,1) != "G" && substr($7,4,1) != "G" && substr($4,8,1) != "G" && substr($7,8,1) != "G")' | wc -l
	EXPECT_EQ(26, test.sites_.size());

	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0 == gsub("N","",$3)){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10), $1}' | awk '((0 == $1 || 3 == $1 || 4 == $1 || 5 == $1) && substr($2,3,1) != "C" && substr($5,3,1) != "C" && substr($2,7,1) != "C" && substr($5,7,1) != "C" && substr($2,7,1) != "G" && substr($5,7,1) != "G" && substr($3,4,1) != "C" && substr($6,4,1) != "C" && substr($3,7,1) != "G" && substr($6,7,1) != "G" && substr($4,2,1) != "C" && substr($7,2,1) != "C" && substr($4,4,1) != "C" && substr($7,4,1) != "C" && substr($4,4,1) != "G" && substr($7,4,1) != "G" && substr($4,8,1) != "G" && substr($7,8,1) != "G"){print $1}' | sort | uniq -c
	EXPECT_EQ(2*2, test.gc_sites_.at(0));
	EXPECT_EQ(3*2, test.gc_sites_.at(30));
	EXPECT_EQ(5*2, test.gc_sites_.at(40));
	EXPECT_EQ(9*2, test.gc_sites_.at(50));
	EXPECT_EQ(26*2, SumVect(test.gc_sites_));
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0 == gsub("N","",$3)){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10), $1}' | awk '(substr($2,3,1) != "C" && substr($5,3,1) != "C" && substr($2,7,1) != "C" && substr($5,7,1) != "C" && substr($2,7,1) != "G" && substr($5,7,1) != "G" && substr($3,4,1) != "C" && substr($6,4,1) != "C" && substr($3,7,1) != "G" && substr($6,7,1) != "G" && substr($4,2,1) != "C" && substr($7,2,1) != "C" && substr($4,4,1) != "C" && substr($7,4,1) != "C" && substr($4,4,1) != "G" && substr($7,4,1) != "G" && substr($4,8,1) != "G" && substr($7,8,1) != "G"){pos=5; print substr($2,pos+1,1); print substr($5,pos+1,1)}' | sort | uniq -c
	EXPECT_EQ(7*2, test.sur_sites_.at(5*4+0));
	EXPECT_EQ(17*2, test.sur_sites_.at(5*4+1));
	EXPECT_EQ(9*2, test.sur_sites_.at(5*4+2));
	EXPECT_EQ(19*2, test.sur_sites_.at(5*4+3));
	// cat <(seqtk seq reference-test-withN.fa | awk '(2==NR)') | awk '{for(i=51;i<=length($0)-50-10; i+=1){print substr($0,i-10,30), i-51}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){forward=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(forward,11,10), forward, $0}' | awk '(0 == gsub("N","",$3)){print gsub(/[GC]/,"",$2), substr($3,1,10), substr($3,11,10), substr($3,21,10), substr($4,1,10), substr($4,11,10), substr($4,21,10), $1}' | awk '(substr($2,3,1) != "C" && substr($5,3,1) != "C" && substr($2,7,1) != "C" && substr($5,7,1) != "C" && substr($2,7,1) != "G" && substr($5,7,1) != "G" && substr($3,4,1) != "C" && substr($6,4,1) != "C" && substr($3,7,1) != "G" && substr($6,7,1) != "G" && substr($4,2,1) != "C" && substr($7,2,1) != "C" && substr($4,4,1) != "C" && substr($7,4,1) != "C" && substr($4,4,1) != "G" && substr($7,4,1) != "G" && substr($4,8,1) != "G" && substr($7,8,1) != "G"){pos=3; print substr($3,pos+1,1); print substr($6,pos+1,1)}' | sort | uniq -c
	EXPECT_EQ(17*2, test.sur_sites_.at(13*4+0));
	EXPECT_EQ(0*2, test.sur_sites_.at(13*4+1));
	EXPECT_EQ(18*2, test.sur_sites_.at(13*4+2));
	EXPECT_EQ(17*2, test.sur_sites_.at(13*4+3));

	vector<double> parameter;
	parameter.resize(test.sur_bias_.size(), 1.0);
	test.DeactivateZeroCounts(parameter, 0.0, 0);

	EXPECT_EQ(1.0, parameter.at(5*4+0));
	EXPECT_EQ(1.0, parameter.at(5*4+1));
	EXPECT_EQ(1.0, parameter.at(5*4+2));
	EXPECT_EQ(1.0, parameter.at(5*4+3));

	EXPECT_EQ(1.0, parameter.at(13*4+0));
	EXPECT_EQ(0.0, parameter.at(13*4+1));
	EXPECT_EQ(1.0, parameter.at(13*4+2));
	EXPECT_EQ(1.0, parameter.at(13*4+3));
}

void FragmentDistributionStatsTest::TestBiasCalculationVectorsNormalizations(){
	BiasCalculationVectors test;

	// NormGC
	test.gc_sites_.fill(0);
	test.gc_bias_.fill(0.0);

	test.gc_sites_.at(12) = 40;
	test.gc_sites_.at(33) = 39;
	test.gc_sites_.at(63) = 20;
	test.gc_sites_.at(50) = 1;
	test.total_sites_ = 100;

	test.gc_bias_.at(12) = 0.8;
	test.gc_bias_.at(33) = 0.2;
	test.gc_bias_.at(63) = 0.2;
	test.gc_bias_.at(50) = 1000.0;

	test.NormGC();

	EXPECT_DOUBLE_EQ(2.0, test.gc_bias_.at(12));
	EXPECT_DOUBLE_EQ(0.5, test.gc_bias_.at(33));
	EXPECT_DOUBLE_EQ(0.5, test.gc_bias_.at(63));
	EXPECT_DOUBLE_EQ(2500.0, test.gc_bias_.at(50));

	test.gc_sites_.at(33) = 20;
	test.gc_sites_.at(77) = 19;
	test.gc_bias_.at(77) = 2500.0;

	test.NormGC();

	EXPECT_DOUBLE_EQ(2.0, test.gc_bias_.at(12));
	EXPECT_DOUBLE_EQ(0.5, test.gc_bias_.at(33));
	EXPECT_DOUBLE_EQ(0.5, test.gc_bias_.at(63));
	EXPECT_DOUBLE_EQ(2500.0, test.gc_bias_.at(77));
	EXPECT_DOUBLE_EQ(2500.0, test.gc_bias_.at(50));

	// NormSurroundings
	vector<double> sur_vect;
	sur_vect.resize( 4*Surrounding::Length()+1 , 200.0 );

	test.sur_count_.fill(20);

	sur_vect.at(0) = -2.0;
	sur_vect.at(1) = -0.2;
	sur_vect.at(2) = -0.3;
	sur_vect.at(3) =  0.1;
	sur_vect.at(4) =  0.2;

	test.sur_count_.at(6) = 0;
	sur_vect.at(5) = -0.2;
	sur_vect.at(6) = -0.3;
	sur_vect.at(7) =  0.1;
	sur_vect.at(8) =  0.2;

	test.sur_count_.at(116) = 0;
	test.sur_count_.at(117) = 0;
	test.sur_count_.at(119) = 0;
	sur_vect.at(117) = -0.2;
	sur_vect.at(118) = -0.3;
	sur_vect.at(119) =  0.1;
	sur_vect.at(120) =  0.2;

	test.NormSurroundings(sur_vect);

	EXPECT_DOUBLE_EQ(-2.15, test.sur_bias_.at(0));
	EXPECT_DOUBLE_EQ(-2.25, test.sur_bias_.at(1));
	EXPECT_DOUBLE_EQ(-1.85, test.sur_bias_.at(2));
	EXPECT_DOUBLE_EQ(-1.75, test.sur_bias_.at(3));

	EXPECT_DOUBLE_EQ(-0.1, test.sur_bias_.at(4));
	EXPECT_DOUBLE_EQ(-0.2, test.sur_bias_.at(5));
	EXPECT_DOUBLE_EQ( 0.0, test.sur_bias_.at(6));
	EXPECT_DOUBLE_EQ( 0.3, test.sur_bias_.at(7));

	EXPECT_DOUBLE_EQ( 0.0, test.sur_bias_.at(116));
	EXPECT_DOUBLE_EQ( 0.0, test.sur_bias_.at(117));
	EXPECT_DOUBLE_EQ( 0.0, test.sur_bias_.at(118));
	EXPECT_DOUBLE_EQ( 0.0, test.sur_bias_.at(119));

	// UnnormSurroundingGradients
	vector<double> grad;
	grad.resize( 4*Surrounding::Length()+1 , 5.0 ); // Something != 0.0, to test if it was properly set to zero

	test.sur_grad_.fill(200.0);

	test.sur_grad_.at(0) = -2.0;
	test.sur_grad_.at(1) = -3.0;
	test.sur_grad_.at(2) =  1.0;
	test.sur_grad_.at(3) =  2.0;

	test.sur_grad_.at(4) = -2.0;
	test.sur_grad_.at(5) = -3.0;
	test.sur_grad_.at(6) =  1.0;
	test.sur_grad_.at(7) =  2.0;

	test.sur_grad_.at(116) = -2.0;
	test.sur_grad_.at(117) = -3.0;
	test.sur_grad_.at(118) =  1.0;
	test.sur_grad_.at(119) =  2.0;

	test.UnnormSurroundingGradients(grad, sur_vect);

	EXPECT_DOUBLE_EQ(-2.0, grad.at(0));
	EXPECT_DOUBLE_EQ(-1.5, grad.at(1));
	EXPECT_DOUBLE_EQ(-2.5, grad.at(2));
	EXPECT_DOUBLE_EQ( 1.5, grad.at(3));
	EXPECT_DOUBLE_EQ( 2.5, grad.at(4));

	EXPECT_DOUBLE_EQ(-1.0, grad.at(5));
	EXPECT_DOUBLE_EQ(-2.0, grad.at(6));
	EXPECT_DOUBLE_EQ( 0.0, grad.at(7));
	EXPECT_DOUBLE_EQ( 3.0, grad.at(8));

	EXPECT_DOUBLE_EQ( 0.0, grad.at(117));
	EXPECT_DOUBLE_EQ( 0.0, grad.at(118));
	EXPECT_DOUBLE_EQ( 0.0, grad.at(119));
	EXPECT_DOUBLE_EQ( 0.0, grad.at(120));
}

void FragmentDistributionStatsTest::TestBiasCalculationVectorsSpline(){
	BiasCalculationVectors test;

	test.gc_sites_.fill(0);
	test.gc_count_.fill(0);

	test.gc_sites_.at(98) = 1;
	test.gc_count_.at(98) = 1;
	test.total_sites_ = 1;

	test.DefineStartingKnots();

	EXPECT_EQ(95, test.gc_knots_.at(0));
	EXPECT_EQ(96, test.gc_knots_.at(1));
	EXPECT_EQ(97, test.gc_knots_.at(2));
	EXPECT_EQ(98, test.gc_knots_.at(3));
	EXPECT_EQ(99, test.gc_knots_.at(4));
	EXPECT_EQ(100, test.gc_knots_.at(5));

	test.gc_count_.at(98) = 0;
	test.gc_sites_.at(98) = 0;
	test.gc_sites_.at(12) = 40;
	test.gc_count_.at(12) = 1;
	test.gc_sites_.at(23) = 60;
	test.gc_count_.at(23) = 1;
	test.total_sites_ = 100;

	test.DefineStartingKnots();

	EXPECT_EQ(12, test.gc_knots_.at(0));
	EXPECT_EQ(13, test.gc_knots_.at(1));
	EXPECT_EQ(14, test.gc_knots_.at(2));
	EXPECT_EQ(15, test.gc_knots_.at(3));
	EXPECT_EQ(16, test.gc_knots_.at(4));
	EXPECT_EQ(23, test.gc_knots_.at(5));

	test.gc_sites_.at(98) = 1;
	test.total_sites_ = 101;
	test.DefineStartingKnots();

	EXPECT_EQ(12, test.gc_knots_.at(0));
	EXPECT_EQ(13, test.gc_knots_.at(1));
	EXPECT_EQ(14, test.gc_knots_.at(2));
	EXPECT_EQ(15, test.gc_knots_.at(3));
	EXPECT_EQ(16, test.gc_knots_.at(4));
	EXPECT_EQ(23, test.gc_knots_.at(5));

	test.gc_count_.at(98) = 1;
	test.DefineStartingKnots();

	EXPECT_EQ(12, test.gc_knots_.at(0));
	EXPECT_EQ(13, test.gc_knots_.at(1));
	EXPECT_EQ(14, test.gc_knots_.at(2));
	EXPECT_EQ(15, test.gc_knots_.at(3));
	EXPECT_EQ(23, test.gc_knots_.at(4));
	EXPECT_EQ(98, test.gc_knots_.at(5));

	for(auto gc=13; gc<23; ++gc){
		test.gc_sites_.at(gc) = 10;
		test.gc_count_.at(gc) = 1;
	}
	for(auto gc=24; gc<34; ++gc){
		test.gc_sites_.at(gc) = 20;
		test.gc_count_.at(gc) = 1;
	}
	test.total_sites_ = 401;
	test.DefineStartingKnots();

	EXPECT_EQ(12, test.gc_knots_.at(0));
	EXPECT_EQ(16, test.gc_knots_.at(1));
	EXPECT_EQ(23, test.gc_knots_.at(2));
	EXPECT_EQ(26, test.gc_knots_.at(3));
	EXPECT_EQ(30, test.gc_knots_.at(4));
	EXPECT_EQ(98, test.gc_knots_.at(5));

	test.gc_count_.at(98) = 0;
	test.DefineStartingKnots();

	EXPECT_EQ(12, test.gc_knots_.at(0));
	EXPECT_EQ(16, test.gc_knots_.at(1));
	EXPECT_EQ(23, test.gc_knots_.at(2));
	EXPECT_EQ(26, test.gc_knots_.at(3));
	EXPECT_EQ(30, test.gc_knots_.at(4));
	EXPECT_EQ(33, test.gc_knots_.at(5));

	//void PrepareSplines();
	//void GetSplineCoefficients(double &a, double &b, double &c, double &d, uintPercent k, const std::vector<double> &spline_pars);
	//void CalculateSpline(const std::vector<double> &spline_pars);
	//void CalculateSplineGrad(const std::vector<double> &spline_pars, std::vector<double> &grad, uintNumFits grad_offset);
	//void GetGCSpline();
}

void FragmentDistributionStatsTest::TestBiasCalculationVectorsLikelihoods(){
	//void GetLogLikeBase();
	//static double LogLikelihoodPoisson(const std::vector<double> &x, std::vector<double> &grad, void* f_data);

	//static double LogLikeGcSpline(const std::vector<double> &x, std::vector<double> &grad, void* f_data);

	//static double LogLikelihoodNbinom(const std::vector<double> &x, std::vector<double> &grad, void* f_data);

	//static double LogLikelihoodConstDispersion(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
}

void FragmentDistributionStatsTest::TestBiasCalculation(){
	// To set how much space is needed
	test_->insert_lengths_[250] = 0;
	test_->abundance_.resize(species_reference_.NumberSequences(), 0);

	// Create data by the model that is used for the bias fitting (so it should be easy)
	test_->fragment_sites_by_ref_seq_bin_.resize(species_reference_.NumberSequences());
	test_->fragment_sites_by_ref_seq_bin_cur_id_.resize(species_reference_.NumberSequences());

	uintSeqLen last_size = 0;
	for(uintSeqLen fragment_length=65; fragment_length <= 95; fragment_length += 5){
		CreateCoverageData(fragment_length);

		test_->insert_lengths_.at(fragment_length) += test_->fragment_sites_by_ref_seq_bin_.at(0).size() - last_size;
		last_size = test_->fragment_sites_by_ref_seq_bin_.at(0).size();
	}
	test_->abundance_.at(0) += test_->fragment_sites_by_ref_seq_bin_.at(0).size();

	// Fit data
	FragmentDuplicationStats duplications;
	duplications.PrepareTmpDuplicationVector(100);
	uintNumThreads num_threads = min(static_cast<uintNumThreads>(4), test_with_num_threads);
	array<FragmentDistributionStats::ThreadData, 4> thread_data = {FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0))};
	mutex print_mutex;

	ReduceVerbosity(1); // Suppress warnings
	test_->AddNewBiasCalculations(1, thread_data.at(0), print_mutex);

	thread threads[num_threads];
	for(auto i = num_threads; i--; ){
		threads[i] = thread(BiasCalculationThread, std::ref(*test_), std::cref(species_reference_), std::ref(duplications), std::ref(thread_data.at(i).bias_calc_vects_), std::ref(print_mutex));
	}
	for(auto i = num_threads; i--; ){
		threads[i].join();
	}

	test_->FinalizeBiasCalculation(species_reference_, num_threads, duplications);
	RestoreTestVerbosity();

	EXPECT_NEAR(0.6, test_->gc_fragment_content_bias_.at(27), 0.1);
	EXPECT_NEAR(0.7, test_->gc_fragment_content_bias_.at(28), 0.1);
	EXPECT_NEAR(0.8, test_->gc_fragment_content_bias_.at(29), 0.1);
	EXPECT_NEAR(0.9, test_->gc_fragment_content_bias_.at(30), 0.1);
	EXPECT_NEAR(1.0, test_->gc_fragment_content_bias_.at(31), 0.1);
	EXPECT_NEAR(0.8, test_->gc_fragment_content_bias_.at(32), 0.2);
	EXPECT_NEAR(0.6, test_->gc_fragment_content_bias_.at(33), 0.2);
	EXPECT_NEAR(0.4, test_->gc_fragment_content_bias_.at(34), 0.2);
	EXPECT_DOUBLE_EQ(1.0, *max_element(test_->gc_fragment_content_bias_.begin(), test_->gc_fragment_content_bias_.end()));

	intSurrounding pos_base = 262144;
	EXPECT_NEAR(-0.4, accumulate(test_->fragment_surroundings_bias_.bias_.at(1).begin(), test_->fragment_surroundings_bias_.bias_.at(1).begin()+pos_base, 0.0)/pos_base, 0.2);
	EXPECT_NEAR(-0.4, accumulate(test_->fragment_surroundings_bias_.bias_.at(1).begin()+pos_base, test_->fragment_surroundings_bias_.bias_.at(1).begin()+2*pos_base, 0.0)/pos_base, 0.2);
	EXPECT_NEAR(0.4, accumulate(test_->fragment_surroundings_bias_.bias_.at(1).begin()+2*pos_base, test_->fragment_surroundings_bias_.bias_.at(1).begin()+3*pos_base, 0.0)/pos_base, 0.2);
	EXPECT_NEAR(0.4, accumulate(test_->fragment_surroundings_bias_.bias_.at(1).begin()+3*pos_base, test_->fragment_surroundings_bias_.bias_.at(1).begin()+4*pos_base, 0.0)/pos_base, 0.2);

	for(uintSeqLen fragment_length=65; fragment_length <= 95; fragment_length += 5){
		EXPECT_NEAR(1.0, test_->insert_lengths_bias_.at(fragment_length), 0.1) << "Fragment length " << fragment_length << " wrong\n";
	}
	EXPECT_DOUBLE_EQ(1.0, *max_element(test_->insert_lengths_bias_.begin(), test_->insert_lengths_bias_.end()));

	EXPECT_DOUBLE_EQ(1.0, test_->ref_seq_bias_.at(0));

	EXPECT_NEAR(0.5, test_->dispersion_parameters_.at(0), 0.4);
	EXPECT_NEAR(1.0, test_->dispersion_parameters_.at(1), 0.1);
}

void FragmentDistributionStatsTest::TestUniformBias(){
	test_->abundance_.clear();
	test_->abundance_.resize(species_reference_.NumberSequences(), 0);

	test_->insert_lengths_.Clear();
	test_->insert_lengths_[50] = 500;
	test_->insert_lengths_[51] = 1000;
	test_->insert_lengths_[52] = 750;

	test_->SetUniformBias();

	EXPECT_EQ(test_->abundance_.size(), test_->ref_seq_bias_.size());
	for(auto bias : test_->ref_seq_bias_){
		EXPECT_EQ(1.0, bias);
	}

	EXPECT_EQ(0, test_->gc_fragment_content_bias_.from());
	EXPECT_EQ(101, test_->gc_fragment_content_bias_.to());
	EXPECT_DOUBLE_EQ(1.0, test_->gc_fragment_content_bias_.Max());
	EXPECT_DOUBLE_EQ(101.0, SumVectD(test_->gc_fragment_content_bias_));

	EXPECT_EQ(50, test_->insert_lengths_bias_.from());
	EXPECT_EQ(53, test_->insert_lengths_bias_.to());
	EXPECT_DOUBLE_EQ(0.5, test_->insert_lengths_bias_[50]);
	EXPECT_DOUBLE_EQ(1.0, test_->insert_lengths_bias_[51]);
	EXPECT_DOUBLE_EQ(0.75, test_->insert_lengths_bias_[52]);

	EXPECT_DOUBLE_EQ(Surrounding::Size(), accumulate(test_->fragment_surroundings_bias_.bias_.at(0).begin(), test_->fragment_surroundings_bias_.bias_.at(0).end(), 0.0));
	EXPECT_DOUBLE_EQ(Surrounding::Size(), accumulate(test_->fragment_surroundings_bias_.bias_.at(1).begin(), test_->fragment_surroundings_bias_.bias_.at(1).end(), 0.0));
	EXPECT_DOUBLE_EQ(Surrounding::Size(), accumulate(test_->fragment_surroundings_bias_.bias_.at(2).begin(), test_->fragment_surroundings_bias_.bias_.at(2).end(), 0.0));
	EXPECT_DOUBLE_EQ(1.0, *max_element(test_->fragment_surroundings_bias_.bias_.at(0).begin(), test_->fragment_surroundings_bias_.bias_.at(0).end()));
	EXPECT_DOUBLE_EQ(1.0, *max_element(test_->fragment_surroundings_bias_.bias_.at(1).begin(), test_->fragment_surroundings_bias_.bias_.at(1).end()));
	EXPECT_DOUBLE_EQ(1.0, *max_element(test_->fragment_surroundings_bias_.bias_.at(2).begin(), test_->fragment_surroundings_bias_.bias_.at(2).end()));
}

void FragmentDistributionStatsTest::TestDrawCounts(){
	// Set biases to 0.5 to verify they are used and only use gc bias to really vary bias
	test_->ref_seq_bias_.clear();
	test_->ref_seq_bias_.resize(1, 0.5);

	test_->insert_lengths_bias_.Clear();
	test_->insert_lengths_bias_[367] = 0.5;
	test_->insert_lengths_bias_.Shrink();

	// Surrounding has to be tuned to result in 0.5 for the total bias each (start and end)
	Surrounding start_sur, end_sur;
	start_sur.sur_.at(0) = 873425;
	start_sur.sur_.at(1) = 34;
	start_sur.sur_.at(2) = 7467;
	end_sur.sur_.at(0) = 364;
	end_sur.sur_.at(1) = 856687;
	end_sur.sur_.at(2) = 34562;

	for( auto &block : test_->fragment_surroundings_bias_.bias_){
		block.clear();
		block.resize(Surrounding::Size(), -1000.0); // Set to big negative number, so that bias will be very close to zero
	}
	test_->fragment_surroundings_bias_.bias_.at(0).at(start_sur.sur_.at(0)) = -0.3;
	test_->fragment_surroundings_bias_.bias_.at(1).at(start_sur.sur_.at(1)) = -0.7;
	test_->fragment_surroundings_bias_.bias_.at(2).at(start_sur.sur_.at(2)) = -log(3)+1.0;
	test_->fragment_surroundings_bias_.bias_.at(0).at(end_sur.sur_.at(0)) = -log(3)+1.0;
	test_->fragment_surroundings_bias_.bias_.at(1).at(end_sur.sur_.at(1)) = -0.7;
	test_->fragment_surroundings_bias_.bias_.at(2).at(end_sur.sur_.at(2)) = -0.3;

	test_->gc_fragment_content_bias_.Clear();
	test_->gc_fragment_content_bias_[43] = 1.0; // Use set it, no meaning so far
	test_->gc_fragment_content_bias_.Shrink();

	// Start with poisson
	test_->dispersion_parameters_.at(0) = 0.0;
	test_->dispersion_parameters_.at(1) = 1e-100;

	CheckDrawnCounts(1.0, array<double, 6>({0.3678794, 0.7357589, 0.9196986, 0.9810118, 0.9963402, 0.9994058}), start_sur, end_sur); // R: ppois(0:5, 1.0)
	CheckDrawnCounts(0.5, array<double, 6>({0.6065307, 0.9097960, 0.9856123, 0.9982484, 0.9998279, 0.9999858}), start_sur, end_sur); // R: ppois(0:5, 0.5)
	CheckDrawnCounts(0.1, array<double, 4>({0.9048374, 0.9953212, 0.9998453, 0.9999962}), start_sur, end_sur); // R: ppois(0:5, 0.1)

	// Constant dispersion
	test_->dispersion_parameters_.at(1) = 5.0;
	CheckDrawnCounts(1.0, array<double, 6>({0.6988271, 0.8152983, 0.8735339, 0.9091223, 0.9328479, 0.9494559}), start_sur, end_sur); // R: pnbinom(0:5, size=0.2, mu=1.0)
	CheckDrawnCounts(0.5, array<double, 6>({0.7783705, 0.8895663, 0.9372217, 0.9621840, 0.9764482, 0.9850067}), start_sur, end_sur); // R: pnbinom(0:5, size=0.2, mu=0.5)
	CheckDrawnCounts(0.1, array<double, 6>({0.9221079, 0.9835818, 0.9958765, 0.9988819, 0.9996834, 0.9999078}), start_sur, end_sur); // R: pnbinom(0:5, size=0.2, mu=0.1)

	// Bias dependent dispersion
	test_->dispersion_parameters_.at(0) = 5000.0;
	test_->dispersion_parameters_.at(1) = 10000.0;

	CheckDrawnCounts(1.0, array<double, 6>({0.9993591, 0.9994258, 0.9994591, 0.9994813, 0.9994979, 0.9995113}), start_sur, end_sur); // R: pnbinom(0:5, size=1.0/(5000+1.0*10000), mu=1.0)
	CheckDrawnCounts(0.5, array<double, 6>({0.9995396, 0.9995896, 0.9996145, 0.9996312, 0.9996437, 0.9996537}), start_sur, end_sur); // R: pnbinom(0:5, size=0.5/(5000+0.5*10000), mu=0.5)
	CheckDrawnCounts(0.1, array<double, 6>({0.9998550, 0.9998717, 0.9998800, 0.9998856, 0.9998897, 0.9998931}), start_sur, end_sur); // R: pnbinom(0:5, size=0.1/(5000+0.1*10000), mu=0.1)
}

void FragmentDistributionStatsTest::TestRefSeqSplitting(){
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0)' | wc -l
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0 && length($0) > 40000)' | wc -l
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0 && length($0) > 80000)' | wc -l
	EXPECT_EQ(1230+17+3 , test_->CreateRefBins(species_reference_, 40000)); // Maximum sequence length is 88768, so 40000 to have sequences with more than one cut)

	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0 && length($0) > 40000){print NR/2-1, length($0)}'
	EXPECT_EQ(0, test_->ref_seq_start_bin_.at(0));
	EXPECT_EQ(1+1, test_->ref_seq_start_bin_.at(1));
	EXPECT_EQ(1198+5, test_->ref_seq_start_bin_.at(1198));
	EXPECT_EQ(1200+6, test_->ref_seq_start_bin_.at(1200));
	EXPECT_EQ(1206+8, test_->ref_seq_start_bin_.at(1206));

	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0){print length($0)}' | head -n5
	EXPECT_EQ(33220, test_->RefSeqSplitLength(0, species_reference_) );
	EXPECT_EQ(33316, test_->RefSeqSplitLength(1, species_reference_) );
	EXPECT_EQ(29122, test_->RefSeqSplitLength(1200, species_reference_) );
	EXPECT_EQ(29589, test_->RefSeqSplitLength(1225, species_reference_) );

	EXPECT_EQ(0, test_->GetRefSeqId(0));
	EXPECT_EQ(0, test_->GetRefSeqId(1));
	EXPECT_EQ(1, test_->GetRefSeqId(2));
	EXPECT_EQ(1198, test_->GetRefSeqId(1203));
	EXPECT_EQ(1198, test_->GetRefSeqId(1204));
	EXPECT_EQ(1200, test_->GetRefSeqId(1206));
	EXPECT_EQ(1200, test_->GetRefSeqId(1207));
	EXPECT_EQ(1200, test_->GetRefSeqId(1208));
	EXPECT_EQ(1201, test_->GetRefSeqId(1209));

	EXPECT_EQ(0, test_->GetRefSeqBin(0, 0, species_reference_));
	EXPECT_EQ(0, test_->GetRefSeqBin(0, 33219, species_reference_));
	EXPECT_EQ(1, test_->GetRefSeqBin(0, 33220, species_reference_));
	EXPECT_EQ(1, test_->GetRefSeqBin(0, 66438, species_reference_));
	EXPECT_EQ(2, test_->GetRefSeqBin(1, 0, species_reference_));
	EXPECT_EQ(2, test_->GetRefSeqBin(1, 33315, species_reference_));
	EXPECT_EQ(1206, test_->GetRefSeqBin(1200, 0, species_reference_));
	EXPECT_EQ(1206, test_->GetRefSeqBin(1200, 29121, species_reference_));
	EXPECT_EQ(1207, test_->GetRefSeqBin(1200, 29122, species_reference_));
	EXPECT_EQ(1207, test_->GetRefSeqBin(1200, 58243, species_reference_));
	EXPECT_EQ(1208, test_->GetRefSeqBin(1200, 58244, species_reference_));
	EXPECT_EQ(1208, test_->GetRefSeqBin(1200, 87364, species_reference_));
	EXPECT_EQ(1240, test_->GetRefSeqBin(1225, 0, species_reference_));
	EXPECT_EQ(1240, test_->GetRefSeqBin(1225, 29588, species_reference_));
	EXPECT_EQ(1241, test_->GetRefSeqBin(1225, 29589, species_reference_));
	EXPECT_EQ(1241, test_->GetRefSeqBin(1225, 59177, species_reference_));
	EXPECT_EQ(1242, test_->GetRefSeqBin(1225, 59178, species_reference_));
	EXPECT_EQ(1242, test_->GetRefSeqBin(1225, 88767, species_reference_));

	vector<bool> used_ref_seqs;
	used_ref_seqs.resize(species_reference_.NumberSequences(), false);
	used_ref_seqs.at(0) = true;
	used_ref_seqs.at(1) = true;
	used_ref_seqs.at(1198) = true;
	used_ref_seqs.at(1200) = true;
	used_ref_seqs.at(1206) = true;
	EXPECT_EQ(10 , species_reference_.NumRefSeqBinsInNxx(used_ref_seqs, 40000) );

	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0){print length($0)/(1+int(length($0)/40000))}' | sort -nr | head -n5
	EXPECT_EQ( 39041, test_->MaxRefSeqBinLength(species_reference_) );
}

void FragmentDistributionStatsTest::TestRefBinProcessing(){
	// Test pipeline for fragments to fit
	test_->tmp_insert_lengths_.resize(1001);
	test_->ref_seq_start_bin_.clear();
	auto num_bins = test_->CreateRefBins(species_reference_, 40000);
	test_->fragment_sites_by_ref_seq_bin_.resize( num_bins );
	std::array<uintRefSeqBin, 9> used_bins({0, 1, 2, 1206, 1207, 1208, 1240, 1241, 1242});
	for( auto ref_bin : used_bins ){
		test_->fragment_sites_by_ref_seq_bin_.at(ref_bin).resize(2);
	}
	test_->fragment_sites_by_ref_seq_bin_cur_id_.resize( num_bins );
	test_->fragment_sites_by_ref_seq_bin_by_insert_length_.resize(num_bins);

	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0 && length($0) > 40000){print NR/2-1, length($0)}'
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk '(NR%2==0){print length($0)}' | head -n5
	std::array<uintRefSeqBin, 9> ref_seq_ids = {0, 0, 1, 1200, 1200, 1200, 1225, 1225, 1225};
	std::array<uintRefSeqBin, 18> frag_length = {250, 500, 750, 1000, 250, 500, 750, 1000, 250, 500, 750, 1000, 250, 500, 750, 1000, 250, 500};
	std::array<uintRefSeqBin, 9> start = {0, 33220, 0, 0, 29122, 58244, 0, 29589, 59178};
	std::array<uintRefSeqBin, 9> end = {33220, 66439, 33316, 29122, 58244, 87365, 29589, 59178, 88768};
	std::array<uintRefSeqBin, 18> positions = {50, 33219, 33220, 66439-50-1000, 50, 33316-50-500, 50, 29121, 29122, 58243, 58244, 87365-50-1000, 50, 29588, 29589, 59177, 59178, 88768-50-500}; // Remove kMinDistToRefSeqEnds and Fragment length from last position in last bin of reference sequence
	std::array<uintRefSeqBin, 18> orientation = {0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0};

	for( auto i=ref_seq_ids.size() ; i--; ){
		for( auto k=2; k--; ){
			test_->AddFragmentSite(ref_seq_ids.at(i), frag_length.at(2*i+k), positions.at(2*i+k), orientation.at(2*i+k), species_reference_);
		}
	}

	for( auto i=ref_seq_ids.size() ; i--; ){
		for( auto k=2; k--; ){
			EXPECT_EQ( 2*positions.at(2*i+k)+orientation.at(2*i+k), test_->fragment_sites_by_ref_seq_bin_.at(used_bins.at(i)).at(!k).first ) << "Bin " << used_bins.at(i) << " k " << k;
			EXPECT_EQ( frag_length.at(2*i+k), test_->fragment_sites_by_ref_seq_bin_.at(used_bins.at(i)).at(!k).second ) << "Bin " << used_bins.at(i) << " k " << k;
		}
	}

	vector<uintSeqLen> tmp_storage;
	vector<FragmentSite> sites;
	for( auto i=ref_seq_ids.size() ; i--; ){
		test_->SortFragmentSites(used_bins.at(i), tmp_storage);

		for( auto k=2; k--; ){
			species_reference_.GetFragmentSites( sites, ref_seq_ids.at(i), frag_length.at(2*i+k), start.at(i), end.at(i) );
			test_->AddFragmentsToSites(sites, test_->fragment_sites_by_ref_seq_bin_by_insert_length_.at(used_bins.at(i)).at(frag_length.at(2*i+k)), max(Reference::kMinDistToRefSeqEnds, start.at(i)) );

			if(0==k){
				EXPECT_EQ( !orientation.at(2*i+k), sites.at(0).count_forward_ ) << "Bin " << used_bins.at(i) << " k " << k;
				EXPECT_EQ( orientation.at(2*i+k), sites.at(0).count_reverse_ ) << "Bin " << used_bins.at(i) << " k " << k;
			}
			else{
				auto pos = positions.at(2*i+k)-max(Reference::kMinDistToRefSeqEnds, start.at(i));
				EXPECT_EQ( !orientation.at(2*i+k), sites.at(pos).count_forward_ ) << "Bin " << used_bins.at(i) << " k " << k;
				EXPECT_EQ( orientation.at(2*i+k), sites.at(pos).count_reverse_ ) << "Bin " << used_bins.at(i) << " k " << k;
			}
		}
	}
}

namespace reseq{
	TEST_F(FragmentDistributionStatsTest, BiasCalculationVectors){
		LoadReference("reference-test-withN.fa");

		TestBiasCalculationVectorsPreprocessing();
		TestBiasCalculationVectorsNormalizations();
		TestBiasCalculationVectorsSpline();
		TestBiasCalculationVectorsLikelihoods();
	}

	TEST_F(FragmentDistributionStatsTest, BiasCalculation){
		// samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:100-199 | awk '(NR==1){print $0 "_1000times"}(NR==2){for(i=0;i<1000;i++){print $0}}' | seqtk seq > reference-bias-calculation.fa
		LoadReference("reference-bias-calculation.fa");
		CreateTestObject(&species_reference_);

		TestBiasCalculation();
		TestUniformBias();
	}

	TEST_F(FragmentDistributionStatsTest, UpdateRefSeqBias){
		LoadReference("reference-test.fa");
		CreateTestObject(&species_reference_);
		std::mt19937_64 rgen;

		ReduceVerbosity(1); // Suppress warnings

		test_->ref_seq_bias_.resize(1, 0.5);
		EXPECT_TRUE(test_-> UpdateRefSeqBias(RefSeqBiasSimulation::kKeep, std::string(), species_reference_, rgen));
		EXPECT_EQ(2, test_->ref_seq_bias_.size()) << "Keeping reference bias does not switch to no bias correctly.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(0)) << "Keeping reference bias does not switch to no bias correctly.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(1)) << "Keeping reference bias does not switch to no bias correctly.";

		RestoreTestVerbosity();

		test_->ref_seq_bias_.at(0) = 0.25;
		test_->ref_seq_bias_.at(1) = 0.5;
		EXPECT_TRUE(test_-> UpdateRefSeqBias(RefSeqBiasSimulation::kKeep, std::string(), species_reference_, rgen));
		EXPECT_EQ(2, test_->ref_seq_bias_.size()) << "Keeping reference bias does not work.";
		EXPECT_EQ(0.25, test_->ref_seq_bias_.at(0)) << "Keeping reference bias does not work.";
		EXPECT_EQ(0.5, test_->ref_seq_bias_.at(1)) << "Keeping reference bias does not work.";

		test_->ref_seq_bias_.resize(1);
		EXPECT_TRUE(test_-> UpdateRefSeqBias(RefSeqBiasSimulation::kNo, std::string(), species_reference_, rgen));
		EXPECT_EQ(2, test_->ref_seq_bias_.size()) << "No reference bias does not work.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(0)) << "No reference bias does not work.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(1)) << "No reference bias does not work.";

		test_->ref_seq_bias_.resize(1);
		EXPECT_TRUE(test_-> UpdateRefSeqBias(RefSeqBiasSimulation::kFile, std::string(PROJECT_SOURCE_DIR)+"/test/ref-bias-test.txt", species_reference_, rgen));
		EXPECT_EQ(2, test_->ref_seq_bias_.size()) << "Reference bias from file does not work.";
		EXPECT_EQ(2.0, test_->ref_seq_bias_.at(0)) << "Reference bias from file does not work.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(1)) << "Reference bias from file does not work.";
	}

	TEST_F(FragmentDistributionStatsTest, Functionality){
		LoadReference("drosophila-GCF_000001215.4_cut.fna"); // So we have many sequences
		ASSERT_TRUE( test_ = new FragmentDistributionStats ) << "Could not allocate memory for FragmentDistributionStats object\n";

		TestDrawCounts();
		TestRefSeqSplitting();
		TestRefBinProcessing();
	}
}
