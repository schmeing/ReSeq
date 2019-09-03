#include "FragmentDistributionStatsTest.h"
using reseq::FragmentDistributionStatsTest;

#include <algorithm>
using std::max_element;
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

#include "utilities.h"
using reseq::utilities::IntPow;
using reseq::utilities::InvLogit2;
using reseq::utilities::Percent;

void FragmentDistributionStatsTest::CreateTestObject(const Reference *ref){
	ASSERT_TRUE( test_ = new FragmentDistributionStats ) << "Could not allocate memory for FragmentDistributionStats object\n";

	test_->CreateRefBins(*ref, 100000000);

	vector<uint64_t> reads_per_ref_seq_;
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
		uint16_t template_segment,
		uint32_t at_pos,
		uint64_t cont_a,
		uint64_t cont_c,
		uint64_t cont_g,
		uint64_t cont_t,
		const char * context ){
	EXPECT_EQ(cont_a, test.outskirt_content_[template_segment][0][at_pos]) << "outskirt_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_c, test.outskirt_content_[template_segment][1][at_pos]) << "outskirt_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_g, test.outskirt_content_[template_segment][2][at_pos]) << "outskirt_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_t, test.outskirt_content_[template_segment][3][at_pos]) << "outskirt_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
}

void FragmentDistributionStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void FragmentDistributionStatsTest::CreateCoverageDataHelper(uint8_t gc_perc, const array<int32_t, 3> &start_sur, const array<int32_t, 3> &end_sur, uint32_t start_pos, uint32_t fragment_length, mt19937_64 &rgen){
	uint16_t ref_seq_id(0);

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
	switch( static_cast<char>(species_reference_.ReferenceSequence(0)[start_pos]) ){
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
	switch( static_cast<char>(species_reference_.ReferenceSequence(0)[start_pos+fragment_length-1]) ){
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
	for(int32_t entries=count; entries--; ){
		test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).push_back({start_pos << 1, fragment_length}); // Forward
		++test_->fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_id);
	}

	count = FragmentDistributionStats::NegativeBinomial(probability, dispersion, zero_to_one(rgen));
	for(int32_t entries=count; entries--; ){
		test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).push_back({(start_pos << 1)+1, fragment_length}); // Reverse
		++test_->fragment_sites_by_ref_seq_bin_cur_id_.at(ref_seq_id);
	}
}

void FragmentDistributionStatsTest::CreateCoverageData(uint32_t fragment_length){
	mt19937_64 rgen;
	rgen.seed( 201907171113 );

	uint16_t ref_seq_id(0);
	//test_->fragment_sites_by_ref_seq_bin_by_insert_length_.at(ref_seq_id).at(fragment_length).reserve(100*2*species_reference_.SequenceLength(ref_seq_id)); // 100: Max entries per site, 2: forward/reverse
	test_->fragment_sites_by_ref_seq_bin_.at(ref_seq_id).reserve(4*species_reference_.SequenceLength(ref_seq_id));
	uint16_t dist_ref_seq_ends(50);

	uint8_t gc_perc;
	uint32_t n_count(0);
	uint64_t gc = species_reference_.GCContentAbsolut(n_count, ref_seq_id, dist_ref_seq_ends, dist_ref_seq_ends+fragment_length);
	if(n_count < fragment_length){
		gc_perc = Percent(gc, fragment_length-n_count);
	}
	else{
		gc_perc = 50;
	}

	array<int32_t, Reference::num_surrounding_blocks_> start_sur, end_sur;
	species_reference_.ForwardSurroundingWithN(start_sur, ref_seq_id, dist_ref_seq_ends);
	species_reference_.ReverseSurroundingWithN(end_sur, ref_seq_id, dist_ref_seq_ends+fragment_length-1);

	CreateCoverageDataHelper(gc_perc, start_sur, end_sur, dist_ref_seq_ends, fragment_length, rgen);

	const Dna5String &ref_seq(species_reference_[ref_seq_id]);
	for( uint32_t start_pos=dist_ref_seq_ends; start_pos < length(ref_seq)-fragment_length-dist_ref_seq_ends; ){
		species_reference_.UpdateGC( gc, n_count, ref_seq_id, start_pos, start_pos+fragment_length );
		if(n_count < fragment_length){
			gc_perc = Percent(gc, fragment_length-n_count);
		}
		else{
			gc_perc = 50;
		}

		species_reference_.UpdateReverseSurroundingWithN( end_sur, ref_seq, start_pos+fragment_length, ref_seq_id );
		species_reference_.UpdateForwardSurroundingWithN( start_sur, ref_seq, ++start_pos, ref_seq_id );

		CreateCoverageDataHelper(gc_perc, start_sur, end_sur, start_pos, fragment_length, rgen);
	}
}

void FragmentDistributionStatsTest::TestSrr490124Equality(const FragmentDistributionStats &test, const char *context){
	EXPECT_EQ(1, test.abundance_.size() ) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.abundance_[0] ) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{if(1==NR%2){store=$4}else{print store, $4, $4-store+length($10)}}'
	EXPECT_EQ(59, test.insert_lengths_.size()) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.insert_lengths_[126]) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.insert_lengths_[184]) << "SRR490124-4pairs insert_lengths_ wrong for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{if(1==NR%2){store=$4}else{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" store "-" $4+length($10)-1)}}' | seqtk seq | awk '(0==NR%2){len=length($1); print (gsub("G","",$1)+gsub("C","",$1))/len}'
	EXPECT_EQ(2, test.gc_fragment_content_.size()) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_fragment_content_[52]) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_fragment_content_[53]) << "SRR490124-4pairs gc_fragment_content_ wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.bam -q 10 -f 144 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq) <(samtools view ecoli-SRR490124-4pairs.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq -r) | awk '(0==NR%2){print ">0", substr($0,1,10); print substr($0,length($0)-9,10); print ">1", substr($0,11,10); print substr($0,length($0)-19,10);print ">2", substr($0,21,10); print substr($0,length($0)-29,10);}' | seqtk seq -r | awk '{if(1==NR%2){sur=substr($0, 2, 1); print "start", sur, $2}else{print "end", sur, $0}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}' | sort
	EXPECT_EQ(1, test.fragment_surroundings_[0][39619]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[0][1009062]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[1][561168]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[1][996771]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[2][307041]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[2][856654]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[0][607462]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[0][983881]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[1][354618]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[1][765926]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[2][333453]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.fragment_surroundings_[2][405510]) << "SRR490124-4pairs fragment_surroundings_ not correct for " << context << '\n';

	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(1==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4-20 "-" $4-1)}(0==NR%2){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4+length($10) "-" $4+length($10)+19)}' | seqtk seq -r | awk '(2==NR%4){store=$0}(0==NR%4){print $0 store}' | awk '{for(pos=1;pos<=length($0);pos+=1){print substr($0,pos,1), pos-1}}' | sort | uniq -c | sort -k2,2 -k3,3n
	EXPECT_EQ(0, test.outskirt_content_[0][0].size()) << "SRR490124-4pairs outskirt_content_[0][0].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_[0][1].size()) << "SRR490124-4pairs outskirt_content_[0][1].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_[0][2].size()) << "SRR490124-4pairs outskirt_content_[0][2].size() not correct for " << context << '\n';
	EXPECT_EQ(0, test.outskirt_content_[0][3].size()) << "SRR490124-4pairs outskirt_content_[0][3].size() not correct for " << context << '\n';
	EXPECT_EQ(37, test.outskirt_content_[1][0].size()) << "SRR490124-4pairs outskirt_content_[1][0].size() not correct for " << context << '\n';
	EXPECT_EQ(36, test.outskirt_content_[1][1].size()) << "SRR490124-4pairs outskirt_content_[1][1].size() not correct for " << context << '\n';
	EXPECT_EQ(40, test.outskirt_content_[1][2].size()) << "SRR490124-4pairs outskirt_content_[1][2].size() not correct for " << context << '\n';
	EXPECT_EQ(37, test.outskirt_content_[1][3].size()) << "SRR490124-4pairs outskirt_content_[1][3].size() not correct for " << context << '\n';
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
	EXPECT_EQ(1, test.abundance_[2] ) << "abundance_ wrong in coverage test\n";
	EXPECT_EQ(1, test.abundance_[3] ) << "abundance_ wrong in coverage test\n";
	EXPECT_EQ(1, test.abundance_[8] ) << "abundance_ wrong in coverage test\n";
}
// Adapter with mapping quality 8 still used for adapters
void FragmentDistributionStatsTest::TestAdapters(const FragmentDistributionStats &test){
	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(177, test.insert_lengths_.size()) << "insert_lengths_.size() wrong with adapters\n";
	EXPECT_EQ(2, test.insert_lengths_[0]) << "insert_lengths_[0] wrong with adapters\n";
	EXPECT_EQ(1, test.insert_lengths_[58]) << "insert_lengths_[58] wrong with adapters\n";
	EXPECT_EQ(1, test.insert_lengths_[80]) << "insert_lengths_[80] wrong with low mapping quality adapters\n";
	EXPECT_EQ(1, test.insert_lengths_[81]) << "insert_lengths_[81] wrong with adapters\n";
	EXPECT_EQ(1, test.insert_lengths_[95]) << "insert_lengths_[95] wrong with adapters\n";
	EXPECT_EQ(1, test.insert_lengths_[176]) << "insert_lengths_[176] wrong with adapters\n";
	EXPECT_EQ(7, SumVect(test.insert_lengths_)) << "Sum of insert_lengths_ wrong with adapters\n";

	// Check if adapter part of read has been ignored for these statistics
	// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 16 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8 "-" $4+length($10)-1)}' | seqtk seq | awk '(0==NR%2){len=length($0); print int(gsub(/[CG]/,"",$0)*100/len+0.5)}' | sort -n
	EXPECT_EQ(7, test.gc_fragment_content_.size() );
	EXPECT_EQ(1, test.gc_fragment_content_[47] );
	EXPECT_EQ(1, test.gc_fragment_content_[53] );

	// cat <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 144 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq) <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-10 "-" $4+length($10)-1+10)}' | seqtk seq -r) | awk '(0==NR%2){print ">0", substr($0,1,10); print substr($0,length($0)-9,10); print ">1", substr($0,11,10); print substr($0,length($0)-19,10);print ">2", substr($0,21,10); print substr($0,length($0)-29,10);}' | seqtk seq -r | awk '{if(1==NR%2){sur=substr($0, 2, 1); print "start", sur, $2}else{print "end", sur, $0}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}' | sort
	EXPECT_EQ(1, test.fragment_surroundings_[0][817646]);
	EXPECT_EQ(1, test.fragment_surroundings_[0][940919]);
	EXPECT_EQ(1, test.fragment_surroundings_[1][10510]);
	EXPECT_EQ(1, test.fragment_surroundings_[1][477197]);
	EXPECT_EQ(1, test.fragment_surroundings_[2][495577]);
	EXPECT_EQ(1, test.fragment_surroundings_[2][718722]);
	EXPECT_EQ(1, test.fragment_surroundings_[0][301332]);
	EXPECT_EQ(1, test.fragment_surroundings_[0][824719]);
	EXPECT_EQ(1, test.fragment_surroundings_[1][339634]);
	EXPECT_EQ(1, test.fragment_surroundings_[1][595486]);
	EXPECT_EQ(1, test.fragment_surroundings_[2][158717]);
	EXPECT_EQ(1, test.fragment_surroundings_[2][289540]);

	// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 80 -F 32 | awk '($4+2000>$8){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8-20 "-" $4+length($10)-1+20)}' | seqtk seq -r | awk '(0==NR%2){print substr($0,1,20) substr($0,length($0)-19,20)}' | awk '{for(pos=1;pos<=length($0);++pos){print 1, substr($0,pos,1), pos-1}}' | sort -k1,2 -k3,3n
	EXPECT_EQ(31, test.outskirt_content_[1][1].size() );
	EXPECT_EQ(1, test.outskirt_content_[1][1][4] );
	EXPECT_EQ(1, test.outskirt_content_[1][1][7] );
	EXPECT_EQ(1, test.outskirt_content_[1][1][31] );
	EXPECT_EQ(1, test.outskirt_content_[1][1][34] );
}

void FragmentDistributionStatsTest::BiasCalculationThread(FragmentDistributionStats &test, const Reference &reference, FragmentDuplicationStats &duplications, BiasCalculationVectors &thread_values, mutex &print_mutex){
	test.ExecuteBiasCalculations( reference, duplications, thread_values, print_mutex );
}

namespace reseq{
	TEST_F(FragmentDistributionStatsTest, SeparatingAndCombiningBias){
		LoadReference("reference-test.fa");
		CreateTestObject(&species_reference_);

		test_->abundance_[0] = 7084; // bias=1.0

		test_->insert_lengths_[10] = 7084; // bias=1.0
		test_->insert_lengths_.Shrink();

		test_->gc_fragment_content_[0] = 5200; // bias=1.0
		test_->gc_fragment_content_[30] = 384; // bias=0.1
		test_->gc_fragment_content_[40] = 1500; // bias=0.3
		test_->gc_fragment_content_.Shrink();

		// Fragment start
		test_->fragment_surroundings_[0][954112] = 1500; // TGGATTAAAA bias=1.0
		test_->fragment_surroundings_[0][771557] = 192; // GTTACCTGCC bias=1.0
		test_->fragment_surroundings_[0][756483] = 3200; // GTGAGTAAAT bias=1.0
		test_->fragment_surroundings_[0][594450] = 2000; // GCACAGACAG bias=0.4
		test_->fragment_surroundings_[0][258785] = 192; // ATTTAGTGAC bias=0.8
		test_->fragment_surroundings_[1][8941] = 1500; // AAAGAGTGTC bias=1.0
		test_->fragment_surroundings_[1][756483] = 192; // GTGAGTAAAT bias=1.0
		test_->fragment_surroundings_[1][787452] = 2*3200; // TAAAATTTTA bias=0.8
		test_->fragment_surroundings_[1][196668] = 2000; // ATAAAAATTA bias=1.0
		test_->fragment_surroundings_[1][461635] = 192; // CTAAGTCAAT bias=1.0
		test_->fragment_surroundings_[2][930377] = 1500; // TGATAGCAGC bias=1.0
		test_->fragment_surroundings_[2][787452] = 192; // TAAAATTTTA bias=0.8
		test_->fragment_surroundings_[2][1017802] = 3200; // TTGACTTAGG bias=1.0
		test_->fragment_surroundings_[2][297745] = 2000; // CAGAGTACAC bias=1.0
		test_->fragment_surroundings_[2][4080] = 192; // AAAATTTTAA bias=0.6

		// Fragment end
		test_->fragment_surroundings_[0][649012] = 1500; // GCTGCTATCA bias=1.0
		test_->fragment_surroundings_[0][787452] = 192; // TAAAATTTTA bias=0.8
		test_->fragment_surroundings_[0][377552] = 3200; // CCTAAGTCAA bias=1.0
		test_->fragment_surroundings_[0][766430] = 2000; // GTGTACTCTG bias=1.0
		test_->fragment_surroundings_[0][983295] = 192; // TTAAAATTTT bias=1.0
		test_->fragment_surroundings_[1][542591] = 1500; // GACACTCTTT bias=1.0
		test_->fragment_surroundings_[1][258513] = 192; // ATTTACTCAC bias=1.0
		//test_->fragment_surroundings_[1][787452] = 2*3200; // TAAAATTTTA bias=0.8; it is also a start surrounding
		test_->fragment_surroundings_[1][802803] = 2000; // TAATTTTTAT bias=1.0
		test_->fragment_surroundings_[1][254450] = 192; // ATTGACTTAG bias=1.0
		test_->fragment_surroundings_[2][1044692] = 1500; // TTTTAATCCA bias=1.0
		test_->fragment_surroundings_[2][674497] = 192; // GGCAGGTAAC bias=0.6
		test_->fragment_surroundings_[2][258513] = 3200; // ATTTACTCAC bias=1.0
		test_->fragment_surroundings_[2][505785] = 2000; // CTGTCTGTGC bias=1.0
		test_->fragment_surroundings_[2][739075] = 192; // GTCACTAAAT bias=0.8

		vector<double> separated_surroundings;
		test_->SeparateSurroundingPositions(separated_surroundings, test_->fragment_surroundings_);
		double mean = (192 + 3200 + 192+3200+2000+1500+2000 + 1500+192+192)/4.0;
		EXPECT_DOUBLE_EQ((192-mean)/IntPow(4,9), separated_surroundings[0]);
		EXPECT_DOUBLE_EQ((3200-mean)/IntPow(4,9), separated_surroundings[1]);
		EXPECT_DOUBLE_EQ((192+3200+2000+1500+2000-mean)/IntPow(4,9), separated_surroundings[2]);
		EXPECT_DOUBLE_EQ((1500+192+192-mean)/IntPow(4,9), separated_surroundings[3]);
		mean = (1500+1500+192+3200 + 192+192 + 2000+2000 + 3200+192)/4.0;
		EXPECT_DOUBLE_EQ((1500+1500+192+3200-mean)/IntPow(4,9), separated_surroundings[36]);
		EXPECT_DOUBLE_EQ((192+192-mean)/IntPow(4,9), separated_surroundings[37]);
		EXPECT_DOUBLE_EQ((2000+2000-mean)/IntPow(4,9), separated_surroundings[38]);
		EXPECT_DOUBLE_EQ((3200+192-mean)/IntPow(4,9), separated_surroundings[39]);

		mean = (1500+2000+192+192 + 192 + 192+1500 + 2*3200+2000)/4.0;
		EXPECT_EQ((1500+2000+192+192-mean)/IntPow(4,9), separated_surroundings[40]);
		EXPECT_EQ((192-mean)/IntPow(4,9), separated_surroundings[41]);
		EXPECT_EQ((192+1500-mean)/IntPow(4,9), separated_surroundings[42]);
		EXPECT_EQ((2*3200+2000-mean)/IntPow(4,9), separated_surroundings[43]);
		mean = (2*3200+2000 + 1500+192 + 192 + 192+192+1500+2000)/4.0;
		EXPECT_EQ((2*3200+2000-mean)/IntPow(4,9), separated_surroundings[76]);
		EXPECT_EQ((1500+192-mean)/IntPow(4,9), separated_surroundings[77]);
		EXPECT_EQ((192-mean)/IntPow(4,9), separated_surroundings[78]);
		EXPECT_EQ((192+192+1500+2000-mean)/IntPow(4,9), separated_surroundings[79]);

		mean = (192+3200 + 2000+2000 + 192+192 + 1500+192+3200+1500)/4.0;
		EXPECT_EQ((192+3200-mean)/IntPow(4,9), separated_surroundings[80]);
		EXPECT_EQ((2000+2000-mean)/IntPow(4,9), separated_surroundings[81]);
		EXPECT_EQ((192+192-mean)/IntPow(4,9), separated_surroundings[82]);
		EXPECT_EQ((1500+192+3200+1500-mean)/IntPow(4,9), separated_surroundings[83]);
		mean = (192+192+1500 + 1500+2000+192+3200+2000 + 3200 + 192)/4.0;
		EXPECT_EQ((192+192+1500-mean)/IntPow(4,9), separated_surroundings[116]);
		EXPECT_EQ((1500+2000+192+3200+2000-mean)/IntPow(4,9), separated_surroundings[117]);
		EXPECT_EQ((3200-mean)/IntPow(4,9), separated_surroundings[118]);
		EXPECT_EQ((192-mean)/IntPow(4,9), separated_surroundings[119]);

		array<double, 4*Reference::num_surrounding_blocks_*Reference::surrounding_range_> separated_bias;
		separated_bias.fill(0.0);
		for(uint16_t block=0; block <= 80; block += 40){
			//ACGTTGCATA: 114252
			separated_bias[block+0] = 0.9;
			separated_bias[block+5] = 1.0;
			separated_bias[block+10] = 1.0;
			separated_bias[block+15] = 1.0;
			separated_bias[block+19] = 1.0;
			separated_bias[block+22] = 1.0;
			separated_bias[block+25] = 1.0;
			separated_bias[block+28] = 1.0;
			separated_bias[block+35] = 1.0;
			separated_bias[block+36] = 1.0;
			//ACGTT[C]CATA: 113996
			separated_bias[block+21] = 0.8;
		}

		array<vector<double>, Reference::num_surrounding_blocks_> combined_bias;
		test_->CombineSurroundingPositions(combined_bias, separated_bias);

		for(uint16_t block=0; block < combined_bias.size(); ++block){
			EXPECT_DOUBLE_EQ(9.9, combined_bias[block][114252]);
			EXPECT_DOUBLE_EQ(9.7, combined_bias[block][113996]);
		}
	}

	TEST_F(FragmentDistributionStatsTest, BiasCalculationVectors){
		LoadReference("reference-test-withN.fa");

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

		//test.GetWeights();
		double coverage = static_cast<double>(5)/(19*2);
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

	TEST_F(FragmentDistributionStatsTest, BiasCalculation){
		// samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:100-199 | awk '(NR==1){print $0 "_1000times"}(NR==2){for(i=0;i<1000;i++){print $0}}' | seqtk seq > reference-bias-calculation.fa
		LoadReference("reference-bias-calculation.fa");
		CreateTestObject(&species_reference_);

		// To set how much space is needed
		test_->insert_lengths_[250] = 0;
		test_->abundance_.resize(species_reference_.NumberSequences(), 0);

		// Create data by the model that is used for the bias fitting (so it should be easy)
		test_->fragment_sites_by_ref_seq_bin_.resize(species_reference_.NumberSequences());
		test_->fragment_sites_by_ref_seq_bin_cur_id_.resize(species_reference_.NumberSequences());

		uint32_t last_size = 0;
		for(uint32_t fragment_length=65; fragment_length <= 95; fragment_length += 5){
			CreateCoverageData(fragment_length);

			test_->insert_lengths_.at(fragment_length) += test_->fragment_sites_by_ref_seq_bin_.at(0).size() - last_size;
			last_size = test_->fragment_sites_by_ref_seq_bin_.at(0).size();
		}
		test_->abundance_.at(0) += test_->fragment_sites_by_ref_seq_bin_.at(0).size();

		// Fit data
		FragmentDuplicationStats duplications;
		duplications.PrepareTmpDuplicationVector(100);
		uint32_t num_threads = 4;
		array<FragmentDistributionStats::ThreadData, 4> thread_data = {FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0)), FragmentDistributionStats::ThreadData(100, species_reference_.SequenceLength(0))};
		mutex print_mutex;

		uint16_t cur_verbosity = kVerbosityLevel;
		kVerbosityLevel = 1; // Suppress warnings
		test_->AddNewBiasCalculations(1, thread_data.at(0), print_mutex);

		thread threads[num_threads];
		for(auto i = num_threads; i--; ){
			threads[i] = thread(BiasCalculationThread, std::ref(*test_), std::cref(species_reference_), std::ref(duplications), std::ref(thread_data.at(i).bias_calc_vects_), std::ref(print_mutex));
		}
		for(auto i = num_threads; i--; ){
			threads[i].join();
		}
		test_->FinalizeBiasCalculation(species_reference_, 4, duplications);
		kVerbosityLevel = cur_verbosity;

		EXPECT_NEAR(0.6, test_->gc_fragment_content_bias_.at(27), 0.1);
		EXPECT_NEAR(0.7, test_->gc_fragment_content_bias_.at(28), 0.1);
		EXPECT_NEAR(0.8, test_->gc_fragment_content_bias_.at(29), 0.1);
		EXPECT_NEAR(0.9, test_->gc_fragment_content_bias_.at(30), 0.1);
		EXPECT_NEAR(1.0, test_->gc_fragment_content_bias_.at(31), 0.1);
		EXPECT_NEAR(0.8, test_->gc_fragment_content_bias_.at(32), 0.1);
		EXPECT_NEAR(0.6, test_->gc_fragment_content_bias_.at(33), 0.2);
		EXPECT_NEAR(0.4, test_->gc_fragment_content_bias_.at(34), 0.2);
		EXPECT_DOUBLE_EQ(1.0, *max_element(test_->gc_fragment_content_bias_.begin(), test_->gc_fragment_content_bias_.end()));

		uint32_t pos_base = 262144;
		EXPECT_NEAR(-0.4, accumulate(test_->fragment_surroundings_bias_.at(1).begin(), test_->fragment_surroundings_bias_.at(1).begin()+pos_base, 0.0)/pos_base, 0.2);
		EXPECT_NEAR(-0.4, accumulate(test_->fragment_surroundings_bias_.at(1).begin()+pos_base, test_->fragment_surroundings_bias_.at(1).begin()+2*pos_base, 0.0)/pos_base, 0.2);
		EXPECT_NEAR(0.4, accumulate(test_->fragment_surroundings_bias_.at(1).begin()+2*pos_base, test_->fragment_surroundings_bias_.at(1).begin()+3*pos_base, 0.0)/pos_base, 0.2);
		EXPECT_NEAR(0.4, accumulate(test_->fragment_surroundings_bias_.at(1).begin()+3*pos_base, test_->fragment_surroundings_bias_.at(1).begin()+4*pos_base, 0.0)/pos_base, 0.2);

		for(uint32_t fragment_length=65; fragment_length <= 95; fragment_length += 5){
			EXPECT_NEAR(1.0, test_->insert_lengths_bias_.at(fragment_length), 0.1) << "Fragment length " << fragment_length << " wrong\n";
		}
		EXPECT_DOUBLE_EQ(1.0, *max_element(test_->insert_lengths_bias_.begin(), test_->insert_lengths_bias_.end()));

		EXPECT_DOUBLE_EQ(1.0, test_->ref_seq_bias_.at(0));

		EXPECT_NEAR(0.5, test_->dispersion_parameters_.at(0), 0.4);
		EXPECT_NEAR(1.0, test_->dispersion_parameters_.at(1), 0.1);
	}

	TEST_F(FragmentDistributionStatsTest, UpdateRefSeqBias){
		LoadReference("reference-test.fa");
		CreateTestObject(&species_reference_);
		std::mt19937_64 rgen;

		uint16_t cur_verbosity = kVerbosityLevel;
		kVerbosityLevel = 1; // Suppress warnings

		test_->ref_seq_bias_.resize(1, 0.5);
		EXPECT_TRUE(test_-> UpdateRefSeqBias(RefSeqBiasSimulation::kKeep, std::string(), species_reference_, rgen));
		EXPECT_EQ(2, test_->ref_seq_bias_.size()) << "Keeping reference bias does not switch to no bias correctly.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(0)) << "Keeping reference bias does not switch to no bias correctly.";
		EXPECT_EQ(1.0, test_->ref_seq_bias_.at(1)) << "Keeping reference bias does not switch to no bias correctly.";

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

		kVerbosityLevel = cur_verbosity;
	}
}
