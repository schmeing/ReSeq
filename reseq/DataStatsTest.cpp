#include "DataStatsTest.h"
using reseq::DataStatsTest;

#include <sstream>
using std::stringstream;
//include <string>
using std::string;
using std::to_string;

#include <seqan/bam_io.h>
using seqan::CharString;
using seqan::toCString;

#include "AdapterStatsTest.h"
using reseq::AdapterStatsTest;
#include "CoverageStatsTest.h"
using reseq::CoverageStatsTest;
#include "ErrorStatsTest.h"
using reseq::ErrorStatsTest;
#include "FragmentDistributionStatsTest.h"
using reseq::FragmentDistributionStatsTest;
#include "FragmentDuplicationStatsTest.h"
using reseq::FragmentDuplicationStatsTest;
#include "QualityStatsTest.h"
using reseq::QualityStatsTest;
#include "TileStatsTest.h"
using reseq::TileStatsTest;

void DataStatsTest::Register(){
	// Guarantees that library is included
}

void DataStatsTest::CreateTestObject(Reference *ref){
	ASSERT_TRUE( test_ = new DataStats(ref) ) << "Could not allocate memory for DataStats object\n";
}

void DataStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void DataStatsTest::LoadStats(const string &stats_file, bool ignore_tiles, bool calculate_biases, const string adapter, const string variants){
	if(ignore_tiles){
		test_->tiles_.IgnoreTiles();
	}

	ASSERT_TRUE( test_->ReadBam( (string(PROJECT_SOURCE_DIR)+"/test/"+stats_file).c_str(), (string(PROJECT_SOURCE_DIR)+"/adapters/"+adapter+".fa").c_str(), (string(PROJECT_SOURCE_DIR)+"/adapters/"+adapter+".mat").c_str(), (variants.empty()?variants:string(PROJECT_SOURCE_DIR)+"/test/"+variants), 100000000, 1, calculate_biases ) );
	test_->PrepareTesting();
}

void DataStatsTest::TestSequenceContent(
		uintTempSeq template_segment,
		uintReadLen at_pos,
		uintNucCount cont_a,
		uintNucCount cont_c,
		uintNucCount cont_g,
		uintNucCount cont_t,
		uintNucCount cont_n,
		const char * context ){
	EXPECT_EQ(cont_a, test_->sequence_content_[template_segment][0][at_pos]) << "sequence_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_c, test_->sequence_content_[template_segment][1][at_pos]) << "sequence_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_g, test_->sequence_content_[template_segment][2][at_pos]) << "sequence_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_t, test_->sequence_content_[template_segment][3][at_pos]) << "sequence_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_n, test_->sequence_content_[template_segment][4][at_pos]) << "sequence_content_[" << template_segment << "] position " << at_pos << " wrong for " << context << '\n';
}

void DataStatsTest::TestSequenceContentReference(
		uintTempSeq template_segment,
		uintTempSeq strand,
		uintReadLen at_pos,
		uintNucCount cont_a,
		uintNucCount cont_c,
		uintNucCount cont_g,
		uintNucCount cont_t,
		const char * context ){
	EXPECT_EQ(cont_a, test_->sequence_content_reference_[template_segment][strand][0][at_pos]) << "sequence_content_reference_[" << template_segment << "][" << strand << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_c, test_->sequence_content_reference_[template_segment][strand][1][at_pos]) << "sequence_content_reference_[" << template_segment << "][" << strand << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_g, test_->sequence_content_reference_[template_segment][strand][2][at_pos]) << "sequence_content_reference_[" << template_segment << "][" << strand << "] position " << at_pos << " wrong for " << context << '\n';
	EXPECT_EQ(cont_t, test_->sequence_content_reference_[template_segment][strand][3][at_pos]) << "sequence_content_reference_[" << template_segment << "][" << strand << "] position " << at_pos << " wrong for " << context << '\n';
}

void DataStatsTest::TestSrr490124Equality(const char *context, bool test_tile_information){
	EXPECT_TRUE( test_->reference_ ) << "SRR490124-4pairs reference_ wrong for " << context << '\n';
	EXPECT_EQ(1, test_->reference_->NumberSequences() ) << "SRR490124-4pairs reference_ wrong for " << context << '\n';
	EXPECT_TRUE(test_->reference_->ReferenceIdFirstPart(0) == "NC_000913.3") << test_->reference_->ReferenceIdFirstPart(0) << " wrong for SRR490124-4pairs reference_ for " << context << '\n';

	CoverageStatsTest::TestSrr490124Equality(test_->coverage_, context);
	EXPECT_EQ(0, test_->first_read_records_.size() ) << "SRR490124-4pairs first_read_records_ wrong for " << context << '\n';

	ErrorStatsTest::TestSrr490124Equality(test_->errors_, context);
	FragmentDistributionStatsTest::TestSrr490124Equality(test_->fragment_distribution_, context);
	FragmentDuplicationStatsTest::TestSrr490124Equality(test_->duplicates_, context);
	QualityStatsTest::TestSrr490124Equality(test_->qualities_, context);
	TileStatsTest::TestSrr490124Equality(test_->tiles_, context, test_tile_information);

	EXPECT_EQ(8, test_->total_number_reads_) << "SRR490124-4pairs total_number_reads_ wrong for " << context << '\n';
	for( int templ_seg=2; templ_seg--; ){
		EXPECT_EQ(1, test_->read_lengths_[templ_seg].size()) << "SRR490124-4pairs read_lengths_[" << templ_seg << "] wrong for " << context << '\n';
		EXPECT_EQ(4, test_->read_lengths_[templ_seg][100]) << "SRR490124-4pairs read_lengths_[" << templ_seg << "] wrong for " << context << '\n';
	}

	EXPECT_EQ(19, test_->proper_pair_mapping_quality_.size()) << "SRR490124-4pairs proper_pair_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[6]) << "SRR490124-4pairs proper_pair_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[23]) << "SRR490124-4pairs proper_pair_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[24]) << "SRR490124-4pairs proper_pair_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(0, test_->improper_pair_mapping_quality_.size()) << "SRR490124-4pairs improper_pair_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(1, test_->single_read_mapping_quality_.size()) << "SRR490124-4pairs single_read_mapping_quality_ wrong for " << context << '\n';
	EXPECT_EQ(1, test_->single_read_mapping_quality_[3]) << "SRR490124-4pairs single_read_mapping_quality_ wrong for " << context << '\n';

	// samtools view ecoli-SRR490124-4pairs.sam | awk '{print int($2%256/128), gsub("G","",$10)+gsub("C","",$10)}' | sort -n
	EXPECT_EQ(12, test_->gc_read_content_[0].size()) << "SRR490124-4pairs gc_read_content_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[0][48]) << "SRR490124-4pairs gc_read_content_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[0][52]) << "SRR490124-4pairs gc_read_content_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[0][53]) << "SRR490124-4pairs gc_read_content_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[0][59]) << "SRR490124-4pairs gc_read_content_[0] wrong for " << context << '\n';
	EXPECT_EQ(20, test_->gc_read_content_[1].size()) << "SRR490124-4pairs gc_read_content_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[1][37]) << "SRR490124-4pairs gc_read_content_[1] wrong for " << context << '\n';
	EXPECT_EQ(2, test_->gc_read_content_[1][47]) << "SRR490124-4pairs gc_read_content_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_[1][56]) << "SRR490124-4pairs gc_read_content_[1] wrong for " << context << '\n';
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{print int($2%256/128), $4}' | sort -n | awk '{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2 "-" $2+99)}' | seqtk seq | awk '(0==NR%2){print gsub("G","",$0)+gsub("C","",$0)}'
	EXPECT_EQ(4, test_->gc_read_content_reference_[0].size()) << "SRR490124-4pairs gc_read_content_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_reference_[0][54]) << "SRR490124-4pairs gc_read_content_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_reference_[0][57]) << "SRR490124-4pairs gc_read_content_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(7, test_->gc_read_content_reference_[1].size()) << "SRR490124-4pairs gc_read_content_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_reference_[1][46]) << "SRR490124-4pairs gc_read_content_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_reference_[1][52]) << "SRR490124-4pairs gc_read_content_reference_[1] wrong for " << context << '\n';
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{print int($2%256/128), gsub("G","",$10)+gsub("C","",$10)}' | sort -n
	EXPECT_EQ(7, test_->gc_read_content_mapped_[0].size()) << "SRR490124-4pairs gc_read_content_mapped_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_mapped_[0][53]) << "SRR490124-4pairs gc_read_content_mapped_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_mapped_[0][59]) << "SRR490124-4pairs gc_read_content_mapped_[0] wrong for " << context << '\n';
	EXPECT_EQ(10, test_->gc_read_content_mapped_[1].size()) << "SRR490124-4pairs gc_read_content_mapped_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_mapped_[1][47]) << "SRR490124-4pairs gc_read_content_mapped_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test_->gc_read_content_mapped_[1][56]) << "SRR490124-4pairs gc_read_content_mapped_[1] wrong for " << context << '\n';
	for( int templ_seg=2; templ_seg--; ){
		EXPECT_EQ(1, test_->n_content_[templ_seg].size()) << "SRR490124-4pairs n_content_[" << templ_seg << "] wrong for " << context << '\n';
		EXPECT_EQ(4, test_->n_content_[templ_seg][0]) << "SRR490124-4pairs n_content_[" << templ_seg << "] wrong for " << context << '\n';
	}

	// samtools view ecoli-SRR490124-4pairs.sam | awk '(0==int($2%256/128)){print ">a\n" $10}' | seqtk seq -r | awk '(0==NR%2){print substr($0,1,1)}' | sort | uniq -c
	TestSequenceContent(0, 0, 2, 1, 1, 0, 0, context);
	TestSequenceContent(0, 1, 1, 2, 0, 1, 0, context);
	TestSequenceContent(0, 2, 1, 2, 1, 0, 0, context);
	TestSequenceContent(0, 3, 0, 0, 3, 1, 0, context);
	TestSequenceContent(0, 4, 2, 0, 1, 1, 0, context);
	TestSequenceContent(0, 95, 2, 0, 1, 1, 0, context);
	TestSequenceContent(0, 96, 2, 0, 1, 1, 0, context);
	TestSequenceContent(0, 97, 2, 0, 1, 1, 0, context);
	TestSequenceContent(0, 98, 1, 0, 1, 2, 0, context);
	TestSequenceContent(0, 99, 2, 0, 1, 1, 0, context);
	// samtools view ecoli-SRR490124-4pairs.sam | awk '(1==int($2%256/128)){print substr($10,1,1)}' | sort | uniq -c
	TestSequenceContent(1, 0, 1, 1, 1, 1, 0, context);
	TestSequenceContent(1, 1, 2, 0, 0, 2, 0, context);
	TestSequenceContent(1, 2, 2, 0, 1, 1, 0, context);
	TestSequenceContent(1, 3, 1, 1, 0, 2, 0, context);
	TestSequenceContent(1, 4, 2, 2, 0, 0, 0, context);
	TestSequenceContent(1, 95, 0, 3, 0, 1, 0, context);
	TestSequenceContent(1, 96, 0, 1, 3, 0, 0, context);
	TestSequenceContent(1, 97, 0, 3, 1, 0, 0, context);
	TestSequenceContent(1, 98, 1, 0, 0, 3, 0, context);
	TestSequenceContent(1, 99, 0, 1, 1, 2, 0, context);

	// Position needs to be updated at the end (starting from 1 not from 0 as in the test function below)
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(0==int($2%256/128)){print $4}' | sort -n | awk '{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $1 "-" $1+99)}' | seqtk seq -r | awk '(0==NR%2){print substr($0,1,1)}' | sort | uniq -c
	TestSequenceContentReference(0, 1, 0, 0, 1, 1, 0, context);
	TestSequenceContentReference(0, 1, 1, 0, 1, 0, 1, context);
	TestSequenceContentReference(0, 1, 2, 0, 1, 1, 0, context);
	TestSequenceContentReference(0, 1, 3, 0, 0, 2, 0, context);
	TestSequenceContentReference(0, 1, 4, 0, 0, 1, 1, context);
	TestSequenceContentReference(0, 1, 95, 1, 0, 1, 0, context);
	TestSequenceContentReference(0, 1, 96, 0, 0, 2, 0, context);
	TestSequenceContentReference(0, 1, 97, 1, 0, 1, 0, context);
	TestSequenceContentReference(0, 1, 98, 1, 0, 0, 1, context);
	TestSequenceContentReference(0, 1, 99, 2, 0, 0, 0, context);
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(1==int($2%256/128)){print $4}' | sort -n | awk '{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $1 "-" $1+99)}' | seqtk seq | awk '(0==NR%2){print substr($0,1,1)}' | sort | uniq -c
	TestSequenceContentReference(1, 0, 0, 0, 0, 1, 1, context);
	TestSequenceContentReference(1, 0, 1, 1, 0, 0, 1, context);
	TestSequenceContentReference(1, 0, 2, 1, 0, 1, 0, context);
	TestSequenceContentReference(1, 0, 3, 0, 1, 0, 1, context);
	TestSequenceContentReference(1, 0, 4, 1, 1, 0, 0, context);
	TestSequenceContentReference(1, 0, 95, 0, 2, 0, 0, context);
	TestSequenceContentReference(1, 0, 96, 1, 0, 1, 0, context);
	TestSequenceContentReference(1, 0, 97, 0, 1, 1, 0, context);
	TestSequenceContentReference(1, 0, 98, 0, 1, 0, 1, context);
	TestSequenceContentReference(1, 0, 99, 1, 0, 0, 1, context);
	for(int base=4; base--; ){
		EXPECT_EQ(0, test_->sequence_content_reference_[0][0][base].size()) << "SRR490124-4pairs sequence_content_reference_[0][0][" << base << "].size() not correct for " << context << '\n';
		EXPECT_EQ(0, test_->sequence_content_reference_[1][1][base].size()) << "SRR490124-4pairs sequence_content_reference_[1][1][" << base << "].size() not correct for " << context << '\n';
	}

	// cat ecoli-SRR490124-4pairs-R1.fq ecoli-SRR490124-4pairs-R2.fq | awk -v FS="" '(2==NR%4){count=0;base="";for (i=1;i<=NF;i++){ if(base == $i){count+=1}else{if("" != base){print base, count};count=1;base=$i}};print base, count}' | sort | uniq -c
	TestVectEquality({1,{107,20,9,5,0,0,3}}, test_->homopolymer_distribution_[0], context, "SRR490124-4pairs homopolymer_distribution_[0]", " not correct for ");
	TestVectEquality({1,{131,30,3,3}}, test_->homopolymer_distribution_[1], context, "SRR490124-4pairs homopolymer_distribution_[1]", " not correct for ");
	TestVectEquality({1,{122,21,5,2}}, test_->homopolymer_distribution_[2], context, "SRR490124-4pairs homopolymer_distribution_[2]", " not correct for ");
	TestVectEquality({1,{93,31,5,4}}, test_->homopolymer_distribution_[3], context, "SRR490124-4pairs homopolymer_distribution_[3]", " not correct for ");
	EXPECT_EQ(0, test_->homopolymer_distribution_[4].size()) << "SRR490124-4pairs homopolymer_distribution_[4].size() not correct for " << context << '\n';
}

void DataStatsTest::TestTiles(){
	TileStatsTest::TestTiles(test_->tiles_);
	QualityStatsTest::TestTiles(test_->qualities_);

	// samtools view ecoli-tiles.bam | awk 'BEGIN{sum=0}(0==int($2%256/128)){sum += gsub("N","",$10)}END{print sum}'
	EXPECT_EQ(5, SumVect(test_->sequence_content_[0][4])) << "sequence_content_[0][4] wrong in tile test\n";
	// samtools view ecoli-tiles.bam | awk '(0==int($2%256/128)){print substr($10,47,1)}' | sort | uniq -c
	EXPECT_EQ(1, test_->sequence_content_[0][4][46]) << "sequence_content_[0][4] wrong in tile test\n";
	EXPECT_EQ(1, test_->sequence_content_[0][4][66]) << "sequence_content_[0][4] wrong in tile test\n";

	// cat ecoli-tiles-R1.fq ecoli-tiles-R2.fq | awk -v FS="" '(2==NR%4){count=0;base="";for (i=1;i<=NF;i++){ if(base == $i){count+=1}else{if("" != base){print base, count};count=1;base=$i}};print base, count}' | sort | uniq -c
	TestVectEquality({2,{1,1}}, test_->homopolymer_distribution_[4], "tile test", "homopolymer_distribution_[4]", " not correct for ");
}

void DataStatsTest::TestDuplicates(){
	CoverageStatsTest::TestDuplicates(test_->coverage_);
	EXPECT_EQ(0, test_->first_read_records_.size() ) << "first_read_records_ wrong in duplicates test\n";

	ErrorStatsTest::TestDuplicates(test_->errors_);
	FragmentDistributionStatsTest::TestDuplicates(test_->fragment_distribution_);
	FragmentDuplicationStatsTest::TestDuplicates(test_->duplicates_);
	QualityStatsTest::TestDuplicates(test_->qualities_);

	// No single mappings
	// samtools view ecoli-duplicates.bam | awk '{print int($2%4/2), $5}' | sort -k1,1n -k2,2n | uniq -c
	EXPECT_EQ(42, test_->proper_pair_mapping_quality_.size()) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[1]) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[23]) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(2, test_->proper_pair_mapping_quality_[24]) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(6, test_->proper_pair_mapping_quality_[40]) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(22, test_->proper_pair_mapping_quality_[42]) << "proper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(1, test_->improper_pair_mapping_quality_.size()) << "improper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(2, test_->improper_pair_mapping_quality_[42]) << "improper_pair_mapping_quality_ wrong in duplicates test\n";
	EXPECT_EQ(0, test_->single_read_mapping_quality_.size()) << "single_read_mapping_quality_ wrong in duplicates test\n";

	// samtools view -q 10 -f 99 ecoli-duplicates.bam | awk '(0 != substr($1,12,1)){print $6; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $4 "-" $4+length($10)-1)}' | awk 'BEGIN{last=""}{if(">" == substr($0,1,1)){print ">" last; last=""}else{print last; last = $0}}END{print last}' | seqtk seq | awk '(1==NR%2){cigar = substr($0,2,length($0)-1)}(0==NR%2){pos=1; print_pos=0; num=0; for(i=1;i<=length(cigar);i+=1){e=substr(cigar,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e){while(0<num){print print_pos, substr($0,pos,1); ++print_pos; ++pos; --num}};if("D"==e){pos+=num;num=0};if("I"==e){print_pos+=num;num=0}}}}' | awk '($1 < 10 && $1 ~ /^[4-7]/)' | sort -k1,1n -k2,2 | uniq -c
	TestSequenceContentReference(0, 0, 4, 7, 0, 1, 5, "duplicates test");
	TestSequenceContentReference(0, 0, 5, 5, 0, 7, 1, "duplicates test");
	TestSequenceContentReference(0, 0, 6, 0, 1, 5, 7, "duplicates test");
	TestSequenceContentReference(0, 0, 7, 0, 0, 3, 10, "duplicates test");
	// samtools view -q 10 -f 83 ecoli-duplicates.bam | awk '(0 != substr($1,12,1)){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $4 "-" $4+length($10)-1)}' | seqtk seq -r | awk '(0==NR%2){print substr($0,5,1)}' | sort | uniq -c
	TestSequenceContentReference(0, 1, 4, 1, 0, 0, 0, "duplicates test");
}

void DataStatsTest::TestVariants(){
	// Manually modified values from TestSrr490124Equality
	CoverageStatsTest::TestVariants(test_->coverage_);
	ErrorStatsTest::TestVariants(test_->errors_);
	QualityStatsTest::TestVariants(test_->qualities_);

	EXPECT_EQ(3, test_->gc_read_content_reference_[0].size());
	EXPECT_EQ(1, test_->gc_read_content_reference_[0][55]);
	EXPECT_EQ(1, test_->gc_read_content_reference_[0][57]);
	EXPECT_EQ(9, test_->gc_read_content_reference_[1].size());
	EXPECT_EQ(1, test_->gc_read_content_reference_[1][46]);
	EXPECT_EQ(1, test_->gc_read_content_reference_[1][54]);

	EXPECT_EQ(6, test_->gc_read_content_mapped_[0].size());
	EXPECT_EQ(1, test_->gc_read_content_mapped_[0][54]);
	EXPECT_EQ(1, test_->gc_read_content_mapped_[0][59]);
	EXPECT_EQ(11, test_->gc_read_content_mapped_[1].size());
	EXPECT_EQ(1, test_->gc_read_content_mapped_[1][47]);
	EXPECT_EQ(1, test_->gc_read_content_mapped_[1][57]);

	TestSequenceContentReference(0, 1, 0, 0, 1, 0, 0, "variant test");
	TestSequenceContentReference(0, 1, 1, 0, 1, 0, 1, "variant test");
	TestSequenceContentReference(0, 1, 2, 0, 1, 1, 0, "variant test");
	TestSequenceContentReference(0, 1, 3, 0, 0, 2, 0, "variant test");
	TestSequenceContentReference(0, 1, 4, 0, 0, 1, 1, "variant test");
	TestSequenceContentReference(0, 1, 95, 1, 0, 1, 0, "variant test");
	TestSequenceContentReference(0, 1, 96, 0, 0, 2, 0, "variant test");
	TestSequenceContentReference(0, 1, 97, 1, 0, 1, 0, "variant test");
	TestSequenceContentReference(0, 1, 98, 1, 0, 0, 1, "variant test");
	TestSequenceContentReference(0, 1, 99, 1, 0, 0, 0, "variant test");

	TestSequenceContentReference(1, 0, 0, 0, 0, 1, 0, "variant test");
	TestSequenceContentReference(1, 0, 1, 1, 0, 0, 1, "variant test");
	TestSequenceContentReference(1, 0, 2, 1, 0, 1, 0, "variant test");
	TestSequenceContentReference(1, 0, 3, 0, 1, 0, 1, "variant test");
	TestSequenceContentReference(1, 0, 4, 1, 1, 0, 0, "variant test");
	TestSequenceContentReference(1, 0, 95, 0, 2, 0, 0, "variant test");
	TestSequenceContentReference(1, 0, 96, 1, 0, 1, 0, "variant test");
	TestSequenceContentReference(1, 0, 97, 0, 1, 1, 0, "variant test");
	TestSequenceContentReference(1, 0, 98, 0, 1, 0, 1, "variant test");
	TestSequenceContentReference(1, 0, 99, 0, 0, 0, 1, "variant test");
	for(int base=4; base--; ){
		EXPECT_EQ(0, test_->sequence_content_reference_[0][0][base].size());
		EXPECT_EQ(0, test_->sequence_content_reference_[1][1][base].size());
	}
}

void DataStatsTest::TestCrossDuplicates(){
	EXPECT_TRUE( test_->reference_ ) << "reference_ wrong in cross duplicates test\n";
	EXPECT_EQ(1230, test_->reference_->NumberSequences() ) << "reference_ wrong in cross duplicates test\n";
	EXPECT_TRUE(test_->reference_->ReferenceIdFirstPart(0) == "NW_007931110.1") << test_->reference_->ReferenceIdFirstPart(0) << " wrong for reference_ in cross duplicates test\n";
	EXPECT_TRUE(test_->reference_->ReferenceIdFirstPart(1) == "NW_007931111.1") << test_->reference_->ReferenceIdFirstPart(1) << " wrong for reference_ in cross duplicates test\n";

	CoverageStatsTest::TestCrossDuplicates(test_->coverage_);
	EXPECT_EQ(0, test_->first_read_records_.size() ) << "first_read_records_ wrong for in cross duplicates test\n";

	FragmentDistributionStatsTest::TestCrossDuplicates(test_->fragment_distribution_);
	FragmentDuplicationStatsTest::TestCrossDuplicates(test_->duplicates_);

	EXPECT_EQ(0, test_->proper_pair_mapping_quality_.size()) << "proper_pair_mapping_quality_ wrong in cross duplicates test\n";
	EXPECT_EQ(30, test_->improper_pair_mapping_quality_.size()) << "improper_pair_mapping_quality_ wrong in cross duplicates test\n";
	EXPECT_EQ(1, test_->improper_pair_mapping_quality_[1]) << "improper_pair_mapping_quality_ wrong in cross duplicates test\n";
	EXPECT_EQ(7, test_->improper_pair_mapping_quality_[30]) << "improper_pair_mapping_quality_ wrong in cross duplicates test\n";
	EXPECT_EQ(0, test_->single_read_mapping_quality_.size()) << "single_read_mapping_quality_ wrong in cross duplicates test\n";
}

void DataStatsTest::TestCoverage(){
	CoverageStatsTest::TestCoverage(test_->coverage_);
	EXPECT_EQ(0, test_->first_read_records_.size() ) << "first_read_records_ wrong in coverage test\n";

	FragmentDistributionStatsTest::TestCoverage(test_->fragment_distribution_);
	QualityStatsTest::TestCoverage(test_->qualities_);
}

void DataStatsTest::TestAdapters(){
	AdapterStatsTest::TestAdapters(test_->adapters_);
	ErrorStatsTest::TestAdapters(test_->errors_);
	FragmentDistributionStatsTest::TestAdapters(test_->fragment_distribution_);
	QualityStatsTest::TestAdapters(test_->qualities_);

	// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 81 -F 32 | awk '{print $6; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $8 "-" $4+length($10)-1)}' | awk 'BEGIN{last=""}{if(">" == substr($0,1,1)){print ">" last; last=""}else{print last; last = $0}}END{print last}' | seqtk seq -r | awk '(1==NR%2){cigar = substr($0,2,length($0)-1)}(0==NR%2){pos=1; print_pos=0; num=0; for(i=1;i<=length(cigar);i+=1){e=substr(cigar,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e){while(0<num && pos<=length($0)){print print_pos, substr($0,pos,1); ++print_pos; ++pos; --num}};if("D"==e){pos+=num};if("I"==e){print_pos+=num};num=0}}}' | sort -k2,2 -k1,1n | uniq -c
	EXPECT_EQ(71, test_->sequence_content_reference_[0][1][0].size() );
	EXPECT_EQ(78, test_->sequence_content_reference_[0][1][1].size() );
	EXPECT_EQ(73, test_->sequence_content_reference_[0][1][2].size() );
	EXPECT_EQ(77, test_->sequence_content_reference_[0][1][3].size() );
	EXPECT_EQ(1, test_->sequence_content_reference_[0][1][3][4] );
	// samtools view ecoli-SRR490124-adapter.bam -q 10 -f 161 -F 16 | awk '{print $6; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $4 "-" $8+length($10)-1)}' | awk 'BEGIN{last=""}{if(">" == substr($0,1,1)){print ">" last; last=""}else{print last; last = $0}}END{print last}' | seqtk seq | awk '(1==NR%2){cigar = substr($0,2,length($0)-1)}(0==NR%2){pos=1; print_pos=0; num=0; for(i=1;i<=length(cigar);i+=1){e=substr(cigar,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e){while(0<num && pos<=length($0)){print print_pos, substr($0,pos,1); ++print_pos; ++pos; --num}};if("D"==e){pos+=num};if("I"==e){print_pos+=num};num=0}}}' | sort -k2,2 -k1,1n | uniq -c
	EXPECT_EQ(77, test_->sequence_content_reference_[1][0][0].size() );
	EXPECT_EQ(73, test_->sequence_content_reference_[1][0][1].size() );
	EXPECT_EQ(78, test_->sequence_content_reference_[1][0][2].size() );
	EXPECT_EQ(71, test_->sequence_content_reference_[1][0][3].size() );
}

void DataStatsTest::TestNexteraAdapters(){
	AdapterStatsTest::TestNexteraAdapters(test_->adapters_);
}

void DataStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

namespace reseq{
	TEST_F(DataStatsTest, Construction){
		CreateTestObject(NULL);

		// Constructor
		string error_msg = "Initializing of class member read_lengths_ failed\n";
		EXPECT_TRUE( test_->read_lengths_[0].empty() ) << error_msg;
		EXPECT_TRUE( test_->read_lengths_[1].empty() ) << error_msg;

		for(int i=2;i--;){
			for(int j=5;j--;){
				EXPECT_TRUE( test_->sequence_content_[i][j].empty() ) << "Initializing of class member sequence_content_ failed\n";
			}
		}

		// Shrink
		EXPECT_NO_THROW( test_->Shrink() ) << "Crashed during shrinking of empty object\n";
	}

	TEST_F(DataStatsTest, Ecoli){
		LoadReference("ecoli-GCF_000005845.2_ASM584v2_genomic.fa");
		CreateTestObject(&species_reference_);

		// zcat ecoli-SRR490124_1.fastq.gz | head -n16 > ecoli-SRR490124-4pairs-R1.fq
		// bowtie2 -X 2000 -x ecoli-GCF_000005845.2_ASM584v2_genomic.fa -1 <(reseq-prepare-names.py ecoli-SRR490124-4pairs-R1.fq ecoli-SRR490124-4pairs-R2.fq) -2 <(reseq-prepare-names.py ecoli-SRR490124-4pairs-R2.fq ecoli-SRR490124-4pairs-R1.fq) | samtools sort -m 10G -@ 4 -T _tmp -o ecoli-SRR490124-4pairs.sam -
		LoadStats("ecoli-SRR490124-4pairs.sam");
		TestSrr490124Equality("samfile");

		DeleteTestObject();
		CreateTestObject(&species_reference_);

		// samtools view -o ecoli-SRR490124-4pairs.bam ecoli-SRR490124-4pairs.sam
		LoadStats("ecoli-SRR490124-4pairs.bam");
		TestSrr490124Equality("bamfile");
		string save_test_file = string(PROJECT_SOURCE_DIR)+"/test/saveTest.reseq";
		ASSERT_TRUE( test_->Save(save_test_file.c_str()) );

		DeleteTestObject();
		CreateTestObject(&species_reference_);

		ASSERT_TRUE( test_->Load(save_test_file.c_str()) );
		test_->PrepareTesting();
		TestSrr490124Equality("save and reload", false);
		EXPECT_EQ( 0, remove(save_test_file.c_str()) ) << "Error deleting file: " << save_test_file << '\n';

		DeleteTestObject();
		CreateTestObject(&species_reference_);
		LoadStats("ecoli-tiles.bam");
		TestTiles();

		DeleteTestObject();
		CreateTestObject(&species_reference_);
		// bowtie2 -X 2000 -x ecoli-GCF_000005845.2_ASM584v2_genomic.fa -1 ecoli-duplicates-R1.fq -2 ecoli-duplicates-R2.fq | samtools sort -m 10G -@ 4 -T _tmp -o ecoli-duplicates.bam -
		LoadStats("ecoli-duplicates.bam", true);
		TestDuplicates();

		DeleteTestObject();
		CreateTestObject(&species_reference_);
		LoadStats("ecoli-SRR490124-4pairs.bam", false, false, "TruSeq_v2", "ecoli-SRR490124-4pairs.vcf");
		TestVariants();
	}

	TEST_F(DataStatsTest, Drosophila){
		LoadReference("drosophila-GCF_000001215.4_cut.fna");
		CreateTestObject(&species_reference_);

		// bowtie2 -X 2000 -x drosophila-GCF_000001215.4_cut.fna -1 drosophila-crossDuplicates-R1.fq -2 drosophila-crossDuplicates-R2.fq | samtools sort -m 10G -@ 4 -T _tmp -o drosophila-crossDuplicates.bam -
		LoadStats("drosophila-crossDuplicates.bam", true);
		TestCrossDuplicates();

		DeleteTestObject();
		CreateTestObject(&species_reference_);

		// samtools faidx drosophila-GCF_000001215.4_cut.fna NW_007931113.1 | seqtk seq | awk '(0==NR%2){print length($0)}' | awk '{for(pos=1;pos<=$0-150;++pos){system("samtools faidx drosophila-GCF_000001215.4_cut.fna NW_007931113.1:" pos "-" pos+49 " NW_007931113.1:" pos+100 "-" pos+149)}}' | awk '(1==NR%2){print $0}(0==NR%2){system("seqtk seq drosophila-GCF_000001215.4_cut.fna | grep -i -o " $0 " | wc -l")}' | awk '(1==NR%4){head1=$0}(2==NR%4){num1=$0}(3==NR%4){head2=$0}(0==NR%4 && 1==num1 && 1==$0){print head1, head2, num1, $0}'
		// bowtie2 -X 2000 -x drosophila-GCF_000001215.4_cut.fna -1 drosophila-coverageTest-R1.fq -2 drosophila-coverageTest-R2.fq | samtools sort -m 10G -@ 4 -T _tmp -o drosophila-coverageTest.bam -
		LoadStats("drosophila-coverageTest.bam", true);
		TestCoverage();
	}

	TEST_F(DataStatsTest, Adapter){
		LoadReference("ecoli-GCF_000005845.2_ASM584v2_genomic.fa");
		CreateTestObject(&species_reference_);

		// Check if adapters are counted properly
		// samtools view -q 10  ../mtest/ecoli-SRR490124-4x.bam | awk '(int($2%32/16) && $8>$4 && $8<$4+200){len=0; for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e || "D"==e){len+=num}; num=0}}; print $4-$8+len, $1}' | sort -n | head
		// zcat ../../../ecoli/data/SRR490124_1.fastq.gz | awk '{if(1==NR%4){if($1=="@SRR490124.10672007" || $1=="@SRR490124.2703087" || $1=="@SRR490124.7718805" || $1=="@SRR490124.7353367" || $1=="@SRR490124.3382600" || $1=="@SRR490124.8368120" || $1=="@SRR490124.10415918" || $1=="@SRR490124.110795" || $1=="@SRR490124.202918"){print $0;pri=1}else{pri=0}}else{if(pri){print $0}}}' > ecoli-SRR490124-adapter-R1.fq
		// Manually added additional information in id and inserted _ to keep it after mapping
		// bowtie2 -X 2000 -x ecoli-GCF_000005845.2_ASM584v2_genomic.fa -1 ecoli-SRR490124-adapter-R1.fq -2 ecoli-SRR490124-adapter-R2.fq | samtools sort -m 10G -@ 4 -T _tmp -o ecoli-SRR490124-adapter.bam -
		LoadStats("ecoli-SRR490124-adapter.bam");

		// Check if adapter information is stored properly
		string save_test_file = string(PROJECT_SOURCE_DIR)+"/test/saveTest.reseq";
		ASSERT_TRUE( test_->Save(save_test_file.c_str()) );

		TestAdapters();

		DeleteTestObject();
		CreateTestObject(&species_reference_);

		ASSERT_TRUE( test_->Load(save_test_file.c_str()) );
		test_->PrepareTesting();
		TestAdapters();
		EXPECT_EQ( 0, remove(save_test_file.c_str()) ) << "Error deleting file: " << save_test_file << '\n';
	}

	TEST_F(DataStatsTest, NexteraAdapter){
		LoadReference("ecoli-GCF_000005845.2_ASM584v2_genomic.fa");
		CreateTestObject(&species_reference_);

		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R1_001.fastq.gz | awk '(2==NR%4){print NR, $0}' | grep --color=always 'CCGAGCCCACGAGAC' | head -n20
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R1_001.fastq.gz | awk '(NR>=25)' | head -n4 > ecoli-S1L001-adapters-R1.fq
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R1_001.fastq.gz | awk '(NR>=33)' | head -n4 >> ecoli-S1L001-adapters-R1.fq
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R1_001.fastq.gz | awk '(NR>=129)' | head -n4 >> ecoli-S1L001-adapters-R1.fq
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R2_001.fastq.gz | awk '(NR>=25)' | head -n4 > ecoli-S1L001-adapters-R2.fq
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R2_001.fastq.gz | awk '(NR>=33)' | head -n4 >> ecoli-S1L001-adapters-R2.fq
		// zcat Nextera-Repeat-600pM-2x151-DI/Nextera_L1/Ecoli1_L001_S5_L001_R2_001.fastq.gz | awk '(NR>=129)' | head -n4 >> ecoli-S1L001-adapters-R2.fq
		// bowtie2 -X 2000 -x ecoli-GCF_000005845.2_ASM584v2_genomic.fa -1 ecoli-S1L001-adapters-R1.fq -2 ecoli-S1L001-adapters-R2.fq | samtools sort -m 10G -@ 4 -T _tmp -o ecoli-S1L001-adapters.bam -

		LoadStats("ecoli-S1L001-adapters.bam",false,false,"Nextera_XTv2");

		TestNexteraAdapters();
	}
}
