#include "CoverageStatsTest.h"
using reseq::CoverageStatsTest;

void CoverageStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new CoverageStats() ) << "Could not allocate memory for CoverageStats object\n";
}

void CoverageStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void CoverageStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void CoverageStatsTest::TestSrr490124Equality(const CoverageStats &test, const char *context){
	EXPECT_EQ(0, test.tmp_errors_forward_.size()) << "SRR490124-4pairs tmp_errors_forward_ not empty for " << context << '\n';
	EXPECT_EQ(0, test.tmp_errors_reverse_.size()) << "SRR490124-4pairs tmp_errors_reverse_ not empty for " << context << '\n';

	// echo "ref prev last5 dist gc err pos seq";samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{pos=0;num=0;for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{pos+=num+1; num=0; print int($2%32/16), pos+$4-1, substr($10,pos,1)}}}' | sort -k1,1 -k2,2nr | awk 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};print $0, start}' |sort -k1,1 -k2,2n | awk 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2}; if(0==$1){dist=$2-start;comp=$3}else{dist=$4-$2;comp="N"; if("A"==$3){comp="T"}; if("C"==$3){comp="G"}; if("G"==$3){comp="C"}; if("T"==$3){comp="A"}}; print $1, $2, int((dist+9)/10), comp; if(0==$1){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2-50 "-" $2 " | seqtk seq")}else{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2 "-" $2+50 " | seqtk seq -r")}}' | awk '(1==NR%3){store=$0}(0==NR%3){print store, substr($0,1,length($0)-1), substr($0,length($0),1)}' | awk '{print $6, substr($5,length($5),1), substr($5,length($5)-4,5), $3, (gsub("G","",$5)+gsub("C","",$5))*2, $4, $2, $1}' | sort
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[1][1]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][1]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(3, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][2]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(1).at(1)[4][1]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[1][3]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[3][0]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(3)[4][0]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[2][1]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[3][1]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(2).at(2)[0][3]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(1).at(3)[5][0]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.dominant_errors_by_distance_.at(3).at(3).at(3)[3][2]) << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';

	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(0).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0))) << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[50][2]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][1]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][2]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[58][2]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[60][1]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(1).at(1)[50][1]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(1).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(1).at(0).at(0))) << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[50][0]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[58][3]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(2).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(2).at(0).at(0))) << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[48][1]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[56][3]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[60][1]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(3).at(1).at(3)), SumVect(test.dominant_errors_by_distance_.at(3).at(1).at(3))) << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context << '\n';
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(3).at(1).at(3)[40][0]) << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';

	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[0]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[0])) << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[1]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[1])) << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[3]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[3])) << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to()) << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';

	// echo "ref prev last5 dist gc err pos seq";samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{pos=0;num=0;for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{pos+=num+1; num=0; print int($2%32/16), pos+$4-1, substr($10,pos,1)}}}' | sort -k1,1 -k2,2nr | awk 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};print $0, start}' |sort -k1,1 -k2,2n | awk 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2}; if(0==$1){dist=$2-start;comp=$3}else{dist=$4-$2;comp="N"; if("A"==$3){comp="T"}; if("C"==$3){comp="G"}; if("G"==$3){comp="C"}; if("T"==$3){comp="A"}}; print $1, $2, int((dist+9)/10), comp; if(0==$1){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2-50 "-" $2 " | seqtk seq")}else{system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2 "-" $2+50 " | seqtk seq -r")}}' | awk '(1==NR%3){store=$0}(0==NR%3){print store, substr($0,1,length($0)-1), substr($0,length($0),1)}' | awk '{print $6, substr($5,length($5),1), substr($5,length($5)-4,5), $3, (gsub("G","",$5)+gsub("C","",$5))*2, $4, $2, $1}' | sort -k1,1 -k6,6 -k5,5n
	EXPECT_EQ(6, SumVect(test.error_rates_by_distance_.at(0).at(1))) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[0][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.error_rates_by_distance_.at(0).at(1)[1][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[2][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[3][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[4][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(1).at(3))) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[1][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[4][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(4, SumVect(test.error_rates_by_distance_.at(2).at(1))) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[1][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.error_rates_by_distance_.at(2).at(1)[2][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[3][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(3, SumVect(test.error_rates_by_distance_.at(3).at(0))) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[0][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[1][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[5][100]) << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';

	EXPECT_EQ(6, SumVect(test.error_rates_by_gc_.at(0).at(1))) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[50][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[52][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[56][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[60][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[62][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[64][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(1).at(3))) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[46][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[58][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(4, SumVect(test.error_rates_by_gc_.at(2).at(1))) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[48][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[50][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[58][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[60][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(3, SumVect(test.error_rates_by_gc_.at(3).at(0))) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[40][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[46][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[54][100]) << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';

	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[0]), SumVect(test.error_rates_by_distance_.at(0).at(1)[0])) << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[1]), SumVect(test.error_rates_by_distance_.at(0).at(1)[1])) << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[2]), SumVect(test.error_rates_by_distance_.at(0).at(1)[2])) << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[3]), SumVect(test.error_rates_by_distance_.at(0).at(1)[3])) << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].from() >= test.error_rates_by_gc_.at(0).at(1).from()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].from() >= test.error_rates_by_gc_.at(0).at(1).from()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[2].from() >= test.error_rates_by_gc_.at(0).at(1).from()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].from() >= test.error_rates_by_gc_.at(0).at(1).from()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].to() <= test.error_rates_by_gc_.at(0).at(1).to()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].to() <= test.error_rates_by_gc_.at(0).at(1).to()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[2].to() <= test.error_rates_by_gc_.at(0).at(1).to()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].to() <= test.error_rates_by_gc_.at(0).at(1).to()) << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n';

	EXPECT_EQ(9, test.error_rates_by_distance_sum_[1][100]) << "SRR490124-4pairs error_rates_by_distance_sum_ wrong for " << context << '\n';
	EXPECT_EQ(400, SumVect(test.error_rates_by_distance_sum_)) << "SRR490124-4pairs error_rates_by_distance_sum_ wrong for " << context << '\n';
	EXPECT_EQ(3, test.error_rates_by_gc_sum_[54][100]) << "SRR490124-4pairs error_rates_by_gc_sum_ wrong for " << context << '\n';
	EXPECT_EQ(400, SumVect(test.error_rates_by_gc_sum_)) << "SRR490124-4pairs error_rates_by_gc_sum_ wrong for " << context << '\n';

	// seqtk seq ecoli-GCF_000005845.2_ASM584v2_genomic.fa | awk '(0==NR%2){print length($0)-100}'
	// 50 first and last bases are ignored due to issues with wrong mappings as real mappings would be partially outside of contig
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{if(1==NR%2){store=$4}else{print store, $4, 100-$4+store}}'
	TestVectEquality({0,{4641242,220,90}}, test.coverage_, context, "SRR490124-4pairs coverage_", " not correct for ");
	for(int strand=2; strand--; ){
		TestVectEquality({0,{4641352,200}}, test.coverage_stranded_.at(strand), context, "SRR490124-4pairs coverage_stranded_", " not correct for ");
		EXPECT_EQ(101, test.coverage_stranded_percent_.at(strand).size() ) << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "].size() wrong for " << context << '\n';
		EXPECT_EQ(110, test.coverage_stranded_percent_.at(strand)[0] ) << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][0] wrong for " << context << '\n';
		EXPECT_EQ(90, test.coverage_stranded_percent_.at(strand)[50] ) << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][50] wrong for " << context << '\n';
		EXPECT_EQ(110, test.coverage_stranded_percent_.at(strand)[100] ) << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][100] wrong for " << context << '\n';
		EXPECT_EQ(0, test.coverage_stranded_percent_min_cov_10_.at(strand).size() ) << "SRR490124-4pairs coverage_stranded_percent_min_cov_10_[" << strand << "].size() wrong for " << context << '\n';
		EXPECT_EQ(0, test.coverage_stranded_percent_min_cov_20_.at(strand).size() ) << "SRR490124-4pairs coverage_stranded_percent_min_cov_20_[" << strand << "].size() wrong for " << context << '\n';
	}
	// samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{pos=$4;num=0;strand=int($2%32/16);for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{pos+=num+1; num=0; print strand, pos-1}}; print strand, $4, "start read", NR; print strand, $4+100, "end read", NR}' | sort -k2,2n | awk 'BEGIN{count0=0;count1=0}{if("" != $3){if(0 != count0 || 0 != count1){print count0, count1, count0+count1, "errors"}; print $0; count0=0; count1=0}else{if(0==$1){count0 += 1}else{count1 += 1}}}'
	TestVectEquality({0,{263,46,1}}, test.error_coverage_, context, "SRR490124-4pairs error_coverage_", " not correct for ");
	EXPECT_EQ(101, test.error_coverage_percent_.size() ) << "SRR490124-4pairs error_coverage_percent_.size() wrong for " << context << '\n';
	EXPECT_EQ(263, test.error_coverage_percent_[0] ) << "SRR490124-4pairs error_coverage_percent_[0] wrong for " << context << '\n';
	EXPECT_EQ(34, test.error_coverage_percent_[50] ) << "SRR490124-4pairs error_coverage_percent_[50] wrong for " << context << '\n';
	EXPECT_EQ(13, test.error_coverage_percent_[100] ) << "SRR490124-4pairs error_coverage_percent_[100] wrong for " << context << '\n';
	EXPECT_EQ(0, test.error_coverage_percent_min_cov_10_.size() ) << "SRR490124-4pairs error_coverage_percent_min_cov_10_.size() wrong for " << context << '\n';
	EXPECT_EQ(0, test.error_coverage_percent_min_cov_20_.size() ) << "SRR490124-4pairs error_coverage_percent_min_cov_20_.size() wrong for " << context << '\n';
	EXPECT_EQ(55, test.error_coverage_percent_stranded_[0][0] ) << "SRR490124-4pairs error_coverage_percent_stranded_[0][0] wrong for " << context << '\n';
	EXPECT_EQ(19, test.error_coverage_percent_stranded_[100][0] ) << "SRR490124-4pairs error_coverage_percent_stranded_[100][0] wrong for " << context << '\n';
	EXPECT_EQ(15, test.error_coverage_percent_stranded_[0][100] ) << "SRR490124-4pairs error_coverage_percent_stranded_[0][100] wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_coverage_percent_stranded_[100][100] ) << "SRR490124-4pairs error_coverage_percent_stranded_[100][100] wrong for " << context << '\n';
	EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_10_.size() ) << "SRR490124-4pairs error_coverage_percent_stranded_min_strand_cov_10_.size() wrong for " << context << '\n';
	EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_20_.size() ) << "SRR490124-4pairs error_coverage_percent_stranded_min_strand_cov_20_.size() wrong for " << context << '\n';

	EXPECT_EQ(NULL, test.first_block_ ) << "SRR490124-4pairs first_block_ wrong for " << context << '\n';
	EXPECT_TRUE(NULL == test.last_block_ ) << "SRR490124-4pairs last_block_ wrong for " << context << '\n';
	EXPECT_EQ(0, test.reusable_blocks_.size() ) << "SRR490124-4pairs reusable_blocks_ wrong for " << context << '\n';
}

void CoverageStatsTest::TestDuplicates(const CoverageStats &test){
	// seqtk seq ecoli-GCF_000005845.2_ASM584v2_genomic.fa | awk '(0==NR%2){print length($0)-100}'
	// samtools view -q 10 -f 3 ecoli-duplicates.bam | awk '(0 != substr($1,12,1)){pos=$4; for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e){while(0<num){print pos; ++pos; --num}};if("D"==e){pos+=num};num=0}}}' | sort -n | uniq -c | awk '{print $1}' | sort -n | uniq -c | awk 'BEGIN{sum=0}{print $0; sum+=$1}END{print 4641552-sum, 0}'
	TestVectEquality({0,{4641252,1,101,98,24,1,0,0,25,3,22,0,0,2,21,0,2}}, test.coverage_, "duplicates test", "coverage_", " not correct for ");

	EXPECT_EQ(NULL, test.first_block_ ) << "SRR490124-4pairs first_block_ wrong in duplicates test\n";
	EXPECT_TRUE(NULL == test.last_block_ ) << "SRR490124-4pairs last_block_ wrong in duplicates test\n";
	EXPECT_EQ(0, test.reusable_blocks_.size() ) << "SRR490124-4pairs reusable_blocks_ wrong in duplicates test\n";
}

void CoverageStatsTest::TestVariants(const CoverageStats &test){
	// Manually modified values from TestSrr490124Equality
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[1][1]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][1]);
	EXPECT_EQ(3, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][2]);
	EXPECT_EQ(0, test.dominant_errors_by_distance_.at(0).at(1).at(1)[4][1]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[1][3]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[3][0]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(3)[4][0]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[2][1]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[3][1]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(2).at(2)[0][3]);
	EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(1).at(3)[5][0]);
	EXPECT_EQ(2, test.dominant_errors_by_distance_.at(3).at(3).at(3)[3][2]);

	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(0).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)));
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[50][2]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][1]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][2]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[58][2]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[60][1]);
	EXPECT_EQ(0, test.dominant_errors_by_gc_.at(0).at(1).at(1)[50][1]);
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(1).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(1).at(0).at(0)));
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[50][0]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[58][3]);
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(2).at(0).at(0)), SumVect(test.dominant_errors_by_distance_.at(2).at(0).at(0)));
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[48][1]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[56][3]);
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[60][1]);
	EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(3).at(1).at(3)), SumVect(test.dominant_errors_by_distance_.at(3).at(1).at(3)));
	EXPECT_EQ(1, test.dominant_errors_by_gc_.at(3).at(1).at(3)[40][0]);

	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[0]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[0]));
	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[1]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[1]));
	EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[3]), SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[3]));
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].from() >= test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to());
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to());
	EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to());

	EXPECT_EQ(5, SumVect(test.error_rates_by_distance_.at(0).at(1)));
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[0][100]);
	EXPECT_EQ(2, test.error_rates_by_distance_.at(0).at(1)[1][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[2][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[3][100]);
	EXPECT_EQ(0, test.error_rates_by_distance_.at(0).at(1)[4][100]);
	EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(1).at(3)));
	EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[1][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[4][100]);
	EXPECT_EQ(4, SumVect(test.error_rates_by_distance_.at(2).at(1)));
	EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[1][100]);
	EXPECT_EQ(2, test.error_rates_by_distance_.at(2).at(1)[2][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[3][100]);
	EXPECT_EQ(3, SumVect(test.error_rates_by_distance_.at(3).at(0)));
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[0][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[1][100]);
	EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[5][100]);

	EXPECT_EQ(5, SumVect(test.error_rates_by_gc_.at(0).at(1)));
	EXPECT_EQ(0, test.error_rates_by_gc_.at(0).at(1)[50][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[52][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[56][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[60][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[62][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[64][100]);
	EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(1).at(3)));
	EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[46][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[58][100]);
	EXPECT_EQ(4, SumVect(test.error_rates_by_gc_.at(2).at(1)));
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[48][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[50][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[58][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[60][100]);
	EXPECT_EQ(3, SumVect(test.error_rates_by_gc_.at(3).at(0)));
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[40][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[46][100]);
	EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[54][100]);

	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[0]), SumVect(test.error_rates_by_distance_.at(0).at(1)[0]));
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[1]), SumVect(test.error_rates_by_distance_.at(0).at(1)[1]));
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[2]), SumVect(test.error_rates_by_distance_.at(0).at(1)[2]));
	EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[3]), SumVect(test.error_rates_by_distance_.at(0).at(1)[3]));
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].from() >= test.error_rates_by_gc_.at(0).at(1).from());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].from() >= test.error_rates_by_gc_.at(0).at(1).from());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[2].from() >= test.error_rates_by_gc_.at(0).at(1).from());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].from() >= test.error_rates_by_gc_.at(0).at(1).from());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].to() <= test.error_rates_by_gc_.at(0).at(1).to());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].to() <= test.error_rates_by_gc_.at(0).at(1).to());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[2].to() <= test.error_rates_by_gc_.at(0).at(1).to());
	EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].to() <= test.error_rates_by_gc_.at(0).at(1).to());

	EXPECT_EQ(9, test.error_rates_by_distance_sum_[1][100]);
	EXPECT_EQ(394, SumVect(test.error_rates_by_distance_sum_));
	EXPECT_EQ(3, test.error_rates_by_gc_sum_[54][100]);
	EXPECT_EQ(394, SumVect(test.error_rates_by_gc_sum_));

	TestVectEquality({0,{260,45,1}}, test.error_coverage_, "variant test", "SRR490124-4pairs error_coverage_", " not correct for ");
	EXPECT_EQ(101, test.error_coverage_percent_.size() );
	EXPECT_EQ(260, test.error_coverage_percent_[0] );
	EXPECT_EQ(33, test.error_coverage_percent_[50] );
	EXPECT_EQ(13, test.error_coverage_percent_[100] );
	EXPECT_EQ(0, test.error_coverage_percent_min_cov_10_.size() );
	EXPECT_EQ(0, test.error_coverage_percent_min_cov_20_.size() );
	EXPECT_EQ(54, test.error_coverage_percent_stranded_[0][0] );
	EXPECT_EQ(18, test.error_coverage_percent_stranded_[100][0] );
	EXPECT_EQ(15, test.error_coverage_percent_stranded_[0][100] );
	EXPECT_EQ(1, test.error_coverage_percent_stranded_[100][100] );
	EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_10_.size() );
	EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_20_.size() );
}

void CoverageStatsTest::TestCrossDuplicates(const CoverageStats &test){
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk 'BEGIN{sum=0}(0==NR%2 && length($0)>2100){sum += length($0)-100}END{print sum}'
	TestVectEquality({0,{3780545}}, test.coverage_, "cross duplicates test", "coverage_", " not correct for ");

	EXPECT_EQ(NULL, test.first_block_ ) << "SRR490124-4pairs first_block_ wrong in cross duplicates test\n";
	EXPECT_TRUE(NULL == test.last_block_ ) << "SRR490124-4pairs last_block_ wrong in cross duplicates test\n";
	EXPECT_EQ(0, test.reusable_blocks_.size() ) << "SRR490124-4pairs reusable_blocks_ wrong in cross duplicates test\n";
}

void CoverageStatsTest::TestCoverage(const CoverageStats &test){
	// seqtk seq drosophila-GCF_000001215.4_cut.fna | awk 'BEGIN{sum=0}(0==NR%2 && length($0)>2100){sum += length($0)-100}END{print sum-300}'
	TestVectEquality({0,{3780245, 300}}, test.coverage_, "coverage test", "coverage_", " not correct for ");

	EXPECT_EQ(NULL, test.first_block_ ) << "SRR490124-4pairs first_block_ wrong in coverage test\n";
	EXPECT_TRUE(NULL == test.last_block_ ) << "SRR490124-4pairs last_block_ wrong in coverage test\n";
	EXPECT_EQ(0, test.reusable_blocks_.size() ) << "SRR490124-4pairs reusable_blocks_ wrong in coverage test\n";
}
