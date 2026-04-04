#include "CoverageStatsTest.h"
using reseq::CoverageStatsTest;

#include <string>
#include <thread>
#include <vector>
using std::string;

void CoverageStatsTest::CreateTestObject() {
    test_ = std::make_unique<CoverageStats>();
    ASSERT_TRUE(test_) << "Could not allocate memory for CoverageStats object\n";
}

void CoverageStatsTest::DeleteTestObject() {
    test_.reset();
}

void CoverageStatsTest::TearDown() {
    BasicTestClass::TearDown();
    DeleteTestObject();
}

void CoverageStatsTest::TestNonSystematicErrorRate() {
    test_->tmp_block_error_rate_.resize(101);
    test_->tmp_block_percent_systematic_.resize(101);
    test_->tmp_systematic_error_p_values_.resize(101);

    // samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:101-110
    // TAAAATTTTA
    CoverageStats::CoverageBlock block(0, 100);
    block.coverage_.resize(10);
    block.coverage_.at(0).coverage_forward_.at(4) = 3;
    block.coverage_.at(0).coverage_forward_.at(1) = 1;
    block.coverage_.at(0).coverage_forward_.at(3) = 23;
    block.coverage_.at(0).coverage_reverse_.at(0) = 23;
    block.coverage_.at(2).coverage_forward_.at(0) = 22;
    block.coverage_.at(2).coverage_reverse_.at(3) = 22;
    block.coverage_.at(5).coverage_forward_.at(3) = 22;
    block.coverage_.at(5).coverage_reverse_.at(0) = 22;

    block.coverage_.at(7).valid_ = false;
    block.coverage_.at(7).coverage_forward_.at(4) = 100;

    CoverageStats::ThreadData thread;
    // R: 1-ppois(2, 27*4/(88+46+4)/3)
    EXPECT_NEAR(0.002436186, test_->GetPositionProbabilities(&block, species_reference_, thread), 0.000000001);
    EXPECT_EQ(1, test_->tmp_block_error_rate_.at(3));
    EXPECT_EQ(1, test_->tmp_block_percent_systematic_.at(17));
    EXPECT_EQ(5, test_->tmp_systematic_error_p_values_.at(100));
    EXPECT_EQ(1, test_->tmp_systematic_error_p_values_.at(0));
}

void CoverageStatsTest::TestSrr490124Equality(const CoverageStats& test, const char* context) {
    // echo "count ref prev last5 dist gc err pos seq";samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk
    // '{pos=0;num=0;for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{pos+=num+1;
    // num=0; print int($2%32/16), pos+$4-1, substr($10,pos,1)}}}' | sort -k1,1 -k2,2nr | awk
    // 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};print $0, start}' |sort -k1,1 -k2,2n | awk
    // 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};
    // if(0==$1){dist=$2-start;comp=$3}else{dist=$4-$2;comp="N"; if("A"==$3){comp="T"}; if("C"==$3){comp="G"};
    // if("G"==$3){comp="C"}; if("T"==$3){comp="A"}}; print $1, $2, int((dist+9)/10), comp; if(0==$1){system("samtools
    // faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2-50 "-" $2 " | seqtk seq")}else{system("samtools
    // faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2 "-" $2+50 " | seqtk seq -r")}}' | awk
    // '(1==NR%3){store=$0}(0==NR%3){print store, substr($0,1,length($0)-1), substr($0,length($0),1)}' | awk '{print $6,
    // substr($5,length($5),1), substr($5,length($5)-4,5), $3, (gsub("G","",$5)+gsub("C","",$5))*2, $4, $2, $1}' | sort
    // | uniq -c
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[1][1])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][1])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(2, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][2])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(1).at(1)[0][1])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[1][3])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(1).at(0).at(0)[3][0])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(3)[4][0])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[2][1])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(2).at(0).at(0)[3][1])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(2).at(2).at(2)[0][3])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(1).at(3)[5][0])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(3).at(3)[3][2])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(2, test.dominant_errors_by_distance_.at(3).at(3).at(3)[5][4])
        << "SRR490124-4pairs dominant_errors_by_distance_ wrong for " << context << '\n';

    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(0).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)))
        << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(0).at(0).at(0)[50][2])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][1])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][2])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[58][2])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[60][1])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(1).at(1)[52][1])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(1).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(1).at(0).at(0)))
        << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(1).at(0).at(0)[50][0])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[58][3])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(2).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(2).at(0).at(0)))
        << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(2).at(0).at(0)[48][1])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(2).at(0).at(0)[56][3])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[60][1])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(3).at(1).at(3)),
              SumVect(test.dominant_errors_by_distance_.at(3).at(1).at(3)))
        << "SRR490124-4pairs dominant_errors_by_distance_ and dominant_errors_by_gc_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(3).at(1).at(3)[40][0])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(2, test.dominant_errors_by_gc_.at(3).at(3).at(3)[42][4])
        << "SRR490124-4pairs dominant_errors_by_gc_ wrong for " << context << '\n';

    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[0]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[0]))
        << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[1]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[1]))
        << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context
        << '\n';
    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[3]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[3]))
        << "SRR490124-4pairs dominant_errors_by_distance_ and gc_by_distance_de_ are inconsistent for " << context
        << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].to() <= test.dominant_errors_by_gc_.at(0).at(0).at(0).to())
        << "SRR490124-4pairs dominant_errors_by_gc_ and gc_by_distance_de_ are inconsistent for " << context << '\n';

    EXPECT_EQ(1, test.dominant_errors_by_start_rates_.at(0).at(1).at(1)[0][1]);
    EXPECT_EQ(2, test.start_rates_by_distance_de_.at(3).at(2).at(2)[1][100]);

    // echo "count ref prev last5 dist gc err pos seq";samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk
    // '{pos=0;num=0;for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{pos+=num+1;
    // num=0; print int($2%32/16), pos+$4-1, substr($10,pos,1)}}}' | sort -k1,1 -k2,2nr | awk
    // 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};print $0, start}' |sort -k1,1 -k2,2n | awk
    // 'BEGIN{start=0}{if(start+100<$2 || start-100>$2){start=$2};
    // if(0==$1){dist=$2-start;comp=$3}else{dist=$4-$2;comp="N"; if("A"==$3){comp="T"}; if("C"==$3){comp="G"};
    // if("G"==$3){comp="C"}; if("T"==$3){comp="A"}}; print $1, $2, int((dist+9)/10), comp; if(0==$1){system("samtools
    // faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2-50 "-" $2 " | seqtk seq")}else{system("samtools
    // faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:" $2 "-" $2+50 " | seqtk seq -r")}}' | awk
    // '(1==NR%3){store=$0}(0==NR%3){print store, substr($0,1,length($0)-1), substr($0,length($0),1)}' | awk '{print $6,
    // substr($5,length($5),1), substr($5,length($5)-4,5), $3, (gsub("G","",$5)+gsub("C","",$5))*2, $4, $2, $1}' | sort
    // -k1,1 -k6,6 -k5,5n | uniq -c
    EXPECT_EQ(4, SumVect(test.error_rates_by_distance_.at(0).at(1)))
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[0][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(2, test.error_rates_by_distance_.at(0).at(1)[1][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_distance_.at(0).at(1)[2][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[3][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_distance_.at(0).at(1)[4][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(1).at(3)))
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[1][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[4][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, SumVect(test.error_rates_by_distance_.at(2).at(1)))
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_distance_.at(2).at(1)[1][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[2][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_distance_.at(2).at(1)[3][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(3).at(0)))
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_distance_.at(3).at(0)[0][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[1][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[5][100])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';
    EXPECT_EQ(4, test.error_rates_by_distance_.at(3).at(4)[5][0])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';

    EXPECT_EQ(4, SumVect(test.error_rates_by_gc_.at(0).at(1)))
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(0).at(1)[50][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[52][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[56][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[60][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[62][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(0).at(1)[64][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(1).at(3)))
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[46][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[58][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, SumVect(test.error_rates_by_gc_.at(2).at(1)))
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[48][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[50][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[58][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[60][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(3).at(0)))
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[40][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_rates_by_gc_.at(3).at(0)[46][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[54][100])
        << "SRR490124-4pairs error_rates_by_gc_ wrong for " << context << '\n';
    EXPECT_EQ(4, test.error_rates_by_gc_.at(3).at(4)[42][0])
        << "SRR490124-4pairs error_rates_by_distance_ wrong for " << context << '\n';

    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[0]), SumVect(test.error_rates_by_distance_.at(0).at(1)[0]))
        << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[1]), SumVect(test.error_rates_by_distance_.at(0).at(1)[1]))
        << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[2]), SumVect(test.error_rates_by_distance_.at(0).at(1)[2]))
        << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[3]), SumVect(test.error_rates_by_distance_.at(0).at(1)[3]))
        << "SRR490124-4pairs error_rates_by_distance_ and gc_by_distance_er_ are inconsistent for " << context << '\n';
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[0].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[1].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(2)[2].from() >= test.error_rates_by_gc_.at(0).at(2).from())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(2)[2].from() << " < " << test.error_rates_by_gc_.at(0).at(2).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[3].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[0].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[1].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(2)[2].to() <= test.error_rates_by_gc_.at(0).at(2).to())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(2)[2].to() << " > " << test.error_rates_by_gc_.at(0).at(2).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << "SRR490124-4pairs error_rates_by_gc_ and gc_by_distance_er_ are inconsistent for " << context << '\n'
        << test.gc_by_distance_er_.at(0).at(1)[3].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();

    EXPECT_EQ(4, test.error_rates_by_distance_sum_[1][100])
        << "SRR490124-4pairs error_rates_by_distance_sum_ wrong for " << context << '\n';
    EXPECT_EQ(398, SumVect(test.error_rates_by_distance_sum_))
        << "SRR490124-4pairs error_rates_by_distance_sum_ wrong for " << context << '\n';
    EXPECT_EQ(2, test.error_rates_by_gc_sum_[54][100])
        << "SRR490124-4pairs error_rates_by_gc_sum_ wrong for " << context << '\n';
    EXPECT_EQ(398, SumVect(test.error_rates_by_gc_sum_))
        << "SRR490124-4pairs error_rates_by_gc_sum_ wrong for " << context << '\n';

    EXPECT_EQ(1, test.start_rates_by_distance_er_.at(0).at(1)[0][0]);
    EXPECT_EQ(1, test.start_rates_by_distance_er_.at(3).at(0)[1][100]);

    // seqtk seq ecoli-GCF_000005845.2_ASM584v2_genomic.fa | awk '(0==NR%2){print length($0)-100}'
    // 50 first and last bases are ignored due to issues with wrong mappings as real mappings would be partially outside
    // of contig samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(NR <= 3 || NR ==
    // 12){if(1==NR%2){store=$4}else{print store, $4, 100-$4+store}}'
    TestVectEquality({0, {4641242, 200 - 2 * 74, 74, 0, 0, 0, 0, 0, 200 - 2 * 16, 0, 0, 0, 0, 0, 0, 0, 16}},
                     test.coverage_, context, "SRR490124-4pairs coverage_", " not correct for ");
    for (int strand = 2; strand--;) {
        TestVectEquality({0, {4641352, 100, 0, 0, 0, 0, 0, 0, 100}}, test.coverage_stranded_.at(strand), context,
                         "SRR490124-4pairs coverage_stranded_", " not correct for ");
        EXPECT_EQ(101, test.coverage_stranded_percent_.at(strand).size())
            << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "].size() wrong for " << context << '\n';
        EXPECT_EQ(110, test.coverage_stranded_percent_.at(strand)[0])
            << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][0] wrong for " << context << '\n';
        EXPECT_EQ(90, test.coverage_stranded_percent_.at(strand)[50])
            << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][50] wrong for " << context << '\n';
        EXPECT_EQ(110, test.coverage_stranded_percent_.at(strand)[100])
            << "SRR490124-4pairs coverage_stranded_percent_[" << strand << "][100] wrong for " << context << '\n';
        EXPECT_EQ(1, test.coverage_stranded_percent_min_cov_10_.at(strand).size())
            << "SRR490124-4pairs coverage_stranded_percent_min_cov_10_[" << strand << "].size() wrong for " << context
            << '\n';
        EXPECT_EQ(16, test.coverage_stranded_percent_min_cov_10_.at(strand)[50])
            << "SRR490124-4pairs coverage_stranded_percent_min_cov_10_[" << strand << "][50] wrong for " << context
            << '\n';
        EXPECT_EQ(0, test.coverage_stranded_percent_min_cov_20_.at(strand).size())
            << "SRR490124-4pairs coverage_stranded_percent_min_cov_20_[" << strand << "].size() wrong for " << context
            << '\n';
    }
    // samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk
    // '{pos=$4;num=0;strand=int($2%32/16);for(i=6;i<=length($18);i+=1){b=substr($18,i,1);if(b ~
    // /^[0-9]/){num=num*10+b}else{pos+=num+1; num=0; print strand, pos-1}}; print strand, $4, "start read", NR; print
    // strand, $4+100, "end read", NR}' | sort -k2,2n | awk 'BEGIN{count0=0;count1=0}{if("" != $3){if(0 != count0 || 0
    // != count1){print count0, count1, count0+count1, "errors"}; print $0; count0=0; count1=0}else{if(0==$1){count0 +=
    // 1}else{count1 += 1}}}'
    //                                                                   GC      GTA   A      TAT  AC AGCA(15)*1
    // TTATCCGGATCAGGTTCGACGGGTATTTTCTCAGCGCACGCGTACGCGTGGCACCCCGTTGAGAACGTGGTTAGTTGGGTGCTTTTGTGTGCAGACGCAC
    //                          TTTTCACTGGAGCGCGTAGGCGTGGCTCCCCGTTGAGAACGGCGTTAGTGTAGTGATTTTGTTATCAACCAGCAATCATGGATCCGGTGGCGCAAACCAC
    //                    1*(10)   CTCAGC CA      C       A
    //
    //                                                  A     T         T        GTT       ACCC     T     (11)*8
    // GAGCAAACAATCACAGCATGTATTAATTGCCCTGCCCCACCCGCTGCTTCCCCTGGACAGTTTAGGCTTAGTCTCTGGTATCTTTGATATTTTCACGCTT
    //                                                                                    TACCACGCTCTCCGTTGATCTCTCGCAAGGTGGAACCCAACTCGCCCCACTGTGGTTCCCGACGTCCATCATGATGGTGGCGTTTTATCGCCATGCCGGG
    //                                                                              8*(12)    CTTT    GC    G  T     GTT
    //                                                                              C
    TestVectEquality({0, {263, 25, 0, 0, 0, 0, 0, 0, 21}}, test.error_coverage_, context,
                     "SRR490124-4pairs error_coverage_", " not correct for ");
    EXPECT_EQ(101, test.error_coverage_percent_.size())
        << "SRR490124-4pairs error_coverage_percent_.size() wrong for " << context << '\n';
    EXPECT_EQ(263, test.error_coverage_percent_[0])
        << "SRR490124-4pairs error_coverage_percent_[0] wrong for " << context << '\n';
    EXPECT_EQ(34, test.error_coverage_percent_[50])
        << "SRR490124-4pairs error_coverage_percent_[50] wrong for " << context << '\n';
    EXPECT_EQ(12, test.error_coverage_percent_[100])
        << "SRR490124-4pairs error_coverage_percent_[100] wrong for " << context << '\n';
    EXPECT_EQ(51, test.error_coverage_percent_min_cov_10_.size())
        << "SRR490124-4pairs error_coverage_percent_min_cov_10_.size() wrong for " << context << '\n';
    EXPECT_EQ(6, test.error_coverage_percent_min_cov_10_[0])
        << "SRR490124-4pairs error_coverage_percent_min_cov_10_[0] wrong for " << context << '\n';
    EXPECT_EQ(9, test.error_coverage_percent_min_cov_10_[50])
        << "SRR490124-4pairs error_coverage_percent_min_cov_10_[50] wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_coverage_percent_min_cov_10_[100])
        << "SRR490124-4pairs error_coverage_percent_min_cov_10_[100] wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_coverage_percent_min_cov_20_.size())
        << "SRR490124-4pairs error_coverage_percent_min_cov_20_.size() wrong for " << context << '\n';
    EXPECT_EQ(55, test.error_coverage_percent_stranded_[0][0])
        << "SRR490124-4pairs error_coverage_percent_stranded_[0][0] wrong for " << context << '\n';
    EXPECT_EQ(19, test.error_coverage_percent_stranded_[100][0])
        << "SRR490124-4pairs error_coverage_percent_stranded_[100][0] wrong for " << context << '\n';
    EXPECT_EQ(15, test.error_coverage_percent_stranded_[0][100])
        << "SRR490124-4pairs error_coverage_percent_stranded_[0][100] wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_coverage_percent_stranded_[100][100])
        << "SRR490124-4pairs error_coverage_percent_stranded_[100][100] wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_10_.size())
        << "SRR490124-4pairs error_coverage_percent_stranded_min_strand_cov_10_.size() wrong for " << context << '\n';
    EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_20_.size())
        << "SRR490124-4pairs error_coverage_percent_stranded_min_strand_cov_20_.size() wrong for " << context << '\n';

    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load())
        << "SRR490124-4pairs first_live_idx_ wrong for " << context << '\n';
    EXPECT_EQ(SIZE_MAX, test.last_live_idx_.load()) << "SRR490124-4pairs last_live_idx_ wrong for " << context << '\n';
    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load()) << "SRR490124-4pairs blocks cleanup wrong for " << context << '\n';
}

void CoverageStatsTest::TestDuplicates(const CoverageStats& test) {
    // seqtk seq ecoli-GCF_000005845.2_ASM584v2_genomic.fa | awk '(0==NR%2){print length($0)-100}'
    // samtools view -q 10 -f 3 ecoli-duplicates.bam | awk '(0 != substr($1,12,1)){pos=$4;
    // for(i=1;i<=length($6);i+=1){e=substr($6,i,1);if(e ~ /^[0-9]/){num=num*10+e}else{if("M"==e){while(0<num){print
    // pos; ++pos; --num}};if("D"==e){pos+=num};num=0}}}' | sort -n | uniq -c | awk '{print $1}' | sort -n | uniq -c |
    // awk 'BEGIN{sum=0}{print $0; sum+=$1}END{print 4641552-sum, 0}'
    TestVectEquality({0, {4641252, 1, 101, 98, 24, 1, 0, 0, 25, 3, 22, 0, 0, 2, 21, 0, 2}}, test.coverage_,
                     "duplicates test", "coverage_", " not correct for ");

    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load()) << "SRR490124-4pairs first_live_idx_ wrong in duplicates test\n";
    EXPECT_EQ(SIZE_MAX, test.last_live_idx_.load()) << "SRR490124-4pairs last_live_idx_ wrong in duplicates test\n";
    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load()) << "SRR490124-4pairs blocks cleanup wrong in duplicates test\n";
}

void CoverageStatsTest::TestVariants(const CoverageStats& test) {
    // Manually modified values from TestSrr490124Equality
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[1][1]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][1]);
    EXPECT_EQ(2, test.dominant_errors_by_distance_.at(0).at(0).at(0)[3][2]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(0).at(1).at(1)[0][1]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(0)[1][3]);
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(1).at(0).at(0)[3][0]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(1).at(0).at(3)[4][0]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(2).at(0).at(0)[2][1]);
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(2).at(0).at(0)[3][1]);
    EXPECT_EQ(0, test.dominant_errors_by_distance_.at(2).at(2).at(2)[0][3]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(1).at(3)[5][0]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(3).at(3)[3][2]);
    EXPECT_EQ(1, test.dominant_errors_by_distance_.at(3).at(3).at(3)[5][4]);

    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(0).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)));
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(0).at(0).at(0)[50][2]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][1]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[56][2]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[58][2]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(0).at(0)[60][1]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(0).at(1).at(1)[52][1]);
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(1).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(1).at(0).at(0)));
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(1).at(0).at(0)[50][0]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(1).at(0).at(0)[58][3]);
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(2).at(0).at(0)),
              SumVect(test.dominant_errors_by_distance_.at(2).at(0).at(0)));
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(2).at(0).at(0)[48][1]);
    EXPECT_EQ(0, test.dominant_errors_by_gc_.at(2).at(0).at(0)[56][3]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(2).at(0).at(0)[60][1]);
    EXPECT_EQ(SumVect(test.dominant_errors_by_gc_.at(3).at(1).at(3)),
              SumVect(test.dominant_errors_by_distance_.at(3).at(1).at(3)));
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(3).at(1).at(3)[40][0]);
    EXPECT_EQ(1, test.dominant_errors_by_gc_.at(3).at(3).at(3)[42][4]);

    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[0]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[0]));
    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[1]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[1]));
    EXPECT_EQ(SumVect(test.gc_by_distance_de_.at(0).at(0).at(0)[3]),
              SumVect(test.dominant_errors_by_distance_.at(0).at(0).at(0)[3]));
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].from() >=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).from());
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[0].to() <=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).to());
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[1].to() <=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).to());
    EXPECT_TRUE(test.gc_by_distance_de_.at(0).at(0).at(0)[3].to() <=
                test.dominant_errors_by_gc_.at(0).at(0).at(0).to());

    EXPECT_EQ(4, SumVect(test.error_rates_by_distance_.at(0).at(1)));
    EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[0][100]);
    EXPECT_EQ(2, test.error_rates_by_distance_.at(0).at(1)[1][100]);
    EXPECT_EQ(0, test.error_rates_by_distance_.at(0).at(1)[2][100]);
    EXPECT_EQ(1, test.error_rates_by_distance_.at(0).at(1)[3][100]);
    EXPECT_EQ(0, test.error_rates_by_distance_.at(0).at(1)[4][100]);
    EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(1).at(3)));
    EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[1][100]);
    EXPECT_EQ(1, test.error_rates_by_distance_.at(1).at(3)[4][100]);
    EXPECT_EQ(1, SumVect(test.error_rates_by_distance_.at(2).at(1)));
    EXPECT_EQ(0, test.error_rates_by_distance_.at(2).at(1)[1][100]);
    EXPECT_EQ(1, test.error_rates_by_distance_.at(2).at(1)[2][100]);
    EXPECT_EQ(0, test.error_rates_by_distance_.at(2).at(1)[3][100]);
    EXPECT_EQ(2, SumVect(test.error_rates_by_distance_.at(3).at(0)));
    EXPECT_EQ(0, test.error_rates_by_distance_.at(3).at(0)[0][100]);
    EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[1][100]);
    EXPECT_EQ(1, test.error_rates_by_distance_.at(3).at(0)[5][100]);
    EXPECT_EQ(3, test.error_rates_by_distance_.at(3).at(4)[5][0]);

    EXPECT_EQ(4, SumVect(test.error_rates_by_gc_.at(0).at(1)));
    EXPECT_EQ(0, test.error_rates_by_gc_.at(0).at(1)[50][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[52][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[56][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[60][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(0).at(1)[62][100]);
    EXPECT_EQ(0, test.error_rates_by_gc_.at(0).at(1)[64][100]);
    EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(1).at(3)));
    EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[46][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(1).at(3)[58][100]);
    EXPECT_EQ(1, SumVect(test.error_rates_by_gc_.at(2).at(1)));
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[48][100]);
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[50][100]);
    EXPECT_EQ(0, test.error_rates_by_gc_.at(2).at(1)[58][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(2).at(1)[60][100]);
    EXPECT_EQ(2, SumVect(test.error_rates_by_gc_.at(3).at(0)));
    EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[40][100]);
    EXPECT_EQ(0, test.error_rates_by_gc_.at(3).at(0)[46][100]);
    EXPECT_EQ(1, test.error_rates_by_gc_.at(3).at(0)[54][100]);
    EXPECT_EQ(3, test.error_rates_by_gc_.at(3).at(4)[42][0]);

    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[0]), SumVect(test.error_rates_by_distance_.at(0).at(1)[0]));
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[1]), SumVect(test.error_rates_by_distance_.at(0).at(1)[1]));
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[2]), SumVect(test.error_rates_by_distance_.at(0).at(1)[2]));
    EXPECT_EQ(SumVect(test.gc_by_distance_er_.at(0).at(1)[3]), SumVect(test.error_rates_by_distance_.at(0).at(1)[3]));
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << test.gc_by_distance_er_.at(0).at(1)[0].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << test.gc_by_distance_er_.at(0).at(1)[1].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(2)[2].from() >= test.error_rates_by_gc_.at(0).at(2).from())
        << test.gc_by_distance_er_.at(0).at(2)[2].from() << " < " << test.error_rates_by_gc_.at(0).at(2).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].from() >= test.error_rates_by_gc_.at(0).at(1).from())
        << test.gc_by_distance_er_.at(0).at(1)[3].from() << " < " << test.error_rates_by_gc_.at(0).at(1).from();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[0].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << test.gc_by_distance_er_.at(0).at(1)[0].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[1].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << test.gc_by_distance_er_.at(0).at(1)[1].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(2)[2].to() <= test.error_rates_by_gc_.at(0).at(2).to())
        << test.gc_by_distance_er_.at(0).at(2)[2].to() << " > " << test.error_rates_by_gc_.at(0).at(2).to();
    EXPECT_TRUE(test.gc_by_distance_er_.at(0).at(1)[3].to() <= test.error_rates_by_gc_.at(0).at(1).to())
        << test.gc_by_distance_er_.at(0).at(1)[3].to() << " > " << test.error_rates_by_gc_.at(0).at(1).to();

    EXPECT_EQ(4, test.error_rates_by_distance_sum_[1][100]);
    EXPECT_EQ(390, SumVect(test.error_rates_by_distance_sum_));
    EXPECT_EQ(2, test.error_rates_by_gc_sum_[54][100]);
    EXPECT_EQ(390, SumVect(test.error_rates_by_gc_sum_));

    TestVectEquality({0, {259, 24, 0, 0, 0, 0, 0, 0, 21}}, test.error_coverage_, "variant test",
                     "SRR490124-4pairs error_coverage_", " not correct for ");
    EXPECT_EQ(101, test.error_coverage_percent_.size());
    EXPECT_EQ(259, test.error_coverage_percent_[0]);
    EXPECT_EQ(33, test.error_coverage_percent_[50]);
    EXPECT_EQ(12, test.error_coverage_percent_[100]);
    EXPECT_EQ(51, test.error_coverage_percent_min_cov_10_.size());
    EXPECT_EQ(5, test.error_coverage_percent_min_cov_10_[0]);
    EXPECT_EQ(9, test.error_coverage_percent_min_cov_10_[50]);
    EXPECT_EQ(0, test.error_coverage_percent_min_cov_10_[100]);
    EXPECT_EQ(0, test.error_coverage_percent_min_cov_20_.size());
    EXPECT_EQ(53, test.error_coverage_percent_stranded_[0][0]);
    EXPECT_EQ(18, test.error_coverage_percent_stranded_[100][0]);
    EXPECT_EQ(15, test.error_coverage_percent_stranded_[0][100]);
    EXPECT_EQ(0, test.error_coverage_percent_stranded_[100][100]);
    EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_10_.size());
    EXPECT_EQ(0, test.error_coverage_percent_stranded_min_strand_cov_20_.size());
}

void CoverageStatsTest::TestCrossDuplicates(const CoverageStats& test) {
    // seqtk seq drosophila-GCF_000001215.4_cut.fna | awk 'BEGIN{sum=0}(0==NR%2 && length($0)>2100){sum +=
    // length($0)-100}END{print sum}'
    TestVectEquality({0, {2940870}}, test.coverage_, "cross duplicates test", "coverage_", " not correct for ");

    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load())
        << "SRR490124-4pairs first_live_idx_ wrong in cross duplicates test\n";
    EXPECT_EQ(SIZE_MAX, test.last_live_idx_.load())
        << "SRR490124-4pairs last_live_idx_ wrong in cross duplicates test\n";
    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load())
        << "SRR490124-4pairs blocks cleanup wrong in cross duplicates test\n";
}

void CoverageStatsTest::TestCoverage(const CoverageStats& test) {
    // seqtk seq drosophila-GCF_000001215.4_cut.fna | awk 'BEGIN{sum=0}(0==NR%2 && length($0)>2100){sum +=
    // length($0)-100}END{print sum-300}'
    TestVectEquality({0, {2940670, 200}}, test.coverage_, "coverage test", "coverage_", " not correct for ");

    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load()) << "SRR490124-4pairs first_live_idx_ wrong in coverage test\n";
    EXPECT_EQ(SIZE_MAX, test.last_live_idx_.load()) << "SRR490124-4pairs last_live_idx_ wrong in coverage test\n";
    EXPECT_EQ(SIZE_MAX, test.first_live_idx_.load()) << "SRR490124-4pairs blocks cleanup wrong in coverage test\n";
}

namespace reseq {

// Helper: bootstrap the very first block into a Prepare()'d CoverageStats.
// Mirrors the "Initialize first block" branch of EnsureSpace. Must be a
// CoverageStatsTest member to exercise friend-class private access.
CoverageStats::CoverageBlock* CoverageStatsTest::BootstrapFirstBlock(CoverageStats& cs, uintRefSeqId seq_id,
                                                                     uintSeqLen start_pos) {
    cs.blocks_.emplace_back(std::make_unique<CoverageStats::CoverageBlock>(seq_id, start_pos));
    size_t new_idx = cs.blocks_.size() - 1;
    CoverageStats::CoverageBlock* blk = cs.blocks_[new_idx].get();
    blk->coverage_.resize(CoverageStats::kBlockSize);
    blk->block_idx_ = new_idx;
    blk->prev_block_idx_ = SIZE_MAX;
    blk->next_block_idx_ = SIZE_MAX;
    cs.first_live_idx_.store(new_idx, std::memory_order_release);
    cs.last_live_idx_.store(new_idx, std::memory_order_release);
    return blk;
}

TEST_F(CoverageStatsTest, NonSystematicErrorRate) {
    string test_dir;
    ASSERT_TRUE(GetTestDir(test_dir));
    LoadReference(test_dir + "ecoli-GCF_000005845.2_ASM584v2_genomic.fa");
    CreateTestObject();

    TestNonSystematicErrorRate();
}

} // namespace reseq

// Implementation of CoverageStatsTest deque unit-test helpers.
// These are CoverageStatsTest methods so they have friend access to CoverageStats
// private members (blocks_, free_indices_, first_live_idx_, last_live_idx_, etc.).

void CoverageStatsTest::TestDequeLifecycle() {
    auto& cs = *test_;
    // Prepare with minimal plausible parameters so internal state is valid.
    cs.Prepare(/*average_coverage=*/10, /*average_read_length=*/100, /*maximum_read_length_on_reference=*/150);

    // Bootstrap the first block (index 0) directly, as EnsureSpace does.
    CoverageStats::CoverageBlock* blk0 = BootstrapFirstBlock(cs, 0, 0);
    ASSERT_EQ(1u, cs.blocks_.size());
    EXPECT_EQ(0u, blk0->block_idx_);
    EXPECT_EQ(SIZE_MAX, blk0->prev_block_idx_);
    EXPECT_EQ(SIZE_MAX, blk0->next_block_idx_);
    EXPECT_EQ(0u, cs.first_live_idx_.load());
    EXPECT_EQ(0u, cs.last_live_idx_.load());

    // Create a second block via CreateBlock (index 1).
    CoverageStats::CoverageBlock* blk1 = cs.CreateBlock(0, CoverageStats::kBlockSize);
    ASSERT_EQ(2u, cs.blocks_.size());
    EXPECT_EQ(1u, blk1->block_idx_);
    EXPECT_EQ(0u, blk1->prev_block_idx_);
    EXPECT_EQ(SIZE_MAX, blk1->next_block_idx_);
    // blk0's next_block_idx_ should have been updated.
    EXPECT_EQ(1u, blk0->next_block_idx_);
    EXPECT_EQ(1u, cs.last_live_idx_.load());
    EXPECT_TRUE(cs.free_indices_.empty());

    // RemoveBlock blk0: its index should be pushed onto free_indices_.
    cs.RemoveBlock(blk0);
    ASSERT_EQ(1u, cs.free_indices_.size());
    EXPECT_EQ(0u, cs.free_indices_.back());

    // Create a third logical block: should reuse index 0 from free_indices_.
    CoverageStats::CoverageBlock* blk2 = cs.CreateBlock(0, 2 * CoverageStats::kBlockSize);
    EXPECT_TRUE(cs.free_indices_.empty()) << "free_indices_ should be empty after reuse";
    EXPECT_EQ(0u, blk2->block_idx_) << "Reused index should be 0";
    EXPECT_EQ(2u, cs.blocks_.size()) << "Deque size must not grow when reusing";
    EXPECT_EQ(2u * CoverageStats::kBlockSize, blk2->start_pos_);
}

void CoverageStatsTest::TestFindByIndex() {
    auto& cs = *test_;
    cs.Prepare(10, 100, 150);

    // Bootstrap first block at start_pos=0, seq_id=0.
    BootstrapFirstBlock(cs, 0, 0);

    // Add a second block at start_pos=kBlockSize.
    CoverageStats::CoverageBlock* blk1 = cs.CreateBlock(0, CoverageStats::kBlockSize);

    // FindBlock with ref_pos inside block 0 should return the block with start_pos=0.
    CoverageStats::CoverageBlock* found0 = cs.FindBlock(0, 0);
    EXPECT_EQ(0u, found0->start_pos_);

    // FindBlock with ref_pos at start of block 1 should return block 1.
    CoverageStats::CoverageBlock* found1 = cs.FindBlock(0, CoverageStats::kBlockSize);
    EXPECT_EQ(blk1->start_pos_, found1->start_pos_);
    EXPECT_EQ(1u, found1->block_idx_);

    // FindBlock with a position near the end of block 1 should still return block 1.
    CoverageStats::CoverageBlock* found1b = cs.FindBlock(0, CoverageStats::kBlockSize + 5);
    EXPECT_EQ(1u, found1b->block_idx_);
}

void CoverageStatsTest::TestCleanupRecyclesIndices() {
    auto& cs = *test_;
    cs.Prepare(10, 100, 150);

    // Build a chain of three blocks: 0 -> 1 -> 2.
    BootstrapFirstBlock(cs, 0, 0);
    CoverageStats::CoverageBlock* blk1 = cs.CreateBlock(0, CoverageStats::kBlockSize);
    CoverageStats::CoverageBlock* blk2 = cs.CreateBlock(0, 2 * CoverageStats::kBlockSize);
    ASSERT_EQ(3u, cs.blocks_.size());
    EXPECT_TRUE(cs.free_indices_.empty());

    // Simulate cleanup: remove blocks 0 and 1 via RemoveBlock (as CleanUp does).
    CoverageStats::CoverageBlock* blk0 = cs.blocks_[0].get();
    cs.RemoveBlock(blk0);
    cs.RemoveBlock(blk1);
    EXPECT_EQ(2u, cs.free_indices_.size());

    // Advance first_live_idx_ to blk2 to reflect cleanup semantics.
    cs.first_live_idx_.store(blk2->block_idx_, std::memory_order_relaxed);

    // Now create two new blocks: they should reuse the freed indices.
    CoverageStats::CoverageBlock* reused_a = cs.CreateBlock(0, 3 * CoverageStats::kBlockSize);
    CoverageStats::CoverageBlock* reused_b = cs.CreateBlock(0, 4 * CoverageStats::kBlockSize);

    EXPECT_TRUE(cs.free_indices_.empty()) << "All freed indices should have been reused";
    EXPECT_EQ(3u, cs.blocks_.size()) << "Deque must not grow beyond 3 when reusing 2 indices";

    // Both reused blocks must have block_idx_ that were previously freed (0 or 1).
    EXPECT_TRUE(reused_a->block_idx_ == 0u || reused_a->block_idx_ == 1u);
    EXPECT_TRUE(reused_b->block_idx_ == 0u || reused_b->block_idx_ == 1u);
    EXPECT_NE(reused_a->block_idx_, reused_b->block_idx_) << "Each reused index must be distinct";
}

void CoverageStatsTest::TestConcurrentCoverageIncrement() {
    auto& cs = *test_;
    cs.Prepare(10, 100, 150);

    CoverageStats::CoverageBlock* blk = BootstrapFirstBlock(cs, 0, 0);

    const uintSeqLen kPos = 5;
    const int kThreads = 8;
    const int kIncrementsPerThread = 1000;
    const uintCovCount expected = static_cast<uintCovCount>(kThreads) * kIncrementsPerThread;

    // Each thread increments coverage_forward_[A] and coverage_reverse_[C] at kPos.
    std::vector<std::thread> threads;
    threads.reserve(kThreads);
    for (int t = 0; t < kThreads; ++t) {
        threads.emplace_back([&cs, blk, kPos, kIncrementsPerThread]() {
            for (int i = 0; i < kIncrementsPerThread; ++i) {
                cs.AddForward(kPos, blk, seqan::Dna5('A')); // ordinal 0 = A
                cs.AddReverse(kPos, blk, seqan::Dna5('C')); // ordinal 1 = C
            }
        });
    }
    for (auto& t : threads) {
        t.join();
    }

    EXPECT_EQ(expected, blk->coverage_.at(kPos).coverage_forward_.at(0).load())
        << "coverage_forward_[A] should equal kThreads * kIncrementsPerThread";
    EXPECT_EQ(expected, blk->coverage_.at(kPos).coverage_reverse_.at(1).load())
        << "coverage_reverse_[C] should equal kThreads * kIncrementsPerThread";
}

namespace reseq {
TEST_F(CoverageStatsTest, DequeLifecycle) {
    CreateTestObject();
    TestDequeLifecycle();
}

TEST_F(CoverageStatsTest, FindByIndex) {
    CreateTestObject();
    TestFindByIndex();
}

TEST_F(CoverageStatsTest, CleanupRecyclesIndices) {
    CreateTestObject();
    TestCleanupRecyclesIndices();
}

TEST_F(CoverageStatsTest, ConcurrentCoverageIncrement) {
    CreateTestObject();
    TestConcurrentCoverageIncrement();
}
} // namespace reseq
