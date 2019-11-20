#include "QualityStatsTest.h"
using reseq::QualityStatsTest;

#include <string>
using std::to_string;

void QualityStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new QualityStats() ) << "Could not allocate memory for QualityStats object\n";
}

void QualityStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void QualityStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void QualityStatsTest::TestSrr490124Equality(const QualityStats &test, const char *context){
	// echo "count seg ref dom tile pos qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, pos-1, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, pos-1, ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(0).at(1)[0][70].size()) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[0] not correctly shrunken for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(0).at(1)[0][70][2]) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(2).at(1)[0][86][2]) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(100, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0].size()) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[1] not correctly shrunken for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0][10].size()) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[1] not correctly shrunken for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0][10][39]) << "SRR490124-4pairs base_quality_stats_per_tile_per_error_reference_[1] not correct for " << context << '\n';

	// echo "count seg ref dom tile pos rate"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, pos-1, rate}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, pos-1, 0}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(0).at(1).at(2)[0][87][100]) << "SRR490124-4pairs error_rate_for_position_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(0).at(2).at(4)[0][97][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(1).at(0).at(4)[0][1][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_per_error_reference_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(1).at(3).at(4)[0][99][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_per_error_reference_[1] not correct for " << context << '\n';
	// echo "count seg ref tile pos rate"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; print int($1%256/128), ref, 0, pos-1, rate}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, pos-1, 0}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(1)[0][0][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(2)[0][0][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.error_rate_for_position_per_tile_reference_.at(0).at(0)[0][99][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(3)[0][93][100]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(1).at(0)[0][99][100]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(1).at(3)[0][99][0]) << "SRR490124-4pairs error_rate_for_position_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref dom tile rate qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, rate, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, 0, ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(0).at(1).at(2)[0][100][2]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(22, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(0).at(2).at(4)[0][0][2]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_per_error_reference_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(1).at(0).at(4)[0][0][32]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_per_error_reference_[1] not correct for " << context << '\n';
	EXPECT_EQ(35, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(1).at(3).at(4)[0][0][2]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_per_error_reference_[1] not correct for " << context << '\n';
	// echo "count seg ref tile rate qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; print int($1%256/128), ref, 0, rate, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, 0, ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(101, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0].size()) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0][100].size()) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(3, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0][100][2]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(101, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0].size()) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(38, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][0].size()) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(35, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][0][2]) << "SRR490124-4pairs base_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';

	// If there is no previous quality(first base of read) previous quality is set to 1
	// echo "count seg ref tile"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk '{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; print int($1%256/128), ref, 0}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0}}' | sort -n | uniq -c
	EXPECT_EQ(51, SumVect(test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0])) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(54, SumVect(test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(1)[0])) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[1] wrong for " << context << '\n';
	// echo "count seg ref tile prev qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(1==pos){prev=1}else{prev=ord[substr($4,pos-1,1)]-33};print int($1%256/128), ref, 0, prev, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, ord[substr($4,pos-1,1)]-33, ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(4, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37].size()) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][35]) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][36]) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][37]) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][38]) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(0, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(1)[0][33].size()) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(2)[0][33].size()) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(2)[0][33][36]) << "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile rate prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; if(1==pos){prev=1}else{prev=ord[substr($4,pos-1,1)]-33};print int($1%256/128), ref, 0, rate, prev}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, 0, ord[substr($4,pos-1,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(38, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0].size()) << "SRR490124-4pairs preceding_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(15, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0][2]) << "SRR490124-4pairs preceding_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(7, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(2)[0][0][38]) << "SRR490124-4pairs preceding_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_error_rate_per_tile_reference_.at(1).at(1)[0][100].size()) << "SRR490124-4pairs preceding_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(8, test.preceding_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][100][2]) << "SRR490124-4pairs preceding_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile pos prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(1==pos){prev=1}else{prev=ord[substr($4,pos-1,1)]-33};print int($1%256/128), ref, 0, pos-1, prev}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, pos-1, ord[substr($4,pos-1,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(0).at(0)[0][6][38]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(93, test.preceding_quality_for_position_per_tile_reference_.at(0).at(1)[0].size()) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(0).at(2)[0][4][38]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(1).at(2)[0][0][1]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(100, test.preceding_quality_for_position_per_tile_reference_.at(1).at(3)[0].size()) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(1).at(3)[0][99][2]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile sq qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($4);pos+=1){sq+=ord[substr($4,pos,1)]-33}; sq=int(sq/length($4)+0.5); pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; print int($1%256/128), ref, 0, sq, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, sq, ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(2, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0].size()) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(37, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0][19].size()) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(6, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0][19][2]) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0].size()) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(38, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13].size()) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(27, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13][2]) << "SRR490124-4pairs base_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile sq prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($4);pos+=1){sq+=ord[substr($4,pos,1)]-33}; sq=int(sq/length($4)+0.5); pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(1==pos){prev=1}else{prev=ord[substr($4,pos-1,1)]-33};print int($1%256/128), ref, 0, sq, prev}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, sq, ord[substr($4,pos-1,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(38, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(1)[0][19].size()) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(15, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(2)[0][19][2]) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(7, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(1)[0][19][38]) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(38, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(1)[0][13].size()) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(27, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13][2]) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(2, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(1)[0][13][39]) << "SRR490124-4pairs preceding_quality_for_sequence_quality_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile rate sq"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($4);pos+=1){sq+=ord[substr($4,pos,1)]-33}; sq=int(sq/length($4)+0.5); pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b;rate=100}else{ref=substr($3,pos,1);rate=0}; print int($1%256/128), ref, 0, rate, sq}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, 0, sq}}' | sort -n | uniq -c
	EXPECT_EQ(2, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0].size()) << "SRR490124-4pairs sequence_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(21, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0][19]) << "SRR490124-4pairs sequence_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(22, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(2)[0][0][20]) << "SRR490124-4pairs sequence_quality_for_error_rate_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_for_error_rate_per_tile_reference_.at(1).at(1)[0][100].size()) << "SRR490124-4pairs sequence_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(5, test.sequence_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][100][13]) << "SRR490124-4pairs sequence_quality_for_error_rate_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg ref tile pos sq"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($4);pos+=1){sq+=ord[substr($4,pos,1)]-33}; sq=int(sq/length($4)+0.5); pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; print int($1%256/128), ref, 0, pos-1, sq}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 0, pos-1, sq}}' | sort -n | uniq -c
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(1)[0][0][19]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(2)[0][0][20]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(0)[0][99][19]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(0)[0][99][20]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(1).at(0)[0][99][12]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(1).at(3)[0][99][13]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_reference_[1] wrong for " << context << '\n';

	// echo "count seg tile gc sq"; samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($11);pos+=1){sq+=ord[substr($11,pos,1)]-33}; sq=int(sq/length($11)+0.5); print int($2%256/128), sq; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1 " | seqtk seq")}' | awk '(1==NR%3){seg=$1;sq=$2}(0==NR%3){print seg, 0, gsub("G","",$0)+gsub("C","",$0), sq}' | sort | uniq -c
	EXPECT_EQ(4, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0][54][20] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[0][0][54] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0][57][19] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[0][0][57] wrong for " << context << '\n';
	EXPECT_EQ(7, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0][46][13] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[1][0][37] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0][52][12] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_reference_[1][0][46] wrong for " << context << '\n';

	// echo "count seg tile rate sq"; samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0;for(pos=1;pos<=length($11);pos+=1){sq+=ord[substr($11,pos,1)]-33}; sq=int(sq/length($11)+0.5); print int($2%256/128), 0, gsub(/[ACGT]+/,"",$18)/length($10)*100, sq}' | sort | uniq -c
	EXPECT_EQ(3, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0][10][20] ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[0][0][52] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0][12][19] ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[0][0][59] wrong for " << context << '\n';
	EXPECT_EQ(5, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0][11][13] ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[1][0][37] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0][15][12] ) << "SRR490124-4pairs sequence_quality_mean_for_mean_error_rate_per_tile_reference_[1][0][47] wrong for " << context << '\n';

	// echo "count seg tile gc rate"; samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '{print int($2%256/128), gsub(/[ACGT]+/,"",$18)/length($10)*100; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1 " | seqtk seq")}' | awk '(1==NR%3){seg=$1;rate=$2}(0==NR%3){print seg, 0, gsub("G","",$0)+gsub("C","",$0), rate}' | sort | uniq -c
	EXPECT_EQ(4, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0][54][10] ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[0][0][52] wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0][57][12] ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[0][0][59] wrong for " << context << '\n';
	EXPECT_EQ(7, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0][46][11] ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[1][0][37] wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0][52][15] ) << "SRR490124-4pairs mean_error_rate_for_gc_per_tile_reference_[1][0][47] wrong for " << context << '\n';

	// Manually from the 3 above + bowtie fragment length
	EXPECT_EQ(59, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0][126][20] ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[0][0][126] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0][184][19] ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[0][0][184] wrong for " << context << '\n';
	EXPECT_EQ(59, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0][184][13] ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[1][0][184] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0][126][12] ) << "SRR490124-4pairs sequence_quality_mean_for_fragment_length_per_tile_reference_[1][0][126] wrong for " << context << '\n';

	EXPECT_EQ(59, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0][126][10] ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[0][0][126] wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0][184][12] ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[0][0][184] wrong for " << context << '\n';
	EXPECT_EQ(59, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0][184][11] ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[1][0][184] wrong for " << context << '\n';
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0][126][15] ) << "SRR490124-4pairs mean_error_rate_for_fragment_length_per_tile_reference_[1][0][126] wrong for " << context << '\n';

	EXPECT_EQ(59, test.gc_for_fragment_length_per_tile_reference_.at(0)[0].size() ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(0)[0][126][54] ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[0][0][126] wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(0)[0][184][57] ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[0][0][184] wrong for " << context << '\n';
	EXPECT_EQ(59, test.gc_for_fragment_length_per_tile_reference_.at(1)[0].size() ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(1)[0][184][46] ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[1][0][184] wrong for " << context << '\n';
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(1)[0][126][52] ) << "SRR490124-4pairs gc_for_fragment_length_per_tile_reference_[1][0][126] wrong for " << context << '\n';

	// Manually
	for( int templ_seg=2; templ_seg--; ){
		auto tmp_comparison( test.base_quality_for_sequence_per_tile_.at(templ_seg).at(0)[0] );
		for(auto i=4; --i; ){
			tmp_comparison += test.base_quality_for_sequence_per_tile_.at(templ_seg).at(i)[0];
		}
		ShrinkVect(tmp_comparison);
		TestVectEquality2d(test.base_quality_for_sequence_.at(templ_seg), tmp_comparison, context, "SRR490124-4pairs base_quality_for_sequence_per_tile_[" + to_string(templ_seg) + ']', " not identical with base_quality_for_sequence_ for ");
	}

	for( int templ_seg=2; templ_seg--; ){
		auto tmp_comparison2( test.base_quality_for_preceding_quality_per_tile_.at(templ_seg).at(0)[0] );
		for(auto i=4; --i; ){
			tmp_comparison2 += test.base_quality_for_preceding_quality_per_tile_.at(templ_seg).at(i)[0];
		}
		ShrinkVect(tmp_comparison2);
		TestVectEquality2d(test.base_quality_for_preceding_quality_.at(templ_seg), tmp_comparison2, context, "SRR490124-4pairs base_quality_for_preceding_quality_per_tile_[" + to_string(templ_seg) + ']', " not identical with base_quality_for_preceding_quality_ for ");
	}

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(p=1;p<=length(poslist);p+=1){pos=poslist[p];print int($1%256/128), substr($2,pos,1), 0, pos-1, ord[substr($3,pos,1)]-33}}' | sort | uniq -c | sort -k2,4 -k5,6n
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(0).at(0)[0], 3, 1, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[0][0][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(0).at(1)[0], 1, 2, 0, 0, context, "SRR490124-4pairs base_quality_stats_per_tile_[0][1][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(0).at(2)[0], 1, 0, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[0][2][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(0).at(3)[0], 0, 1, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[0][3][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(1).at(0)[0], 1, 10, 1, 0, context, "SRR490124-4pairs base_quality_stats_per_tile_[1][0][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(1).at(1)[0], 1, 0, 0, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[1][1][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(1).at(2)[0], 1, 0, 0, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[1][2][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_.at(1).at(3)[0], 1, 4, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_[1][3][0] not correctly shrunken for ");

	// echo "count seg base tile sq prev"; cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0; sq=0;for(pos=1;pos<=length($3);pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/length($3)+0.5); for(pos=1;pos<=length($2);pos+=1){if(1==pos){prev=1}else{prev=ord[substr($3,pos-1,1)]-33};print int($1%256/128), substr($2,pos,1), 0, sq, prev}}' | sort | uniq -c | sort -k2,4 -k5,6n
	EXPECT_EQ(38, test.preceding_quality_for_sequence_per_tile_.at(0).at(1)[0][19].size()) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(7, test.preceding_quality_for_sequence_per_tile_.at(0).at(1)[0][19][38]) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(16, test.preceding_quality_for_sequence_per_tile_.at(0).at(2)[0][19][2]) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(38, test.preceding_quality_for_sequence_per_tile_.at(1).at(1)[0][21].size()) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_sequence_per_tile_.at(1).at(1)[0][21][39]) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[1] wrong for " << context << '\n';
	EXPECT_EQ(5, test.preceding_quality_for_sequence_per_tile_.at(1).at(3)[0][21][2]) << "SRR490124-4pairs preceding_quality_for_sequence_per_tile_[1] wrong for " << context << '\n';

	// echo "count seg base tile pos prev"; cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(pos=1;pos<=length($2);pos+=1){if(1==pos){prev=1}else{prev=ord[substr($3,pos-1,1)]-33};print int($1%256/128), substr($2,pos,1), 0, pos-1, prev}}' | sort | uniq -c | sort -k2,4 -k5,6n
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_.at(0).at(0)[0][4][39]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(94, test.preceding_quality_for_position_per_tile_.at(0).at(1)[0].size()) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_.at(0).at(2)[0][4][38]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_.at(1).at(2)[0][0][1]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[1] wrong for " << context << '\n';
	EXPECT_EQ(100, test.preceding_quality_for_position_per_tile_.at(1).at(3)[0].size()) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[1] wrong for " << context << '\n';
	EXPECT_EQ(2, test.preceding_quality_for_position_per_tile_.at(1).at(3)[0][99][2]) << "SRR490124-4pairs preceding_quality_for_position_per_tile_[1] wrong for " << context << '\n';

	// echo "count seg base tile pos sq"; cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0; sq=0;for(pos=1;pos<=length($3);pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/length($3)+0.5); for(pos=1;pos<=length($2);pos+=1){print int($1%256/128), substr($2,pos,1), 0, pos-1, sq}}' | sort | uniq -c | sort -k2,4 -k5,6n
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_.at(0).at(0)[0][0][20]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_.at(0).at(1)[0][0][19]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_for_position_per_tile_.at(0).at(3)[0][29][20]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_.at(0).at(0)[0][99][19]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_.at(1).at(3)[0][99][13]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_.at(1).at(3)[0][99][21]) << "SRR490124-4pairs sequence_quality_for_position_per_tile_[1] wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(p=1;p<=length(poslist);p+=1){pos=poslist[p];print pos-1, ord[substr($3,pos,1)]-33}}' | sort | uniq -c | sort -k2,4 -k5,6n
	// All reads are on the reverse strand
	EXPECT_EQ(0, test.base_quality_stats_per_strand_.at(0).size()) << "SRR490124-4pairs base_quality_stats_per_strand_[0] wrong for " << context << '\n';
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_strand_.at(1), 7, 10, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_strand_[1] not correctly shrunken for ");

	// echo "count seg nuc tile nuc_perc sq"; cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("A,C,G,T", baselist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{sq=0; for(pos=1;pos<=length($3);pos+=1){sq+=ord[substr($3,pos,1)]-33}; len=length($3); sq=int(sq/len+0.5); for(plist=1;plist<=length(baselist);plist+=1){base=baselist[plist]; print int($1%256/128), base, 0, int(gsub(base,"",$2)*100/len+0.5), sq}}' | sort | uniq -c
	EXPECT_EQ(9, test.sequence_quality_for_base_per_tile_.at(0).at(0)[0].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(0)[0][24].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][0][0][24].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(0)[0][24][20] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][0][0][24] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(0)[0][25].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][0][0][25].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(0)[0][25][19] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][0][0][25] wrong for " << context << '\n';
	EXPECT_EQ(8, test.sequence_quality_for_base_per_tile_.at(0).at(1)[0].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(1)[0][32].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][32].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(1)[0][32][20] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][32] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(1)[0][25].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][25].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(1)[0][25][19] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][25] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(2)[0][20][20] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][20] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(2)[0][34][19] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][34] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(3)[0][24][20] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][24] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(0).at(3)[0][16][19] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[0][1][0][16] wrong for " << context << '\n';
	EXPECT_EQ(13, test.sequence_quality_for_base_per_tile_.at(1).at(3)[0].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[1][3][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(1).at(3)[0][25].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[1][3][0][25].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(1).at(3)[0][25][21] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[1][3][0][25] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(1).at(3)[0][32].size() ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[1][3][0][32].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_for_base_per_tile_.at(1).at(3)[0][32][13] ) << "SRR490124-4pairs sequence_quality_for_base_per_tile_[1][3][0][32] wrong for " << context << '\n';

	// echo "count seg tile read_gc sq"; samtools view ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($10); sq=0; for(pos=1;pos<=len;pos+=1){sq+=ord[substr($11,pos,1)]-33}; sq=int(sq/len+0.5); print int($2%256/128), 0, int(gsub(/[GC]/,"",$10)*100/len), sq}' | sort | uniq -c
	EXPECT_EQ(12, test.sequence_quality_mean_for_gc_per_tile_.at(0)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[0][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_.at(0)[0][52][20] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[0][0][52] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_.at(0)[0][59][19] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[0][0][59] wrong for " << context << '\n';
	EXPECT_EQ(20, test.sequence_quality_mean_for_gc_per_tile_.at(1)[0].size() ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[1][0].size() wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_.at(1)[0][37][21] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[1][0][37] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_.at(1)[0][47][13] ) << "SRR490124-4pairs sequence_quality_mean_for_gc_per_tile_[1][0][47] wrong for " << context << '\n';

	// echo "count seg sq"; samtools view ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($10); sq=0; for(pos=1;pos<=len;pos+=1){sq+=10**(-(ord[substr($11,pos,1)]-33)/10)}; sq=int(-10*log(sq/length($11))/log(10)+0.5); print int($2%256/128), sq}' | sort | uniq -c
	EXPECT_EQ(3, test.sequence_quality_probability_mean_.at(0).size() ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_probability_mean_.at(0)[6] ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_probability_mean_.at(0)[5] ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';
	EXPECT_EQ(5, test.sequence_quality_probability_mean_.at(1).size() ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_probability_mean_.at(1)[6] ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_probability_mean_.at(1)[4] ) << "SRR490124-4pairs sequence_quality_probability_mean_ wrong for " << context << '\n';

	// samtools view ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($10); seg=int($2%256/128); for(pos=1;pos<=len;pos+=1){print NR, seg, ord[substr($11,pos,1)]-33, len}}' | sort -k2,2 -k1,1n -k3,3n | uniq -c | awk 'BEGIN{n=0}{if($2==n){sum+=$1;max=$4}else{if(0!=min){print "1Minimum:", seg, min; print "2First Quartile:", seg, first; print "3Median:", seg, median; print "4Third Quartile:", seg, third; print "5Max:", seg, max;}; n=$2; seg=$3; sum=$1; min=$4; first=0; median=0; third=0; max=$4};if(0==first&&$5/4<sum){first=$4};if(0==median&&$5/2<sum){median=$4};if(0==third&&$5*3/4<sum){third=$4}}END{print "1Minimum:", seg, min; print "2First Quartile:", seg, first; print "3Median:", seg, median; print "4Third Quartile:", seg, third; print "5Max:", seg, max;}' | sort | uniq -c
	for( int templ_seg=2; templ_seg--; ){
		EXPECT_EQ(1, test.sequence_quality_minimum_.at(templ_seg).size()) << "SRR490124-4pairs sequence_quality_minimum_[" << templ_seg << "] wrong for " << context << '\n';
		EXPECT_EQ(4, test.sequence_quality_minimum_.at(templ_seg)[2]) << "SRR490124-4pairs sequence_quality_minimum_[" << templ_seg << "] wrong for " << context << '\n';
		EXPECT_EQ(1, test.sequence_quality_first_quartile_.at(templ_seg).size()) << "SRR490124-4pairs sequence_quality_first_quartile_[" << templ_seg << "] wrong for " << context << '\n';
		EXPECT_EQ(4, test.sequence_quality_first_quartile_.at(templ_seg)[2]) << "SRR490124-4pairs sequence_quality_first_quartile_[" << templ_seg << "] wrong for " << context << '\n';
	}
	EXPECT_EQ(26, test.sequence_quality_median_.at(0).size()) << "SRR490124-4pairs sequence_quality_median_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_median_.at(0)[26]) << "SRR490124-4pairs sequence_quality_median_[0] wrong for " << context << '\n';
	EXPECT_EQ(26, test.sequence_quality_median_.at(1).size()) << "SRR490124-4pairs sequence_quality_median_[1] wrong for " << context << '\n';
	EXPECT_EQ(3, test.sequence_quality_median_.at(1)[2]) << "SRR490124-4pairs sequence_quality_median_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_median_.at(1)[27]) << "SRR490124-4pairs sequence_quality_median_[1] wrong for " << context << '\n';
	EXPECT_EQ(6, test.sequence_quality_third_quartile_.at(0).size()) << "SRR490124-4pairs sequence_quality_third_quartile_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_third_quartile_.at(0)[36]) << "SRR490124-4pairs sequence_quality_third_quartile_[0] wrong for " << context << '\n';
	EXPECT_EQ(35, test.sequence_quality_third_quartile_.at(1).size()) << "SRR490124-4pairs sequence_quality_third_quartile_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_third_quartile_.at(1)[29]) << "SRR490124-4pairs sequence_quality_third_quartile_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_third_quartile_.at(1)[36]) << "SRR490124-4pairs sequence_quality_third_quartile_[1] wrong for " << context << '\n';
	EXPECT_EQ(3, test.sequence_quality_maximum_.at(0).size()) << "SRR490124-4pairs sequence_quality_maximum_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_maximum_.at(0)[38]) << "SRR490124-4pairs sequence_quality_maximum_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_maximum_.at(0)[39]) << "SRR490124-4pairs sequence_quality_maximum_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_maximum_.at(1).size()) << "SRR490124-4pairs sequence_quality_maximum_[1] wrong for " << context << '\n';
	EXPECT_EQ(3, test.sequence_quality_maximum_.at(1)[39]) << "SRR490124-4pairs sequence_quality_maximum_[1] wrong for " << context << '\n';

	// cat "count seg qual occurances"; samtools view ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($10); seg=int($2%256/128); for(pos=1;pos<=len;pos+=1){print NR, seg, ord[substr($11,pos,1)]-33}}' | sort -k2,2 -k1,1n -k3,3n | uniq -c | awk '{print $3, $4, $1}' | sort -k1,1n -k2,2n -k3,3n | uniq -c
	for( int templ_seg=2; templ_seg--; ){
		EXPECT_EQ(38, test.sequence_quality_content_.at(templ_seg).size()) << "SRR490124-4pairs sequence_quality_content_[" << templ_seg << "].size() wrong for " << context << '\n';
	}
	EXPECT_EQ(21, test.sequence_quality_content_.at(0)[2].size()) << "SRR490124-4pairs sequence_quality_content_[0][2].size() wrong for " << context << '\n';
	EXPECT_EQ(0, test.sequence_quality_content_.at(0)[4].size()) << "SRR490124-4pairs sequence_quality_content_[0][4].size() wrong for " << context << '\n';
	TestVectEquality({0,{2,1,1}}, test.sequence_quality_content_.at(0)[20], context, "SRR490124-4pairs sequence_quality_content_[0][20]", " wrong for ");
	EXPECT_EQ(14, test.sequence_quality_content_.at(0)[39].size()) << "SRR490124-4pairs sequence_quality_content_[0][39].size() wrong for " << context << '\n';
	EXPECT_EQ(60, test.sequence_quality_content_.at(1)[2].size()) << "SRR490124-4pairs sequence_quality_content_[1][2].size() wrong for " << context << '\n';
	TestVectEquality({0,{3,1}}, test.sequence_quality_content_.at(1)[4], context, "SRR490124-4pairs sequence_quality_content_[1][4]", " wrong for ");
	TestVectEquality({0,{1,1,2}}, test.sequence_quality_content_.at(1)[20], context, "SRR490124-4pairs sequence_quality_content_[1][20]", " wrong for ");
	TestVectEquality({0,{1,0,0,0,0,1,0,1,1}}, test.sequence_quality_content_.at(1)[39], context, "SRR490124-4pairs sequence_quality_content_[1][39]", " wrong for ");

	// samtools view ecoli-SRR490124-4pairs.sam | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{last_qual=0; for(pos=1;pos<=length($11);pos+=1){qual=ord[substr($11,pos,1)]-33;if(qual==last_qual){len+=1}else{if(0!=last_qual){print last_qual, len}; last_qual=qual;len=1}}; print last_qual, len}' | sort -k1,1n -k2,2n | uniq -c
	EXPECT_EQ(60, test.homoquality_distribution_[2].size()) << "SRR490124-4pairs homoquality_distribution_[2].size() wrong for " << context << '\n';
	TestVectEquality({1,{1}}, test.homoquality_distribution_[4], context, "SRR490124-4pairs homoquality_distribution_[4]", " not correct for ");
	TestVectEquality({1,{8}}, test.homoquality_distribution_[20], context, "SRR490124-4pairs homoquality_distribution_[20]", " not correct for ");
	TestVectEquality({1,{12,5,2,0,1}}, test.homoquality_distribution_[39], context, "SRR490124-4pairs homoquality_distribution_[39]", " not correct for ");

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{seg=int($1%256/128); for(pos=1;pos<=length($3);pos+=1){print seg, substr($2,pos,1), ord[substr($3,pos,1)]-33}}' | sort -k1,1 -k2,2 -k3,3n | uniq -c | awk 'BEGIN{base=0}{if($3==base){sum+=$1;mean+=$1*$4;max=$4;quals[$4]=$1; while(medsum<=sum/2){medsum+=quals[++median]}}else{if(0!=min){print seg, base, "Minimum:", min; print seg, base, "Mean:", int(mean/sum+0.5); print seg, base, "Median:", median; print seg, base, "Max:", max}; base=$3; seg=$2; sum=$1; mean=$1*$4; min=$4; medsum=$1; median=$4; max=$4;delete quals; quals[$4]=$1}}END{print seg, base, "Minimum:", min; print seg, base, "Mean:", int(mean/sum+0.5); print seg, base, "Median:", median; print seg, base, "Max:", max}'
	TestSeqQualityStats(test.nucleotide_quality_.at(0).at(0), 38, 18, 2, 16, 39, context, "of SRR490124-4pairs nucleotide_quality_[0][0] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(0).at(1), 38, 19, 2, 22, 39, context, "of SRR490124-4pairs nucleotide_quality_[0][1] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(0).at(2), 38, 16, 2, 2, 39, context, "of SRR490124-4pairs nucleotide_quality_[0][2] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(0).at(3), 38, 19, 2, 26, 39, context, "of SRR490124-4pairs nucleotide_quality_[0][3] wrong for ");
	EXPECT_EQ(0, test.nucleotide_quality_.at(0).at(4).size()) << "size of SRR490124-4pairs nucleotide_quality_[0][4] wrong for " << context << '\n';
	TestSeqQualityStats(test.nucleotide_quality_.at(1).at(0), 38, 16, 2, 2, 39, context, "of SRR490124-4pairs nucleotide_quality_[1][0] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(1).at(1), 38, 10, 2, 2, 39, context, "of SRR490124-4pairs nucleotide_quality_[1][1] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(1).at(2), 38, 12, 2, 2, 39, context, "of SRR490124-4pairs nucleotide_quality_[1][2] wrong for ");
	TestSeqQualityStats(test.nucleotide_quality_.at(1).at(3), 38, 12, 2, 2, 39, context, "of SRR490124-4pairs nucleotide_quality_[1][3] wrong for ");
	EXPECT_EQ(0, test.nucleotide_quality_.at(1).at(4).size()) << "size of SRR490124-4pairs nucleotide_quality_[1][4] wrong for " << context << '\n';

	// cat <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(p=1;p<=length(poslist);p+=1){pos=poslist[p];print int($1%256/128), pos-1, ord[substr($3,pos,1)]-33}}' | sort | uniq -c | sort -k2,2 -k3,4n
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_reference_.at(0), 1, 2, 1, 1, context, "SRR490124-4pairs base_quality_stats_reference_[0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_reference_.at(1), 1, 10, 1, 1, context, "SRR490124-4pairs base_quality_stats_reference_[1] not correctly shrunken for ");
	// cat <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print $2; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq -r) <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq) | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}' | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(p=1;p<=length(poslist);p+=1){pos=poslist[p];print int($1%256/128), substr($2,pos,1), 0, pos-1, ord[substr($3,pos,1)]-33}}' | sort | uniq -c | sort -k2,4 -k5,6n
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_reference_.at(0).at(0)[0], 0, 0, 1, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_reference_[0][0][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_reference_.at(1).at(0)[0], 0, 1, 0, 1, context, "SRR490124-4pairs base_quality_stats_per_tile_reference_[1][0][0] not correctly shrunken for ");

	// cat <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print $2; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq -r) <(samtools view -q10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2; system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq) | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}' | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); sq=0; for(pos=1;pos<=len;pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/len+0.5); print int($1%256/128), int(gsub(/[GC]/,"",$2)*100/len), sq}'
	EXPECT_EQ(20, test.average_sequence_quality_for_gc_.at(0)[54] ) << "SRR490124-4pairs average_sequence_quality_for_gc_[0] wrong for " << context << '\n';
	EXPECT_EQ(19, test.average_sequence_quality_for_gc_.at(0)[57] ) << "SRR490124-4pairs average_sequence_quality_for_gc_[0] wrong for " << context << '\n';
	EXPECT_EQ(13, test.average_sequence_quality_for_gc_.at(1)[46] ) << "SRR490124-4pairs average_sequence_quality_for_gc_[1] wrong for " << context << '\n';
	EXPECT_EQ(12, test.average_sequence_quality_for_gc_.at(1)[52] ) << "SRR490124-4pairs average_sequence_quality_for_gc_[1] wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); for(p=1;p<=length(poslist);p+=1){pos=poslist[p]; print pos-1, seg, ord[substr($3,pos,1)]-33, 4}}' | sort -k2,2 -k1,1n -k3,3n | uniq -c | awk 'BEGIN{n=-1}{if($2==n){sum+=$1;mean+=$1*$4;max=$4}else{if(0!=min){print n, seg, "Mean:", int(mean/size+0.5); print n, seg, "Minimum:", min; print n, seg, "First Quartile:", first; print n, seg, "Median:", median; print n, seg, "Third Quartile:", third; print n, seg, "Max:", max;}; n=$2; seg=$3; sum=$1; size=$5; mean=$1*$4; min=$4; first=0; median=0; third=0; max=$4};if(0==first&&$5/4<=sum){first=$4};if(0==median&&$5/2<sum){median=$4};if(0==third&&$5*3/4<sum){third=$4}}END{print n, seg, "Mean:", int(mean/size+0.5); print n, seg, "Minimum:", min; print n, seg, "First Quartile:", first; print n, seg, "Median:", median; print n, seg, "Third Quartile:", third; print n, seg, "Max:", max;}' | sort -k3,3 -k2,2 -k1,1n
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_.at(0), 3, 3, 1, 1, context, "SRR490124-4pairs base_quality_stats_[0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_.at(1), 7, 10, 1, 1, context, "SRR490124-4pairs base_quality_stats_[1] not correctly shrunken for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_mean_.at(0), 38, 38, 2, 2, context, "SRR490124-4pairs base_quality_mean_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_mean_.at(1), 37, 36, 2, 2, context, "SRR490124-4pairs base_quality_mean_[1] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_minimum_.at(0), 37, 37, 2, 2, context, "SRR490124-4pairs base_quality_minimum_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_minimum_.at(1), 33, 30, 2, 2, context, "SRR490124-4pairs base_quality_minimum_[1] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_first_quartile_.at(0), 37, 37, 2, 2, context, "SRR490124-4pairs base_quality_first_quartile_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_first_quartile_.at(1), 33, 30, 2, 2, context, "SRR490124-4pairs base_quality_first_quartile_[1] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_median_.at(0), 38, 38, 2, 2, context, "SRR490124-4pairs base_quality_median_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_median_.at(1), 39, 39, 2, 2, context, "SRR490124-4pairs base_quality_median_[1] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_third_quartile_.at(0), 39, 39, 2, 2, context, "SRR490124-4pairs base_quality_third_quartile_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_third_quartile_.at(1), 39, 39, 2, 2, context, "SRR490124-4pairs base_quality_third_quartile_[1] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_maximum_.at(0), 39, 39, 2, 2, context, "SRR490124-4pairs base_quality_maximum_[0] wrong for ");
	TestValueAtIndex_0_1_98_99(test.base_quality_maximum_.at(1), 39, 39, 2, 2, context, "SRR490124-4pairs base_quality_maximum_[1] wrong for ");

	// Only one tile in the test data
	EXPECT_EQ(0, test.tile_quality_mean_difference_.at(0).size()) << "SRR490124-4pairs tile_quality_mean_difference_ wrong for " << context << '\n';
	EXPECT_EQ(0, test.tile_quality_mean_difference_.at(1).size()) << "SRR490124-4pairs tile_quality_mean_difference_ wrong for " << context << '\n';

	// Only reverse strand in the test data
	EXPECT_EQ(0, test.base_quality_mean_per_strand_.at(0).size()) << "SRR490124-4pairs base_quality_mean_per_strand_[0] wrong for " << context << '\n';
	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("1,2,99,100", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); for(p=1;p<=length(poslist);p+=1){pos=poslist[p]; print pos-1, ord[substr($3,pos,1)]-33}}' | sort -k1,1n -k2,2n | uniq -c | awk 'BEGIN{n=-1}{if($2==n){sum+=$1;mean+=$1*$3}else{if(0!=mean){print n, "Mean:", int(mean/sum+0.5);}; n=$2; sum=$1; mean=$1*$3}}END{print n, seg, "Mean:", int(mean/sum+0.5)}' | sort -k1,1n
	TestValueAtIndex_0_1_98_99(test.base_quality_mean_per_strand_.at(1), 38, 37, 2, 2, context, "SRR490124-4pairs base_quality_mean_per_strand_[1] wrong for ");

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); sq=0;for(pos=1;pos<=len;pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/len+0.5);for(pos=1;pos<=len;pos+=1){print seg, sq, ord[substr($3,pos,1)]-33}}' | sort -k1,1n -k2,2n -k3,3n | uniq -c
	EXPECT_EQ(8, test.base_quality_for_sequence_.at(0).size()) << "SRR490124-4pairs base_quality_for_sequence_[0] wrong for " << context << '\n';
	EXPECT_EQ(37, test.base_quality_for_sequence_.at(0)[19].size()) << "SRR490124-4pairs base_quality_for_sequence_[0] wrong for " << context << '\n';
	EXPECT_EQ(45, test.base_quality_for_sequence_.at(0)[19][2]) << "SRR490124-4pairs base_quality_for_sequence_[0] wrong for " << context << '\n';
	EXPECT_EQ(18, test.base_quality_for_sequence_.at(1).size()) << "SRR490124-4pairs base_quality_for_sequence_[1] wrong for " << context << '\n';
	EXPECT_EQ(38, test.base_quality_for_sequence_.at(1)[21].size()) << "SRR490124-4pairs base_quality_for_sequence_[1] wrong for " << context << '\n';
	EXPECT_EQ(34, test.base_quality_for_sequence_.at(1)[21][2]) << "SRR490124-4pairs base_quality_for_sequence_[1] wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); prev=1; for(pos=1;pos<=len;pos+=1){print seg, prev, ord[substr($3,pos,1)]-33; prev=ord[substr($3,pos,1)]-33}}' | sort -k1,1n -k2,2n -k3,3n | uniq -c
	EXPECT_EQ(400, SumVect(test.base_quality_for_preceding_quality_.at(0))) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(13, test.base_quality_for_preceding_quality_.at(0)[37].size()) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_.at(0)[37][27]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_.at(0)[37][34]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(4, test.base_quality_for_preceding_quality_.at(0)[37][35]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(5, test.base_quality_for_preceding_quality_.at(0)[37][36]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(9, test.base_quality_for_preceding_quality_.at(0)[37][37]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(4, test.base_quality_for_preceding_quality_.at(0)[37][38]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(3, test.base_quality_for_preceding_quality_.at(0)[37][39]) << "SRR490124-4pairs base_quality_for_preceding_quality_[0] wrong for " << context << '\n';
	EXPECT_EQ(400, SumVect(test.base_quality_for_preceding_quality_.at(1))) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(38, test.base_quality_for_preceding_quality_.at(1)[33].size()) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_.at(1)[33][2]) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_.at(1)[33][24]) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(5, test.base_quality_for_preceding_quality_.at(1)[33][36]) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_.at(1)[33][37]) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_.at(1)[33][39]) << "SRR490124-4pairs base_quality_for_preceding_quality_[1] wrong for " << context << '\n';

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); sq=0;for(pos=1;pos<=len;pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/len+0.5);print seg, sq}' | sort -k1,1n -k2,2n | uniq -c
	EXPECT_EQ(8, test.sequence_quality_mean_.at(0).size()) << "SRR490124-4pairs sequence_quality_mean_[0] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_.at(0)[19]) << "SRR490124-4pairs sequence_quality_mean_[0] wrong for " << context << '\n';
	EXPECT_EQ(2, test.sequence_quality_mean_.at(0)[20]) << "SRR490124-4pairs sequence_quality_mean_[0] wrong for " << context << '\n';
	EXPECT_EQ(18, test.sequence_quality_mean_.at(1).size()) << "SRR490124-4pairs sequence_quality_mean_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_.at(1)[13]) << "SRR490124-4pairs sequence_quality_mean_[1] wrong for " << context << '\n';
	EXPECT_EQ(1, test.sequence_quality_mean_.at(1)[21]) << "SRR490124-4pairs sequence_quality_mean_[1] wrong for " << context << '\n';
	for( int templ_seg=2; templ_seg--; ){
		TestVectEquality(test.sequence_quality_mean_.at(templ_seg), test.sequence_quality_mean_per_tile_.at(templ_seg)[0], context, "SRR490124-4pairs sequence_quality_mean_per_tile_[" + to_string(templ_seg) + ']', " not identical with sequence_quality_mean_ for ");
	}

	// cat <(samtools view ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" NR, $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}') <(samtools view ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, $10, $11}') | awk 'BEGIN{split("A,C,G,T", baselist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); sq=0;for(pos=1;pos<=len;pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/len+0.5); for(p=1;p<=length(baselist);p+=1){base = baselist[p];print seg, base, gsub(base,"",$2), sq}}' | sort -k1,2 -k3,3n -k4,4n | uniq -c
	EXPECT_EQ(8, test.average_sequence_quality_for_base_.at(0).at(1).size() ) << "SRR490124-4pairs average_sequence_quality_for_base_[0][1].size() wrong for " << context << '\n';
	EXPECT_EQ(19, test.average_sequence_quality_for_base_.at(0).at(1)[25] ) << "SRR490124-4pairs average_sequence_quality_for_base_[0][1][25] wrong for " << context << '\n';
	EXPECT_EQ(20, test.average_sequence_quality_for_base_.at(0).at(1)[32] ) << "SRR490124-4pairs average_sequence_quality_for_base_[0][1][32] wrong for " << context << '\n';
	EXPECT_EQ(13, test.average_sequence_quality_for_base_.at(1).at(3).size() ) << "SRR490124-4pairs average_sequence_quality_for_base_[1][3].size() wrong for " << context << '\n';
	EXPECT_EQ(21, test.average_sequence_quality_for_base_.at(1).at(3)[25] ) << "SRR490124-4pairs average_sequence_quality_for_base_[1][3][25] wrong for " << context << '\n';
	EXPECT_EQ(13, test.average_sequence_quality_for_base_.at(1).at(3)[32] ) << "SRR490124-4pairs average_sequence_quality_for_base_[1][3][32] wrong for " << context << '\n';
}

void QualityStatsTest::TestTiles(const QualityStats &test){
	// cat <(samtools view ecoli-tiles.bam | awk '(int($2%32/16)){print "@" NR, $2, $1; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){flag=$2;name=$3}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0, name}') <(samtools view ecoli-tiles.bam | awk '(!int($2%32/16)){print $2, $10, $11, $1}') | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{len=length($2); seg=int($1%256/128); sq=0;for(pos=1;pos<=len;pos+=1){sq+=ord[substr($3,pos,1)]-33}; sq=int(sq/len+0.5);for(pos=1;pos<=len;pos+=1){print seg, substr($2,pos,1), substr($4,14,1)-1, sq, ord[substr($3,pos,1)]-33}}' | sort -k1,3 -k2,2n -k4,4n -k5,5n | uniq -c
	EXPECT_EQ(39, test.base_quality_for_sequence_per_tile_.at(0).at(0)[0][37][37]) << "base_quality_for_sequence_per_tile_[0][0][0][37][37] wrong in tile test\n";
	EXPECT_EQ(30, test.base_quality_for_sequence_per_tile_.at(0).at(0)[1][7][7]) << "base_quality_for_sequence_per_tile_[0][0][1][7][7] wrong in tile test\n";
	EXPECT_EQ(31, test.base_quality_for_sequence_per_tile_.at(0).at(1)[0][37][37]) << "base_quality_for_sequence_per_tile_[0][1][0][37][37] wrong in tile test\n";
	EXPECT_EQ(21, test.base_quality_for_sequence_per_tile_.at(0).at(1)[1][7][7]) << "base_quality_for_sequence_per_tile_[0][1][1][7][7] wrong in tile test\n";
	EXPECT_EQ(15, test.base_quality_for_sequence_per_tile_.at(0).at(2)[0][37][37]) << "base_quality_for_sequence_per_tile_[0][2][0][37][37] wrong in tile test\n";
	EXPECT_EQ(27, test.base_quality_for_sequence_per_tile_.at(0).at(2)[1][7][7]) << "base_quality_for_sequence_per_tile_[0][2][1][7][7] wrong in tile test\n";
	EXPECT_EQ(15, test.base_quality_for_sequence_per_tile_.at(0).at(3)[0][37][37]) << "base_quality_for_sequence_per_tile_[0][3][0][37][37] wrong in tile test\n";
	EXPECT_EQ(19, test.base_quality_for_sequence_per_tile_.at(0).at(3)[1][7][7]) << "base_quality_for_sequence_per_tile_[0][3][1][7][7] wrong in tile test\n";

	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_.at(0).at(0)[0][3][39]) << "base_quality_for_preceding_quality_per_tile_[0][0][0][3][39] wrong in tile test\n";
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_.at(0).at(0)[1][3][7]) << "base_quality_for_preceding_quality_per_tile_[0][0][1][3][7] wrong in tile test\n";

	EXPECT_EQ(1, test.base_quality_stats_per_tile_.at(0).at(1)[0][20][39]) << "base_quality_stats_per_tile_[0][1][0][20][39] wrong in tile test\n";
	EXPECT_EQ(1, test.base_quality_stats_per_tile_.at(1).at(3)[0][20][39]) << "base_quality_stats_per_tile_[1][3][0][20][39] wrong in tile test\n";
	SeqQualityStats<uintNucCount> stats = test.base_quality_stats_per_tile_[0][0][0][18];
	stats.Calculate();
	EXPECT_EQ(38, stats.mean_) << "base_quality_stats_per_tile_[0][0][0][18].mean_ wrong in tile test\n";
	EXPECT_EQ(1, test.base_quality_stats_per_tile_.at(0).at(0)[1][20][7]) << "base_quality_stats_per_tile_[0][0][1][20][7] wrong in tile test\n";
	EXPECT_EQ(1, test.base_quality_stats_per_tile_.at(1).at(1)[1][20][7]) << "base_quality_stats_per_tile_[1][1][20][7] wrong in tile test\n";
	stats = test.base_quality_stats_per_tile_.at(0).at(1)[1][21];
	stats.Calculate();
	EXPECT_EQ(8, stats.mean_) << "base_quality_stats_per_tile_[0][1][1][21].mean_ wrong in tile test\n";

	EXPECT_EQ(15, test.tile_quality_mean_difference_.at(0)[0][0] ) << "tile_quality_mean_difference_[0] wrong in tile test\n";
	EXPECT_EQ(-15, test.tile_quality_mean_difference_.at(0)[1][0] ) << "tile_quality_mean_difference_[0] wrong in tile test\n";
	EXPECT_EQ(-6, test.tile_quality_mean_difference_.at(0)[1][46] ) << "tile_quality_mean_difference_[0] wrong in tile test\n";
	EXPECT_EQ(-16, test.tile_quality_mean_difference_.at(0)[1][66] ) << "tile_quality_mean_difference_[0] wrong in tile test\n";
	EXPECT_EQ(15, test.tile_quality_mean_difference_.at(1)[0][0] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";
	EXPECT_EQ(-15, test.tile_quality_mean_difference_.at(1)[1][0] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";
	EXPECT_EQ(0, test.tile_quality_mean_difference_.at(1)[0][46] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";
	EXPECT_EQ(0, test.tile_quality_mean_difference_.at(1)[1][46] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";
	EXPECT_EQ(-15, test.tile_quality_mean_difference_.at(1)[0][99] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";
	EXPECT_EQ(15, test.tile_quality_mean_difference_.at(1)[1][99] ) << "tile_quality_mean_difference_[1] wrong in tile test\n";

	EXPECT_EQ(1, test.sequence_quality_mean_per_tile_.at(0)[0][37]) << "sequence_quality_mean_per_tile_[0][0][37] wrong in tile test\n";
	EXPECT_EQ(1, test.sequence_quality_mean_per_tile_.at(0)[0][38]) << "sequence_quality_mean_per_tile_[0][0][38] wrong in tile test\n";
	EXPECT_EQ(1, test.sequence_quality_mean_per_tile_.at(0)[1][7]) << "sequence_quality_mean_per_tile_[0][1][7] wrong in tile test\n";
	EXPECT_EQ(1, test.sequence_quality_mean_per_tile_.at(0)[1][9]) << "sequence_quality_mean_per_tile_[0][1][9] wrong in tile test\n";
}

void QualityStatsTest::TestDuplicates(const QualityStats &test){
	uintNucCount val_q, val_pq;
	for( uintReadLen pos = 1; pos < 50; ++pos ){
		for( uintBaseCall ref_base = 4; ref_base--; ){
			val_q = 0;
			val_pq = 0;
			for( auto num_bases : test.base_quality_stats_per_tile_reference_.at(1).at(ref_base)[0][pos] ){
				val_q += num_bases;
			}
			for( auto num_bases : test.preceding_quality_for_position_per_tile_reference_.at(1).at(ref_base)[0][pos] ){
				val_pq += num_bases;
			}
			EXPECT_EQ(val_q, val_pq) << "base_quality_stats_per_tile_reference_ and preceding_quality_for_position_per_tile_reference_ do not have the same number of bases at position " << pos << ", base " << ref_base << ". A problem with indels?\n";
		}
	}

	// The tests below were done manually and are not updated yet with verification commands
	EXPECT_EQ(1, test.preceding_quality_for_sequence_per_tile_.at(0).at(1)[0][39][40]) << "preceding_quality_for_sequence_per_tile_[0] wrong if previous quality is larger than max quality for nucleotide\n";
}

void QualityStatsTest::TestVariants(const QualityStats &test){
	// Manually modified values from TestSrr490124Equality
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(0).at(1)[0][70].size());
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(0).at(1)[0][70][2]);
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(0).at(2).at(1)[0][86][2]);
	EXPECT_EQ(99, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0].size());
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0][10].size());
	EXPECT_EQ(1, test.base_quality_stats_per_tile_per_error_reference_.at(1).at(3).at(4)[0][10][39]);

	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(0).at(1).at(2)[0][87][100]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(0).at(2).at(4)[0][97][0]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(1).at(0).at(4)[0][1][0]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_per_error_reference_.at(1).at(3).at(4)[0][99][0]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(1)[0][0][0]);
	EXPECT_EQ(0, test.error_rate_for_position_per_tile_reference_.at(0).at(2)[0][0][0]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(0)[0][99][0]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(0).at(3)[0][93][100]);
	EXPECT_EQ(0, test.error_rate_for_position_per_tile_reference_.at(1).at(0)[0][99][100]);
	EXPECT_EQ(1, test.error_rate_for_position_per_tile_reference_.at(1).at(3)[0][99][0]);

	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(0).at(1).at(2)[0][100][2]);
	EXPECT_EQ(22, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(0).at(2).at(4)[0][0][2]);
	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(1).at(0).at(4)[0][0][32]);
	EXPECT_EQ(35, test.base_quality_for_error_rate_per_tile_per_error_reference_.at(1).at(3).at(4)[0][0][2]);
	EXPECT_EQ(101, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0].size());
	EXPECT_EQ(1, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0][100].size());
	EXPECT_EQ(3, test.base_quality_for_error_rate_per_tile_reference_.at(0).at(3)[0][100][2]);
	EXPECT_EQ(101, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0].size());
	EXPECT_EQ(38, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][0].size());
	EXPECT_EQ(35, test.base_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][0][2]);

	EXPECT_EQ(51, SumVect(test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0]));
	EXPECT_EQ(54, SumVect(test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(1)[0]));
	EXPECT_EQ(4, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37].size());
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][35]);
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][36]);
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][37]);
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(0).at(1)[0][37][38]);
	EXPECT_EQ(0, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(1)[0][33].size());
	EXPECT_EQ(1, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(2)[0][33].size());
	EXPECT_EQ(2, test.base_quality_for_preceding_quality_per_tile_reference_.at(1).at(2)[0][33][36]);

	EXPECT_EQ(38, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0].size());
	EXPECT_EQ(15, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0][2]);
	EXPECT_EQ(7, test.preceding_quality_for_error_rate_per_tile_reference_.at(0).at(2)[0][0][38]);
	EXPECT_EQ(1, test.preceding_quality_for_error_rate_per_tile_reference_.at(1).at(1)[0][100].size());
	EXPECT_EQ(8, test.preceding_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][100][2]);

	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(0).at(0)[0][6][38]);
	EXPECT_EQ(93, test.preceding_quality_for_position_per_tile_reference_.at(0).at(1)[0].size());
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(0).at(2)[0][4][38]);
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(1).at(2)[0][0][1]);
	EXPECT_EQ(99, test.preceding_quality_for_position_per_tile_reference_.at(1).at(3)[0].size());
	EXPECT_EQ(1, test.preceding_quality_for_position_per_tile_reference_.at(1).at(3)[0][99][2]);

	EXPECT_EQ(2, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0].size());
	EXPECT_EQ(37, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0][19].size());
	EXPECT_EQ(6, test.base_quality_for_sequence_quality_per_tile_reference_.at(0).at(3)[0][19][2]);
	EXPECT_EQ(2, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0].size());
	EXPECT_EQ(38, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13].size());
	EXPECT_EQ(27, test.base_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13][2]);

	EXPECT_EQ(38, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(1)[0][19].size());
	EXPECT_EQ(15, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(2)[0][19][2]);
	EXPECT_EQ(7, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(0).at(1)[0][19][38]);
	EXPECT_EQ(38, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(1)[0][13].size());
	EXPECT_EQ(27, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(3)[0][13][2]);
	EXPECT_EQ(2, test.preceding_quality_for_sequence_quality_per_tile_reference_.at(1).at(1)[0][13][39]);

	EXPECT_EQ(2, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0].size());
	EXPECT_EQ(21, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(1)[0][0][19]);
	EXPECT_EQ(21, test.sequence_quality_for_error_rate_per_tile_reference_.at(0).at(2)[0][0][20]);
	EXPECT_EQ(2, test.sequence_quality_for_error_rate_per_tile_reference_.at(1).at(1)[0][100].size());
	EXPECT_EQ(5, test.sequence_quality_for_error_rate_per_tile_reference_.at(1).at(3)[0][100][13]);

	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(1)[0][0][19]);
	EXPECT_EQ(0, test.sequence_quality_for_position_per_tile_reference_.at(0).at(2)[0][0][20]);
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(0).at(0)[0][99][19]);
	EXPECT_EQ(0, test.sequence_quality_for_position_per_tile_reference_.at(0).at(0)[0][99][20]);
	EXPECT_EQ(0, test.sequence_quality_for_position_per_tile_reference_.at(1).at(0)[0][99][12]);
	EXPECT_EQ(1, test.sequence_quality_for_position_per_tile_reference_.at(1).at(3)[0][99][13]);

	EXPECT_EQ(3, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0][55][20] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(0)[0][57][19] );
	EXPECT_EQ(9, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0][46][13] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_gc_per_tile_reference_.at(1)[0][54][12] );

	EXPECT_EQ(3, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0][10][20] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(0)[0][12][19] );
	EXPECT_EQ(4, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0][11][13] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_mean_error_rate_per_tile_reference_.at(1)[0][14][12] );

	EXPECT_EQ(3, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0][55][10] );
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(0)[0][57][12] );
	EXPECT_EQ(9, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0][46][11] );
	EXPECT_EQ(1, test.mean_error_rate_for_gc_per_tile_reference_.at(1)[0][54][14] );

	EXPECT_EQ(59, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0][126][20] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(0)[0][184][19] );
	EXPECT_EQ(59, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0][184][13] );
	EXPECT_EQ(1, test.sequence_quality_mean_for_fragment_length_per_tile_reference_.at(1)[0][126][12] );

	EXPECT_EQ(59, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0][126][10] );
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(0)[0][184][12] );
	EXPECT_EQ(59, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0][184][11] );
	EXPECT_EQ(1, test.mean_error_rate_for_fragment_length_per_tile_reference_.at(1)[0][126][14] );

	EXPECT_EQ(59, test.gc_for_fragment_length_per_tile_reference_.at(0)[0].size() );
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(0)[0][126][55] );
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(0)[0][184][57] );
	EXPECT_EQ(59, test.gc_for_fragment_length_per_tile_reference_.at(1)[0].size() );
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(1)[0][184][46] );
	EXPECT_EQ(1, test.gc_for_fragment_length_per_tile_reference_.at(1)[0][126][54] );

	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_reference_.at(0), 1, 2, 1, 1, "variant test", "SRR490124-4pairs base_quality_stats_reference_[0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_reference_.at(1), 1, 10, 1, 1, "variant test", "SRR490124-4pairs base_quality_stats_reference_[1] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_reference_.at(0).at(0)[0], 0, 0, 1, 1, "variant test", "SRR490124-4pairs base_quality_stats_per_tile_reference_[0][0][0] not correctly shrunken for ");
	TestSizeAtIndex_0_1_98_99(test.base_quality_stats_per_tile_reference_.at(1).at(0)[0], 0, 1, 0, 0, "variant test", "SRR490124-4pairs base_quality_stats_per_tile_reference_[1][0][0] not correctly shrunken for ");
}

void QualityStatsTest::TestCoverage(const QualityStats &test){
	// cat <(samtools view -q10 drosophila-coverageTest.bam | awk '(int($2%32/16)){print $2; system("samtools faidx drosophila-GCF_000001215.4_cut.fna " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq -Ur) <(samtools view -q10 drosophila-coverageTest.bam | awk '(!int($2%32/16)){print $2; system("samtools faidx drosophila-GCF_000001215.4_cut.fna " $3 ":" $4 "-" $4+length($10)-1); print "+"; print $11}' | awk '{if(">"==substr($0,1,1)){print "@" substr($0,2,length($0)-1), store}else{if(2==length($1)||3==length($1)){store=$0}else{print $0}}}' | seqtk seq -U) | awk '(1==NR%4){flag=$2}(2==NR%4){seq=$0}(0==NR%4){print flag, seq, $0}' | awk 'BEGIN{split("2,3", poslist, ","); for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{for(p=1;p<=length(poslist);p+=1){pos=poslist[p];print int($1%256/128), substr($2,pos,1), 0, pos-1, ord[substr($3,pos,1)]-33}}' | sort | uniq -c | sort -k2,4 -k5,6n
	EXPECT_EQ(1, test.base_quality_stats_per_tile_reference_.at(0).at(0)[0][1].size()) << "base_quality_stats_per_tile_reference_[0][0][0][1] not correctly shrunken in coverage test\n";
	EXPECT_EQ(1, test.base_quality_stats_per_tile_reference_.at(1).at(2)[0][2].size()) << "base_quality_stats_per_tile_reference_[1][3][0][2] not correctly shrunken in coverage test\n";
}

void QualityStatsTest::TestAdapters(const QualityStats &test, const char *context, bool bwa){
	uintTileId tile_id;
	if(bwa){
		tile_id = 5; // Different sorting of reads in bam file
	}
	else{
		tile_id = 4;
	}

	// The other read pair that would count is in another tile and therefore is ignored by requiring correct segment-reverseness pairing
	// echo "count seg ref tile pos prev"; cat <(samtools view -q 10 -f 81 -F 32 ecoli-SRR490124-adapter.bam | awk '{print "@" substr($18,6,length($18)-5), $2; print substr($10,$8-$4+1,$4+length($10)-$8); print "+"; print substr($11,$8-$4+1,$4+length($11)-$8)}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 -f 161 -F 16 ecoli-SRR490124-adapter.bam | awk '{print $2, substr($18,6,length($18)-5), substr($10,1,$8+length($10)-$4), substr($11,1,$8+length($11)-$4)}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num && pos<length($3)){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(1==pos){prev=1}else{prev=ord[substr($4,pos-1,1)]-33};print int($1%256/128), ref, 4, pos-1, prev}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), 4, pos-1, ord[substr($4,pos-1,1)]-33}}' | sort -k1,2 -k4,4n -k5,5n | uniq -c
	EXPECT_EQ(77, test.preceding_quality_for_position_per_tile_reference_.at(0).at(3)[tile_id].size() ) << "for " << context;
	EXPECT_EQ(78, test.preceding_quality_for_position_per_tile_reference_.at(1).at(2)[tile_id].size() ) << "for " << context;

	EXPECT_EQ(81, test.error_rate_for_position_per_tile_per_error_reference_.at(0).at(3).at(4)[tile_id].to() ) << "for " << context;
	EXPECT_EQ(81, test.error_rate_for_position_per_tile_per_error_reference_.at(1).at(2).at(4)[tile_id].to() ) << "for " << context;

	if(bwa){
		// cat <(samtools view ecoli-SRR490124-adapter-bwa.bam -q 10 -f 16 -F 32 | awk '($4+2000>$8){print "@" $1, $2; if(-$9 < length($10)){print substr($10,length($10)+$9+1,-$9)}else{print $10}; print "+"; if(-$9 < length($11)){print substr($11,length($11)+$9+1,-$9)}else{print $11}}' | seqtk seq -r) <(samtools view ecoli-SRR490124-adapter-bwa.bam -q 10 -f 32 -F 16 | awk '($8+2000>$4){print "@" $1, $2; if($9 < length($10)){print substr($10,1,$9)}else{print $10}; print "+"; if($9 < length($11)){print substr($11,1,$9)}else{print $11}}' | seqtk seq) | awk 'BEGIN{split("2,80,81,95",poslist,",");for(n=0;n<256;n++){ord[sprintf("%c",n)]=n}}(1==NR%4){flag=$2}(0==NR%4){seg=int(flag%256/128);for(p=1;p<=length(poslist);++p){pos=poslist[p]; if(pos<length($0)){print seg, pos, ord[substr($0,pos+1,1)]-33}}}' | sort -k1,1 -k2,2n | awk 'BEGIN{pos=-1}{if($2 == pos){sum+=$3;++count}else{if(-1!=pos){print seg, pos, int(sum/count+0.5)}; seg=$1; pos=$2; sum=$3; count=1}}END{print seg, pos, int(sum/count+0.5)}'
		EXPECT_EQ(33, test.base_quality_mean_reference_.at(1)[2] ) << "for " << context;
	}
	else{
		// cat <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 16 -F 32 | awk '($4+2000>$8){print "@" $1, $2; if($4 < $8){print substr($10,$8-$4+1,length($10)+$4-$8)}else{print $10}; print "+"; if($4 < $8){print substr($11,$8-$4+1,length($11)+$4-$8)}else{print $11}}' | seqtk seq -r) <(samtools view ecoli-SRR490124-adapter.bam -q 10 -f 32 -F 16 | awk '($8+2000>$4){print "@" $1, $2; if($8 < $4){print substr($10,1,length($10)+$8-$4)}else{print $10}; print "+"; if($8 < $4){print substr($11,1,length($11)+$8-$4)}else{print $11}}' | seqtk seq) | awk 'BEGIN{split("2,80,81,95",poslist,",");for(n=0;n<256;n++){ord[sprintf("%c",n)]=n}}(1==NR%4){flag=$2}(0==NR%4){seg=int(flag%256/128);for(p=1;p<=length(poslist);++p){pos=poslist[p]; if(pos<length($0)){print seg, pos, ord[substr($0,pos+1,1)]-33}}}' | sort -k1,1 -k2,2n | awk 'BEGIN{pos=-1}{if($2 == pos){sum+=$3;++count}else{if(-1!=pos){print seg, pos, int(sum/count+0.5)}; seg=$1; pos=$2; sum=$3; count=1}}END{print seg, pos, int(sum/count+0.5)}'
		EXPECT_EQ(27, test.base_quality_mean_reference_.at(1)[2] ) << "for " << context;
	}
	EXPECT_EQ(36, test.base_quality_mean_reference_.at(0)[2] ) << "for " << context;
	EXPECT_EQ(20, test.base_quality_mean_reference_.at(0)[80] ) << "for " << context;
	EXPECT_EQ(33, test.base_quality_mean_reference_.at(0)[81] ) << "for " << context;
	EXPECT_EQ(33, test.base_quality_mean_reference_.at(0)[95] ) << "for " << context;
	EXPECT_EQ(36, test.base_quality_mean_reference_.at(1)[80] ) << "for " << context;
	EXPECT_EQ(34, test.base_quality_mean_reference_.at(1)[81] ) << "for " << context;
	EXPECT_EQ(29, test.base_quality_mean_reference_.at(1)[95] ) << "for " << context;
}
