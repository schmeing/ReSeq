#include "ErrorStatsTest.h"
using reseq::ErrorStatsTest;

void ErrorStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new ErrorStats() ) << "Could not allocate memory for ErrorStats object\n";
}

void ErrorStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void ErrorStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void ErrorStatsTest::TestSrr490124Equality(const ErrorStats &test, const char *context){
	// echo "count seg ref dom tile called qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, substr($3,pos,1), ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, substr($3,pos,1), ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(16, test.called_bases_by_base_quality_per_tile_.at(0).at(0).at(4)[0][0][2]) << "SRR490124-4pairs called_bases_by_base_quality_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(10, test.called_bases_by_base_quality_per_tile_.at(0).at(3).at(4)[0][3][2]) << "SRR490124-4pairs called_bases_by_base_quality_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(35, test.called_bases_by_base_quality_per_tile_.at(1).at(3).at(4)[0][3][2]) << "SRR490124-4pairs called_bases_by_base_quality_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(38, test.called_bases_by_base_quality_per_tile_.at(1).at(3).at(4)[0][3].size()) << "SRR490124-4pairs called_bases_by_base_quality_per_tile_[1] not correctly shrunken for " << context << '\n';

	// echo "count seg ref dom tile called pos"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, substr($3,pos,1), pos-1}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, substr($3,pos,1), pos-1}}' | sort -n | uniq -c | sort -k2,2n -k7,7n -k3,6
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(0).at(1).at(4)[0][1][0]) << "SRR490124-4pairs called_bases_by_position_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(0).at(2).at(4)[0][2][0]) << "SRR490124-4pairs called_bases_by_position_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(2, test.called_bases_by_position_per_tile_.at(0).at(0).at(4)[0][0][99]) << "SRR490124-4pairs called_bases_by_position_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(1).at(2).at(4)[0][2][0]) << "SRR490124-4pairs called_bases_by_position_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(1).at(3).at(4)[0][3][0]) << "SRR490124-4pairs called_bases_by_position_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(1).at(0).at(1)[0][1][99]) << "SRR490124-4pairs called_bases_by_position_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(100, test.called_bases_by_position_per_tile_.at(1).at(3).at(4)[0][3].size()) << "SRR490124-4pairs called_bases_by_position_per_tile_[1] not correctly shrunken for " << context << '\n';

	// echo "count seg ref dom tile called nerr"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; nerr=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, substr($3,pos,1), nerr}; num=0; nerr+=1}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, substr($3,pos,1), nerr}}' | sort | uniq -c | sort -k2,2n -k7,7n -k3,6
	EXPECT_EQ(1, test.called_bases_by_error_num_per_tile_.at(0).at(0).at(1)[0][1][1]) << "SRR490124-4pairs called_bases_by_error_num_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(2, test.called_bases_by_error_num_per_tile_.at(0).at(2).at(4)[0][2][12]) << "SRR490124-4pairs called_bases_by_error_num_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(25, test.called_bases_by_error_num_per_tile_.at(1).at(0).at(4)[0][0][0]) << "SRR490124-4pairs called_bases_by_error_num_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(10, test.called_bases_by_error_num_per_tile_.at(1).at(3).at(4)[0][3][6]) << "SRR490124-4pairs called_bases_by_error_num_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(12, test.called_bases_by_error_num_per_tile_.at(1).at(3).at(4)[0][3].size()) << "SRR490124-4pairs called_bases_by_error_num_per_tile_[1] not correctly shrunken for " << context << '\n';

	// echo "count seg ref dom tile called rate"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1);rate=100}else{dom="N";rate=0}; print int($1%256/128), ref, dom, 0, substr($3,pos,1), rate}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, substr($3,pos,1), 0}}' | sort | uniq -c
	EXPECT_EQ(1, test.called_bases_by_error_rate_per_tile_.at(0).at(1).at(2)[0][2][100]) << "SRR490124-4pairs called_bases_by_error_rate_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_error_rate_per_tile_.at(0).at(1).at(2)[0][2].size()) << "SRR490124-4pairs called_bases_by_error_rate_per_tile_[0] not correctly shrunken for " << context << '\n';
	EXPECT_EQ(52, test.called_bases_by_error_rate_per_tile_.at(0).at(2).at(4)[0][2][0]) << "SRR490124-4pairs called_bases_by_error_rate_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(31, test.called_bases_by_error_rate_per_tile_.at(1).at(0).at(4)[0][0][0]) << "SRR490124-4pairs called_bases_by_error_rate_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(2, test.called_bases_by_error_rate_per_tile_.at(1).at(3).at(0)[0][0][100]) << "SRR490124-4pairs called_bases_by_error_rate_per_tile_[1] not correct for " << context << '\n';

	// echo "count seg ref dom tile nerr qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; nerr=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, nerr, ord[substr($4,pos,1)]-33}; num=0; nerr+=1}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, nerr, ord[substr($4,pos,1)]-33}}' | sort | uniq -c | sort -k2,2n -k7,7n -k6,6n -k3,5
	EXPECT_EQ(1, test.error_num_by_quality_per_tile_.at(0).at(2).at(1)[0][1][2]) << "SRR490124-4pairs error_num_by_quality_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(2, test.error_num_by_quality_per_tile_.at(0).at(2).at(4)[0][12][2]) << "SRR490124-4pairs error_num_by_quality_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_num_by_quality_per_tile_.at(1).at(0).at(4)[0][0][32]) << "SRR490124-4pairs error_num_by_quality_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(10, test.error_num_by_quality_per_tile_.at(1).at(3).at(4)[0][6][2]) << "SRR490124-4pairs error_num_by_quality_per_tile_[1] not correct for " << context << '\n';

	// echo "count seg ref dom tile nerr pos"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; nerr=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 0, nerr, pos-1}; num=0; nerr+=1}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, nerr, pos-1}}' | sort | uniq -c | sort -k2,2n -k7,7n -k6,6n -k3,5
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(0).at(0).at(1)[0][1][70]) << "SRR490124-4pairs error_num_by_position_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(0).at(2).at(4)[0][12][96]) << "SRR490124-4pairs error_num_by_position_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(1).at(3).at(4)[0][0][0]) << "SRR490124-4pairs error_num_by_position_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(1).at(3).at(4)[0][11][99]) << "SRR490124-4pairs error_num_by_position_per_tile_[1] not correct for " << context << '\n';

	// echo "count seg ref dom tile nerr rate"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; nerr=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1);rate=100}else{dom="N";rate=0}; print int($1%256/128), ref, dom, 0, nerr, rate}; num=0; nerr+=1}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 0, nerr, 0}}' | sort | uniq -c | sort -k2,2n -k7,7n -k6,6n -k3,5
	EXPECT_EQ(1, test.error_num_by_error_rate_per_tile_.at(0).at(2).at(1)[0][1][100]) << "SRR490124-4pairs error_num_by_error_rate_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(2, test.error_num_by_error_rate_per_tile_.at(0).at(2).at(4)[0][12][0]) << "SRR490124-4pairs error_num_by_error_rate_per_tile_[0] not correct for " << context << '\n';
	EXPECT_EQ(2, test.error_num_by_error_rate_per_tile_.at(1).at(1).at(4)[0][1][0]) << "SRR490124-4pairs error_num_by_error_rate_per_tile_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.error_num_by_error_rate_per_tile_.at(1).at(2).at(4)[0][6][0]) << "SRR490124-4pairs error_num_by_error_rate_per_tile_[1] not correct for " << context << '\n';

	// There are no indels in these reads
	EXPECT_EQ(4, test.indel_by_indel_pos_.at(0).at(5)[0][0]) << "SRR490124-4pairs indel_by_indel_pos_ not correct for " << context << '\n';
	// echo "count prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print ">" NR; print $10}' | seqtk seq -r | awk '(0==NR%2)') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $10}') | awk '{for(pos=1;pos<length($1);pos+=1){print substr($1,pos,1)}}' | sort | uniq -c
	EXPECT_EQ(106, test.indel_by_indel_pos_.at(0).at(1)[0][0]) << "SRR490124-4pairs indel_by_indel_pos_ not correct for " << context << '\n';

	EXPECT_EQ(4, test.indel_by_position_.at(0).at(5)[0][0]) << "SRR490124-4pairs indel_by_position_ not correct for " << context << '\n';
	// echo "count prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print ">" NR; print $10}' | seqtk seq -r | awk '(0==NR%2)') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $10}') | awk '{for(pos=1;pos<length($1);pos+=1){print substr($1,pos,1), pos}}' | sort | uniq -c | sort -k3,3n
	EXPECT_EQ(1, test.indel_by_position_.at(0).at(1)[72][0]) << "SRR490124-4pairs indel_by_position_ not correct for " << context << '\n';
	EXPECT_EQ(0, test.indel_by_position_.at(0).at(1)[95][0]) << "SRR490124-4pairs indel_by_position_ not correct for " << context << '\n';

	// echo "count prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+99); print ">" NR; print $10}' | seqtk seq -r) <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+99 " | seqtk seq");print ">" NR; print $10}') | awk '(2==NR%4){gc=gsub("G","",$1)+gsub("C","",$1)}(0==NR%4){for(pos=1;pos<length($1);pos+=1){print substr($1,pos,1), gc}}' | sort | uniq -c
	EXPECT_EQ(1, test.indel_by_gc_.at(0).at(5)[46][0]) << "SRR490124-4pairs indel_by_gc_ not correct for " << context << '\n';
	EXPECT_EQ(29, test.indel_by_gc_.at(0).at(1)[46][0]) << "SRR490124-4pairs indel_by_gc_ not correct for " << context << '\n';

	EXPECT_EQ(4, test.indel_pos_by_position_.at(0).at(5)[0][0]) << "SRR490124-4pairs indel_pos_by_position_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.indel_pos_by_position_.at(0).at(1)[72][0]) << "SRR490124-4pairs indel_pos_by_position_ not correct for " << context << '\n';
	EXPECT_EQ(0, test.indel_pos_by_position_.at(0).at(1)[95][0]) << "SRR490124-4pairs indel_pos_by_position_ not correct for " << context << '\n';

	EXPECT_EQ(1, test.indel_pos_by_gc_.at(0).at(5)[46][0]) << "SRR490124-4pairs indel_pos_by_gc_ not correct for " << context << '\n';
	EXPECT_EQ(29, test.indel_pos_by_gc_.at(0).at(1)[46][0]) << "SRR490124-4pairs indel_pos_by_gc_ not correct for " << context << '\n';

	// echo "count prev"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+99); print ">" NR; print $10}' | seqtk seq -r) <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){system("samtools faidx ecoli-GCF_000005845.2_ASM584v2_genomic.fa " $3 ":" $4 "-" $4+99 " | seqtk seq");print ">" NR; print $10}') | awk '(2==NR%4){gc=gsub("G","",$1)+gsub("C","",$1)}(0==NR%4){for(pos=1;pos<length($1);pos+=1){print substr($1,pos,1), gc, pos}}' | awk '(72==$3 || 95==$3)' | sort -k3,3n
	EXPECT_EQ(1, test.gc_by_position_.at(0).at(3)[72][46]) << "SRR490124-4pairs gc_by_position_ not correct for " << context << '\n';
	EXPECT_EQ(1, test.gc_by_position_.at(0).at(0)[95][46]) << "SRR490124-4pairs gc_by_position_ not correct for " << context << '\n';

	// NM tags
	TestVectEquality({10,{1,0,1}}, test.errors_per_read_.at(0), context, "SRR490124-4pairs errors_per_read_[0]", " not correct for ");
	TestVectEquality({11,{1,0,0,0,1}}, test.errors_per_read_.at(1), context, "SRR490124-4pairs errors_per_read_[1]", " not correct for ");

	// echo "count seg ref called prev qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(1==pos){prev="N"}else{prev=substr($3,pos-1,1)}; print int($1%256/128), ref, substr($3,pos,1), prev, ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), substr($3,pos,1), substr($3,pos-1,1), ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c
	EXPECT_EQ(2, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(0).at(0).at(3)[2]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(1).at(1).at(5)[38]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(2).at(2).at(5)[38]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[0] not correct for " << context << '\n';
	EXPECT_EQ(4, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(3).at(3).at(2)[2]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[0] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(2).at(2).at(5)[39]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[1] not correct for " << context << '\n';
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(5)[39]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[1] not correct for " << context << '\n';
	EXPECT_EQ(15, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(3)[2]) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[1] not correct for " << context << '\n';
	EXPECT_EQ(38, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(3).size()) << "SRR490124-4pairs called_bases_by_base_quality_per_previous_called_base_[1] not correctly shrunken for " << context << '\n';

	// echo "count seg ref called qual"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; print int($1%256/128), ref, substr($3,pos,1), ord[substr($4,pos,1)]-33}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), substr($3,pos,1), ord[substr($4,pos,1)]-33}}' | sort -n | uniq -c  | sort -k2,4 -k5n
	EXPECT_TRUE( test.called_bases_by_base_quality_.at(0).at(3).at(4).empty() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 37, test.called_bases_by_base_quality_.at(0).at(3).at(3).size() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 10, test.called_bases_by_base_quality_.at(0).at(3).at(3)[2] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 4, test.called_bases_by_base_quality_.at(0).at(3).at(3)[38] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(2).size() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(2)[2] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(1).size() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(1)[2] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(0).size() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(0)[2] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(1).at(1).at(0).size() ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';
	EXPECT_EQ( 4, test.called_bases_by_base_quality_.at(1).at(1).at(0)[2] ) << "SRR490124-4pairs called_bases_by_base_quality_ wrong for " << context << '\n';

	// echo "count seg ref called pos"; cat <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(int($2%32/16)){print "@" substr($18,6,length($18)-5), $2; print $10; print "+";print $11}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 ecoli-SRR490124-4pairs.sam | awk '(!int($2%32/16)){print $2, substr($18,6,length($18)-5), $10, $11}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; print int($1%256/128), ref, substr($3,pos,1), pos-1}; num=0}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), substr($3,pos,1), pos-1}}' | sort -n | uniq -c | sort -k2,4 -k5n
	EXPECT_EQ( 24, test.called_bases_by_position_.at(0).at(0).at(1).size() ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[70] ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[71] ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[93] ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 100, test.called_bases_by_position_.at(1).at(3).at(3).size() ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_position_.at(1).at(3).at(3)[3] ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';
	EXPECT_EQ( 1, test.called_bases_by_position_.at(1).at(3).at(3)[99] ) << "SRR490124-4pairs called_bases_by_position_ wrong for " << context << '\n';

	EXPECT_EQ( SumVect(test.called_bases_by_error_num_per_tile_), SumVect(test.called_bases_by_base_quality_per_previous_called_base_) ) << "SRR490124-4pairs bases counted in DataStats and CoverageStats not identical " << context << '\n';
}

void ErrorStatsTest::TestDuplicates(const ErrorStats &test){
	// samtools view -q 10 -f 3 ecoli-duplicates.bam | awk '(0 != substr($1,12,1)){print int($2%256/128), substr($17,6,length($17)-5)}' | sort -k1,1n -k2,2n | uniq -c
	TestVectEquality({0,{6,2,0,4,1,1}}, test.errors_per_read_.at(0), "indels and pcr errors", "errors_per_read_[0]", " wrong with ");
	TestVectEquality({0,{9,1,1,1,0,1,1}}, test.errors_per_read_.at(1), "indels and pcr errors", "errors_per_read_[1]", " wrong with ");

	EXPECT_EQ( SumVect(test.called_bases_by_error_num_per_tile_), SumVect(test.called_bases_by_base_quality_per_previous_called_base_) ) << "Duplicates test bases counted in DataStats and CoverageStats not identical\n";
}

void ErrorStatsTest::TestVariants(const ErrorStats &test){
	// Manually modified values from TestSrr490124Equality
	EXPECT_EQ(15, test.called_bases_by_base_quality_per_tile_.at(0).at(0).at(4)[0][0][2]);
	EXPECT_EQ(10, test.called_bases_by_base_quality_per_tile_.at(0).at(3).at(4)[0][3][2]);
	EXPECT_EQ(35, test.called_bases_by_base_quality_per_tile_.at(1).at(3).at(4)[0][3][2]);
	EXPECT_EQ(38, test.called_bases_by_base_quality_per_tile_.at(1).at(3).at(4)[0][3].size());

	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(0).at(1).at(4)[0][1][0]);
	EXPECT_EQ(0, test.called_bases_by_position_per_tile_.at(0).at(2).at(4)[0][2][0]);
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(0).at(0).at(4)[0][0][99]);
	EXPECT_EQ(1, test.called_bases_by_position_per_tile_.at(1).at(2).at(4)[0][2][0]);
	EXPECT_EQ(0, test.called_bases_by_position_per_tile_.at(1).at(3).at(4)[0][3][0]);
	EXPECT_EQ(0, test.called_bases_by_position_per_tile_.at(1).at(0).at(1)[0][1][99]);
	EXPECT_EQ(99, test.called_bases_by_position_per_tile_.at(1).at(3).at(4)[0][3].size());

	EXPECT_EQ(1, test.called_bases_by_error_num_per_tile_.at(0).at(0).at(1)[0][1][1]);
	EXPECT_EQ(2, test.called_bases_by_error_num_per_tile_.at(0).at(2).at(4)[0][2][12]);
	EXPECT_EQ(25, test.called_bases_by_error_num_per_tile_.at(1).at(0).at(4)[0][0][0]);
	EXPECT_EQ(10, test.called_bases_by_error_num_per_tile_.at(1).at(3).at(4)[0][3][6]);
	EXPECT_EQ(12, test.called_bases_by_error_num_per_tile_.at(1).at(3).at(4)[0][3].size());

	EXPECT_EQ(1, test.called_bases_by_error_rate_per_tile_.at(0).at(1).at(2)[0][2][100]);
	EXPECT_EQ(1, test.called_bases_by_error_rate_per_tile_.at(0).at(1).at(2)[0][2].size());
	EXPECT_EQ(51, test.called_bases_by_error_rate_per_tile_.at(0).at(2).at(4)[0][2][0]);
	EXPECT_EQ(31, test.called_bases_by_error_rate_per_tile_.at(1).at(0).at(4)[0][0][0]);
	EXPECT_EQ(2, test.called_bases_by_error_rate_per_tile_.at(1).at(3).at(0)[0][0][100]);

	EXPECT_EQ(1, test.error_num_by_quality_per_tile_.at(0).at(2).at(1)[0][1][2]);
	EXPECT_EQ(2, test.error_num_by_quality_per_tile_.at(0).at(2).at(4)[0][12][2]);
	EXPECT_EQ(1, test.error_num_by_quality_per_tile_.at(1).at(0).at(4)[0][0][32]);
	EXPECT_EQ(10, test.error_num_by_quality_per_tile_.at(1).at(3).at(4)[0][6][2]);

	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(0).at(0).at(1)[0][1][70]);
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(0).at(2).at(4)[0][12][96]);
	EXPECT_EQ(0, test.error_num_by_position_per_tile_.at(1).at(3).at(4)[0][0][0]);
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(1).at(3).at(4)[0][11][99]);

	EXPECT_EQ(1, test.error_num_by_error_rate_per_tile_.at(0).at(2).at(1)[0][1][100]);
	EXPECT_EQ(2, test.error_num_by_error_rate_per_tile_.at(0).at(2).at(4)[0][12][0]);
	EXPECT_EQ(2, test.error_num_by_error_rate_per_tile_.at(1).at(1).at(4)[0][1][0]);
	EXPECT_EQ(1, test.error_num_by_error_rate_per_tile_.at(1).at(2).at(4)[0][6][0]);

	EXPECT_EQ(2, test.indel_by_indel_pos_.at(0).at(5)[0][0]);
	EXPECT_EQ(106, test.indel_by_indel_pos_.at(0).at(1)[0][0]);

	EXPECT_EQ(2, test.indel_by_position_.at(0).at(5)[0][0]);
	EXPECT_EQ(1, test.indel_by_position_.at(0).at(1)[72][0]);
	EXPECT_EQ(0, test.indel_by_position_.at(0).at(1)[95][0]);

	EXPECT_EQ(1, test.indel_by_gc_.at(0).at(5)[46][0]);
	EXPECT_EQ(29, test.indel_by_gc_.at(0).at(1)[46][0]);

	EXPECT_EQ(2, test.indel_pos_by_position_.at(0).at(5)[0][0]);
	EXPECT_EQ(1, test.indel_pos_by_position_.at(0).at(1)[72][0]);
	EXPECT_EQ(0, test.indel_pos_by_position_.at(0).at(1)[95][0]);

	EXPECT_EQ(1, test.indel_pos_by_gc_.at(0).at(5)[46][0]);
	EXPECT_EQ(29, test.indel_pos_by_gc_.at(0).at(1)[46][0]);

	EXPECT_EQ(1, test.gc_by_position_.at(0).at(3)[72][46]);
	EXPECT_EQ(1, test.gc_by_position_.at(0).at(0)[95][46]);

	TestVectEquality({10,{1,0,1}}, test.errors_per_read_.at(0), "adapter test", "SRR490124-4pairs errors_per_read_[0]", " not correct for ");
	TestVectEquality({11,{1,0,0,1}}, test.errors_per_read_.at(1), "adapter test", "SRR490124-4pairs errors_per_read_[1]", " not correct for ");

	EXPECT_EQ(2, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(0).at(0).at(3)[2]);
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(1).at(1).at(5)[38]);
	EXPECT_EQ(0, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(2).at(2).at(5)[38]);
	EXPECT_EQ(4, test.called_bases_by_base_quality_per_previous_called_base_.at(0).at(3).at(3).at(2)[2]);
	EXPECT_EQ(1, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(2).at(2).at(5)[39]);
	EXPECT_EQ(0, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(5)[39]);
	EXPECT_EQ(15, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(3)[2]);
	EXPECT_EQ(38, test.called_bases_by_base_quality_per_previous_called_base_.at(1).at(3).at(3).at(3).size());

	EXPECT_TRUE( test.called_bases_by_base_quality_.at(0).at(3).at(4).empty() );
	EXPECT_EQ( 37, test.called_bases_by_base_quality_.at(0).at(3).at(3).size() );
	EXPECT_EQ( 10, test.called_bases_by_base_quality_.at(0).at(3).at(3)[2] );
	EXPECT_EQ( 4, test.called_bases_by_base_quality_.at(0).at(3).at(3)[38] );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(2).size() );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(2)[2] );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(1).size() );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(1)[2] );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(0).size() );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(0).at(3).at(0)[2] );
	EXPECT_EQ( 1, test.called_bases_by_base_quality_.at(1).at(1).at(0).size() );
	EXPECT_EQ( 4, test.called_bases_by_base_quality_.at(1).at(1).at(0)[2] );

	EXPECT_EQ( 24, test.called_bases_by_position_.at(0).at(0).at(1).size() );
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[70] );
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[71] );
	EXPECT_EQ( 1, test.called_bases_by_position_.at(0).at(0).at(1)[93] );
	EXPECT_EQ( 99, test.called_bases_by_position_.at(1).at(3).at(3).size() );
	EXPECT_EQ( 1, test.called_bases_by_position_.at(1).at(3).at(3)[3] );
	EXPECT_EQ( 1, test.called_bases_by_position_.at(1).at(3).at(3)[99] );

	EXPECT_EQ( SumVect(test.called_bases_by_error_num_per_tile_), SumVect(test.called_bases_by_base_quality_per_previous_called_base_) );
}

void ErrorStatsTest::TestAdapters(const ErrorStats &test, const char *context, bool bwa){
	uintTileId tile_id;
	if(bwa){
		tile_id = 5; // Different sorting of reads in bam file
	}
	else{
		tile_id = 4;
	}

	// The other read pair that would count is in another tile and therefore is ignored by requiring correct segment-reverseness pairing
	// echo "count seg ref dom tile nerr pos"; cat <(samtools view -q 10 -f 81 -F 32 ecoli-SRR490124-adapter.bam | awk '{print "@" substr($18,6,length($18)-5), $2; print substr($10,$8-$4+1,$4+length($10)-$8); print "+"; print substr($11,$8-$4+1,$4+length($11)-$8)}' | seqtk seq -r | awk '(1==NR%4){md=substr($1,2,length($0)-1); flag=$2}(2==NR%4){seq=$0}(0==NR%4){printf("%i ", flag);num=0;mult=1;for(i=length(md);0<i;i-=1){b=substr(md,i,1);if(b ~ /^[0-9]/){num+=b*mult;mult*=10}else{base=N;if("A"==b){base="T"}; if("C"==b){base="G"}; if("G"==b){base="C"}; if("T"==b){base="A"}; printf("%i%s", num, base); num=0; mult=1}}; print num, seq, $0}') <(samtools view -q 10 -f 161 -F 16 ecoli-SRR490124-adapter.bam | awk '{print $2, substr($18,6,length($18)-5), substr($10,1,$8+length($10)-$4), substr($11,1,$8+length($11)-$4)}' ) | awk 'BEGIN{for(n=0;n<256;n++)ord[sprintf("%c",n)]=n}{pos=0; num=0; nerr=0; for(i=1;i<=length($2);i+=1){b=substr($2,i,1);if(b ~ /^[0-9]/){num=num*10+b}else{while(0 <= num && pos < length($3)){pos += 1; num -= 1; if(0>num){ref=b}else{ref=substr($3,pos,1)}; if(0>num){dom=substr($3,pos,1)}else{dom="N"}; print int($1%256/128), ref, dom, 4, nerr, pos-1}; num=0; nerr+=1}}; for(pos+=1;pos<=length($3);pos+=1){print int($1%256/128), substr($3,pos,1), "N", 4, nerr, pos-1}}' | sort | uniq -c | sort -k2,2n -k3,4 -k6,6n -k7,7n
	EXPECT_EQ(2, test.error_num_by_position_per_tile_.at(0).at(3).at(4)[tile_id].to() ) << "for " << context;
	EXPECT_EQ(81, test.error_num_by_position_per_tile_.at(0).at(3).at(4)[tile_id][1].to() ) << "for " << context;
	EXPECT_EQ(1, test.error_num_by_position_per_tile_.at(1).at(2).at(4)[tile_id].to() ) << "for " << context;
	EXPECT_EQ(81, test.error_num_by_position_per_tile_.at(1).at(2).at(4)[tile_id][0].to() ) << "for " << context;

	EXPECT_EQ( SumVect(test.called_bases_by_error_num_per_tile_), SumVect(test.called_bases_by_base_quality_per_previous_called_base_) ) << "Adapter test bases counted in DataStats and CoverageStats not identical for " << context;

	// No InDels should be detected
	for( uintBaseCall call = 6; call--; ){
		if(4 == call){
			// We don't have N's in this testset
			EXPECT_EQ(0, test.indel_by_indel_pos_.at(0).at(call).to()) << " Call " << call << " for " << context;
		}
		else{
			EXPECT_EQ(1, test.indel_by_indel_pos_.at(0).at(call).to()) << " Call " << call << " for " << context;
			EXPECT_EQ(1, test.indel_by_indel_pos_.at(0).at(call)[0].to()) << " Call " << call << " for " << context;
		}
		EXPECT_EQ(0, test.indel_by_indel_pos_.at(1).at(call).size()) << " Call " << call << " for " << context;
	}
}
