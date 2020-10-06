#include "SurroundingTest.h"
using reseq::SurroundingTest;

#include <array>
using std::array;
#include <bitset>
using std::bitset;
#include <string>
using std::string;

//include <seqan/seq_io.h>
using seqan::DnaString;

//include "utilities.hpp"
using reseq::utilities::at;

void SurroundingTest::Register(){
	// Guarantees that library is included
}

void SurroundingTest::TearDown(){
	BasicTestClass::TearDown();
}

void SurroundingTest::TestBasics(){
	DnaString seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	Surrounding test;
	test.Set(seq, 14);
	EXPECT_EQ(0, test.BaseAt(0));
	EXPECT_EQ(1, test.BaseAt(1));
	EXPECT_EQ(2, test.BaseAt(2));
	EXPECT_EQ(3, test.BaseAt(3));
	EXPECT_EQ(1, test.BaseAt(9));
	EXPECT_EQ(2, test.BaseAt(10));
	EXPECT_EQ(3, test.BaseAt(19));
	EXPECT_EQ(0, test.BaseAt(20));
	EXPECT_EQ(1, test.BaseAt(29));
}

void SurroundingTest::TestSettersAndUpdaters(const string &test_dir){
	Reference ref;
	ASSERT_TRUE( ref.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	// cat reference-test.fa | seqtk seq | awk '(2==NR){print substr($0,length($0)-9,10) $0 substr($0,1,20)}' | awk 'BEGIN{split("0,1,2,15,497,498,499",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print i, pos, substr($0,i*10+pos+1,10)}}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}'
	string error_msg = "The function ForwardSurrounding returns wrong results\n";
	Surrounding surrounding;
	ref.ForwardSurrounding(surrounding,0,0);
	EXPECT_EQ(83, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(163795, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(909796, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,1);
	EXPECT_EQ(332, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(655183, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(493456, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,2);
	EXPECT_EQ(1330, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(523581, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(925249, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,15);
	EXPECT_EQ(1003384, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(495722, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(275383, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,497);
	EXPECT_EQ(1015809, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(313855, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(325511, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,498);
	EXPECT_EQ(917509, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(206845, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(253470, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,499);
	EXPECT_EQ(524308, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(827380, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1013881, surrounding.sur_.at(2)) << error_msg;
	// cat reference-test.fa | seqtk seq | awk '(4==NR){print substr($0,length($0)-9,10) $0 substr($0,1,20)}' | awk 'BEGIN{split("17",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print i, pos, substr($0,i*10+pos+1,10)}}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}'
	ref.ForwardSurrounding(surrounding,1,17);
	EXPECT_EQ(500825, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(144533, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(54422, surrounding.sur_.at(2)) << error_msg;

	error_msg = "The function UpdateForwardSurrounding returns wrong results\n";
	ref.ForwardSurrounding(surrounding,0,0);
	surrounding.UpdateForward(ref.ReferenceSequence(0), 1);
	EXPECT_EQ(332, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(655183, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(493456, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateForward(ref.ReferenceSequence(0), 2);
	EXPECT_EQ(1330, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(523581, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(925249, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,14);
	surrounding.UpdateForward(ref.ReferenceSequence(0), 15);
	EXPECT_EQ(1003384, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(495722, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(275383, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,0,496);
	surrounding.UpdateForward(ref.ReferenceSequence(0), 497);
	EXPECT_EQ(1015809, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(313855, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(325511, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateForward(ref.ReferenceSequence(0), 498);
	EXPECT_EQ(917509, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(206845, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(253470, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateForward(ref.ReferenceSequence(0), 499);
	EXPECT_EQ(524308, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(827380, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1013881, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurrounding(surrounding,1,16);
	surrounding.UpdateForward(ref.ReferenceSequence(1), 17);
	EXPECT_EQ(500825, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(144533, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(54422, surrounding.sur_.at(2)) << error_msg;

	// cat reference-test.fa | seqtk seq | awk '(2==NR){print substr($0,length($0)-19,20) $0 substr($0,1,10)}' | awk 'BEGIN{split("0,1,2,15,497,498,499",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print ">" i, pos; print substr($0,pos+22-i*10,10)}}}' | seqtk seq -r | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}(NR%2==1){pos=substr($0,2,length($0)-1)}(NR%2==0){mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print pos, $0, sur}'
	error_msg = "The function ReverseSurrounding returns wrong results\n";
	ref.ReverseSurrounding(surrounding,0,0);
	EXPECT_EQ(57353, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(846847, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(855350, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,1);
	EXPECT_EQ(538626, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(473855, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1000269, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,2);
	EXPECT_EQ(134656, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(642751, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1036499, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,15);
	EXPECT_EQ(613348, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(739384, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(10042, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,497);
	EXPECT_EQ(524915, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(720884, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(216468, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,498);
	EXPECT_EQ(917660, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(966653, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(54117, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,499);
	EXPECT_EQ(229415, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(241663, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(275673, surrounding.sur_.at(2)) << error_msg;
	// cat reference-test.fa | seqtk seq | awk '(4==NR){print substr($0,length($0)-19,20) $0 substr($0,1,10)}' | awk 'BEGIN{split("17",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p);for(i=0;i<3;++i){print ">" i, pos; print substr($0,pos+22-i*10,10)}}}' | seqtk seq -r | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}(NR%2==1){pos=substr($0,2,length($0)-1)}(NR%2==0){mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print pos, $0, sur}'
	ref.ReverseSurrounding(surrounding,1,17);
	EXPECT_EQ(960397, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(945044, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(626906, surrounding.sur_.at(2)) << error_msg;

	error_msg = "The function UpdateReverseSurrounding returns wrong results\n";
	auto error_msg2 = "The function RollBackReverseSurrounding returns wrong results\n";
	ref.ReverseSurrounding(surrounding,0,0);
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 1);
	EXPECT_EQ(538626, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(473855, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1000269, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 2);
	EXPECT_EQ(134656, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(642751, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(1036499, surrounding.sur_.at(2)) << error_msg;
	surrounding.RollBackReverse(ref.ReferenceSequence(0), 1);
	EXPECT_EQ(538626, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(473855, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(1000269, surrounding.sur_.at(2)) << error_msg2;
	surrounding.RollBackReverse(ref.ReferenceSequence(0), 0);
	EXPECT_EQ(57353, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(846847, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(855350, surrounding.sur_.at(2)) << error_msg2;
	ref.ReverseSurrounding(surrounding,0,14);
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 15);
	EXPECT_EQ(613348, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(739384, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(10042, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,0,16);
	surrounding.RollBackReverse(ref.ReferenceSequence(0), 15);
	EXPECT_EQ(613348, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(739384, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(10042, surrounding.sur_.at(2)) << error_msg2;
	ref.ReverseSurrounding(surrounding,0,496);
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 497);
	EXPECT_EQ(524915, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(720884, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(216468, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 498);
	EXPECT_EQ(917660, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(966653, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(54117, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateReverse(ref.ReferenceSequence(0), 499);
	EXPECT_EQ(229415, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(241663, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(275673, surrounding.sur_.at(2)) << error_msg;
	surrounding.RollBackReverse(ref.ReferenceSequence(0), 498);
	EXPECT_EQ(917660, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(966653, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(54117, surrounding.sur_.at(2)) << error_msg2;
	surrounding.RollBackReverse(ref.ReferenceSequence(0), 497);
	EXPECT_EQ(524915, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(720884, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(216468, surrounding.sur_.at(2)) << error_msg2;
	ref.ReverseSurrounding(surrounding,1,16);
	surrounding.UpdateReverse(ref.ReferenceSequence(1), 17);
	EXPECT_EQ(960397, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(945044, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(626906, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurrounding(surrounding,1,18);
	surrounding.RollBackReverse(ref.ReferenceSequence(1), 17);
	EXPECT_EQ(960397, surrounding.sur_.at(0)) << error_msg2;
	EXPECT_EQ(945044, surrounding.sur_.at(1)) << error_msg2;
	EXPECT_EQ(626906, surrounding.sur_.at(2)) << error_msg2;
}

void SurroundingTest::TestSettersAndUpdatersWithN(const string &test_dir){
	Reference ref;
	ASSERT_TRUE( ref.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );
	at(at(ref.reference_sequences_, 0), 253) = 'N';
	at(at(ref.reference_sequences_, 0), 256) = 'N';
	at(at(ref.reference_sequences_, 0), 263) = 'N';

	// cat reference-test.fa | seqtk seq | awk '(2==NR)' | awk 'BEGIN{split("250,253,254,260",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print i, pos, substr($0,i*10+pos-9,10)}}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}'
	string error_msg = "The function ForwardSurroundingWithN returns wrong results\n";
	Surrounding surrounding;
	ref.ForwardSurroundingWithN(surrounding,0,250);
	EXPECT_EQ(27546, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-7, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurroundingWithN(surrounding,0,260);
	EXPECT_EQ(-7, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(294914, surrounding.sur_.at(2)) << error_msg;
	// cat reference-test.fa | seqtk seq | awk '(4==NR)' | awk 'BEGIN{split("250",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print i, pos, substr($0,i*10+pos-9,10)}}}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{mult=1;sur=0;for(i=length($3);i>0;i-=1){sur+=mult*d[substr($3,i,1)];mult*=4}; print $1, $2, $3, sur}'
	ref.ForwardSurroundingWithN(surrounding,1,250);
	EXPECT_EQ(315976, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(222228, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(636052, surrounding.sur_.at(2)) << error_msg;
	error_msg = "The function UpdateForwardSurroundingWithN returns wrong results\n";
	ref.ForwardSurroundingWithN(surrounding,0,249);
	surrounding.UpdateForwardWithN(ref.ReferenceSequence(0), 250);
	EXPECT_EQ(27546, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-7, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(2)) << error_msg;
	for(auto pos=251; pos < 254; ++pos){
		surrounding.UpdateForwardWithN(ref.ReferenceSequence(0), pos);
	}
	EXPECT_EQ(714407, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-1, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateForwardWithN(ref.ReferenceSequence(0), 254);
	EXPECT_EQ(-10, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-10, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(525384, surrounding.sur_.at(2)) << error_msg;
	for(auto pos=255; pos < 261; ++pos){
		surrounding.UpdateForwardWithN(ref.ReferenceSequence(0), pos);
	}
	EXPECT_EQ(-7, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(294914, surrounding.sur_.at(2)) << error_msg;
	ref.ForwardSurroundingWithN(surrounding,1,249);
	surrounding.UpdateForwardWithN(ref.ReferenceSequence(1), 250);
	EXPECT_EQ(315976, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(222228, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(636052, surrounding.sur_.at(2)) << error_msg;
	// cat reference-test.fa | seqtk seq | awk '(2==NR)' | awk 'BEGIN{split("260,262,263,270",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print ">" i, pos; print substr($0,pos+2-i*10,10)}}}' | seqtk seq -r | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}(NR%2==1){pos=substr($0,2,length($0)-1)}(NR%2==0){mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print pos, $0, sur}'
	error_msg = "The function ReverseSurroundingWithN returns wrong results\n";
	ref.ReverseSurroundingWithN(surrounding,0,260);
	EXPECT_EQ(-3, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-6, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(353371, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurroundingWithN(surrounding,0,270);
	EXPECT_EQ(655351, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-3, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-6, surrounding.sur_.at(2)) << error_msg;
	// cat reference-test.fa | seqtk seq | awk '(4==NR)' | awk 'BEGIN{split("250",poslist,",")}{for(p=1;p<=length(poslist);++p){pos=poslist[p];for(i=0;i<3;++i){print ">" i, pos; print substr($0,pos+2-i*10,10)}}}' | seqtk seq -r | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}(NR%2==1){pos=substr($0,2,length($0)-1)}(NR%2==0){mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print pos, $0, sur}'
	ref.ReverseSurroundingWithN(surrounding,1,250);
	EXPECT_EQ(503704, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(1014243, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(614430, surrounding.sur_.at(2)) << error_msg;
	error_msg = "The function UpdateReverseSurroundingWithN returns wrong results\n";
	ref.ReverseSurroundingWithN(surrounding,0,259);
	surrounding.UpdateReverseWithN(ref.ReferenceSequence(0), 260);
	EXPECT_EQ(-3, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-6, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(353371, surrounding.sur_.at(2)) << error_msg;
	for(auto pos=261; pos < 263; ++pos){
		surrounding.UpdateReverseWithN(ref.ReferenceSequence(0), pos);
	}
	EXPECT_EQ(-1, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-4, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(153157, surrounding.sur_.at(2)) << error_msg;
	surrounding.UpdateReverseWithN(ref.ReferenceSequence(0), 263);
	EXPECT_EQ(913149, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-10, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-10, surrounding.sur_.at(2)) << error_msg;
	for(auto pos=264; pos < 271; ++pos){
		surrounding.UpdateReverseWithN(ref.ReferenceSequence(0), pos);
	}
	EXPECT_EQ(655351, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(-3, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(-6, surrounding.sur_.at(2)) << error_msg;
	ref.ReverseSurroundingWithN(surrounding,1,249);
	surrounding.UpdateReverseWithN(ref.ReferenceSequence(1), 250);
	EXPECT_EQ(503704, surrounding.sur_.at(0)) << error_msg;
	EXPECT_EQ(1014243, surrounding.sur_.at(1)) << error_msg;
	EXPECT_EQ(614430, surrounding.sur_.at(2)) << error_msg;
}

void SurroundingTest::TestModifiers(){
	Reference test_ref;
	resize( test_ref.reference_sequences_, 1 );
	Surrounding surrounding, comp_surrounding;

	// samtools faidx ../test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:1001-1020
	// GTTGCGAGATTTGGACGGAC

	// ChangeSurroundingBase
	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 2000);
	at(at(test_ref.reference_sequences_, 0), 1004) = 'C';
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.ChangeSurroundingBase(10, 'C');
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 2000);
	at(at(test_ref.reference_sequences_, 0), 1003) = 'C';
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.ChangeSurroundingBase(9, 'C');
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	// DeleteSurroundingBaseShiftingOnRightSide
	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 1004);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1005, 2000);
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.DeleteSurroundingBaseShiftingOnRightSide(10, at(at(species_reference_.reference_sequences_, 0), 1024));
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 1008);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1009, 2000);
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.DeleteSurroundingBaseShiftingOnRightSide(14, at(at(species_reference_.reference_sequences_, 0), 1024));
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 1015);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1016, 2000);
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.DeleteSurroundingBaseShiftingOnRightSide(21, at(at(species_reference_.reference_sequences_, 0), 1024));
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	// DeleteSurroundingBaseShiftingOnLeftSide
	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 1003);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1004, 2000);
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1003);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.DeleteSurroundingBaseShiftingOnLeftSide(9, at(at(species_reference_.reference_sequences_, 0), 993));
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 1000);
	at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1001, 2000);
	test_ref.ForwardSurrounding(comp_surrounding, 0, 1003);
	species_reference_.ForwardSurrounding(surrounding, 0, 1004);
	surrounding.DeleteSurroundingBaseShiftingOnLeftSide(6, at(at(species_reference_.reference_sequences_, 0), 993));
	EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0));
	EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1));
	EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2));

	// InsertSurroundingBasesShiftingOnRightSide
	for( uintSurPos ins_pos : array<int, 5>({9,10,18,19,28})){
		for( auto ins_bases : array<const char *, 3>({"AC","ACGT","ACGTACGTACGT"})){
			at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 994+ins_pos);
			at(test_ref.reference_sequences_, 0) += ins_bases;
			at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 994+ins_pos, 2000);
			test_ref.ForwardSurrounding(comp_surrounding, 0, 1004);
			species_reference_.ForwardSurrounding(surrounding, 0, 1004);
			surrounding.InsertSurroundingBasesShiftingOnRightSide(ins_pos, ins_bases);
			EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
		}
	}

	// InsertSurroundingBasesShiftingOnLeftSide
	for( uintSurPos ins_pos : array<int, 6>({20,19,10,9,5,0})){
		for( auto ins_bases : array<const string, 2>({"ACGT","ACGTACGTACGT"})){
			at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 995+ins_pos);
			at(test_ref.reference_sequences_, 0) += ins_bases;
			at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 995+ins_pos, 2000);
			test_ref.ForwardSurrounding(comp_surrounding, 0, 1004+ins_bases.size());
			species_reference_.ForwardSurrounding(surrounding, 0, 1004);
			surrounding.InsertSurroundingBasesShiftingOnLeftSide(ins_pos, ins_bases);
			EXPECT_EQ(comp_surrounding.sur_.at(0), surrounding.sur_.at(0)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(1), surrounding.sur_.at(1)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
			EXPECT_EQ(comp_surrounding.sur_.at(2), surrounding.sur_.at(2)) << "ins_pos: " << ins_pos << " ins_bases: " << ins_bases << std::endl;
		}
	}
}

void SurroundingTest::TestModifiersExtremCases(){
	Surrounding surrounding;
	Surrounding::intType comp_surrounding_full_t, comp_surrounding_full_a, comp_surrounding_start_a, comp_surrounding_end_a;

	comp_surrounding_full_t = Surrounding::Size()-1;
	comp_surrounding_full_a = 0;
	comp_surrounding_start_a = (Surrounding::Size()>>2)-1;
	comp_surrounding_end_a = comp_surrounding_full_t-3;

	surrounding.sur_.at(0) = comp_surrounding_full_t;
	surrounding.sur_.at(1) = comp_surrounding_full_t;
	surrounding.sur_.at(2) = comp_surrounding_full_t;

	// ChangeSurroundingBase
	for( auto block=Surrounding::kNumBlocks; block--; ){
		surrounding.ChangeSurroundingBase(block*Surrounding::kRange, 'A');
		EXPECT_EQ(comp_surrounding_start_a, surrounding.sur_.at(block)) << "Block " << block << std::endl;
		surrounding.ChangeSurroundingBase(block*Surrounding::kRange, 'T');
		EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(block)) << "Block " << block << std::endl;
		surrounding.ChangeSurroundingBase((block+1)*Surrounding::kRange-1, 'A');
		EXPECT_EQ(comp_surrounding_end_a, surrounding.sur_.at(block)) << "Block " << block << std::endl;
		surrounding.ChangeSurroundingBase((block+1)*Surrounding::kRange-1, 'T');
		EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(block)) << "Block " << block << std::endl;
	}

	// DeleteSurroundingBaseShiftingOnRightSide + InsertSurroundingBasesShiftingOnRightSide
	DnaString full_sur_sized_t, full_sur_sized_a;
	for(auto i=Surrounding::Length(); i--; ){
		full_sur_sized_t += 'T';
		full_sur_sized_a += 'A';
	}

	for( auto block=Surrounding::kNumBlocks; block--; ){
		for(auto pos : array<uintSurPos, 2>{static_cast<uintSurPos>(block*Surrounding::kRange), static_cast<uintSurPos>((block+1)*Surrounding::kRange-1)}){
			surrounding.DeleteSurroundingBaseShiftingOnRightSide(pos, 'A');
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_end_a, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnRightSide(pos, "T");
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;

			surrounding.DeleteSurroundingBaseShiftingOnRightSide(pos, 'T');
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnRightSide(pos, "T");
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;

			surrounding.InsertSurroundingBasesShiftingOnRightSide(pos, full_sur_sized_t);
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnRightSide(pos, full_sur_sized_a);
			for( auto comp_block=Surrounding::kNumBlocks; comp_block--; ){
				if(comp_block<block){
					EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
				}
				else if(comp_block==block){
					if(block*Surrounding::kRange == pos){
						EXPECT_EQ(comp_surrounding_full_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
					}
					else{
						EXPECT_EQ(comp_surrounding_end_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
					}
				}
				else{
					EXPECT_EQ(comp_surrounding_full_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
				}
			}
			surrounding.InsertSurroundingBasesShiftingOnRightSide(pos, full_sur_sized_t);
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
		}
	}

	// DeleteSurroundingBaseShiftingOnLeftSide + InsertSurroundingBasesShiftingOnLeftSide
	for( auto block=Surrounding::kNumBlocks; block--; ){
		for(auto pos : array<uintSurPos, 2>{static_cast<uintSurPos>(block*Surrounding::kRange), static_cast<uintSurPos>((block+1)*Surrounding::kRange-1)}){
			surrounding.DeleteSurroundingBaseShiftingOnLeftSide(pos, 'A');
			EXPECT_EQ(comp_surrounding_start_a, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnLeftSide(pos, "T");
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;

			surrounding.DeleteSurroundingBaseShiftingOnLeftSide(pos, 'T');
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnLeftSide(pos, "T");
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;

			surrounding.InsertSurroundingBasesShiftingOnLeftSide(pos, full_sur_sized_t);
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
			surrounding.InsertSurroundingBasesShiftingOnLeftSide(pos, full_sur_sized_a);
			for( auto comp_block=Surrounding::kNumBlocks; comp_block--; ){
				if(comp_block<block){
					EXPECT_EQ(comp_surrounding_full_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
				}
				else if(comp_block==block){
					if(block*Surrounding::kRange == pos){
						EXPECT_EQ(comp_surrounding_start_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
					}
					else{
						EXPECT_EQ(comp_surrounding_full_a, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
					}
				}
				else{
					EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(comp_block)) << "Block " << block << " CompBlock " << comp_block << " Pos " << pos << std::endl;
				}
			}
			surrounding.InsertSurroundingBasesShiftingOnLeftSide(pos, full_sur_sized_t);
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(0)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(1)) << "Block " << block << " Pos " << pos << std::endl;
			EXPECT_EQ(comp_surrounding_full_t, surrounding.sur_.at(2)) << "Block " << block << " Pos " << pos << std::endl;
		}
	}
}

void SurroundingTest::TestCombiningBias(){
	array<double, 4*Surrounding::Length()> separated_bias;
	separated_bias.fill(0.0);
	for(uintSeqLen block=0; block <= 80; block += 40){
		//ACGTTGCATA: 114252
		separated_bias.at(block+0) = 0.9;
		separated_bias.at(block+5) = 1.0;
		separated_bias.at(block+10) = 1.0;
		separated_bias.at(block+15) = 1.0;
		separated_bias.at(block+19) = 1.0;
		separated_bias.at(block+22) = 1.0;
		separated_bias.at(block+25) = 1.0;
		separated_bias.at(block+28) = 1.0;
		separated_bias.at(block+35) = 1.0;
		separated_bias.at(block+36) = 1.0;
		//ACGTT[C]CATA: 113996
		separated_bias.at(block+21) = 0.8;
	}

	SurroundingBias combined_bias;
	combined_bias.CombinePositions(separated_bias);

	for(uintSurBlockId block=0; block < Surrounding::kNumBlocks; ++block){
		EXPECT_DOUBLE_EQ(9.9, combined_bias.bias_.at(block).at(114252));
		EXPECT_DOUBLE_EQ(9.7, combined_bias.bias_.at(block).at(113996));
	}
}

void SurroundingTest::TestSeparatingBias(){
	SurroundingBias combined;

    // Fragment start
    combined.bias_.at(0).at(954112) = 1500; // TGGATTAAAA bias=1.0
    combined.bias_.at(0).at(771557) = 192; // GTTACCTGCC bias=1.0
    combined.bias_.at(0).at(756483) = 3200; // GTGAGTAAAT bias=1.0
    combined.bias_.at(0).at(594450) = 2000; // GCACAGACAG bias=0.4
    combined.bias_.at(0).at(258785) = 192; // ATTTAGTGAC bias=0.8
    combined.bias_.at(1).at(8941) = 1500; // AAAGAGTGTC bias=1.0
    combined.bias_.at(1).at(756483) = 192; // GTGAGTAAAT bias=1.0
    combined.bias_.at(1).at(787452) = 2*3200; // TAAAATTTTA bias=0.8
    combined.bias_.at(1).at(196668) = 2000; // ATAAAAATTA bias=1.0
    combined.bias_.at(1).at(461635) = 192; // CTAAGTCAAT bias=1.0
    combined.bias_.at(2).at(930377) = 1500; // TGATAGCAGC bias=1.0
    combined.bias_.at(2).at(787452) = 192; // TAAAATTTTA bias=0.8
    combined.bias_.at(2).at(1017802) = 3200; // TTGACTTAGG bias=1.0
    combined.bias_.at(2).at(297745) = 2000; // CAGAGTACAC bias=1.0
    combined.bias_.at(2).at(4080) = 192; // AAAATTTTAA bias=0.6

    // Fragment end
    combined.bias_.at(0).at(649012) = 1500; // GCTGCTATCA bias=1.0
    combined.bias_.at(0).at(787452) = 192; // TAAAATTTTA bias=0.8
    combined.bias_.at(0).at(377552) = 3200; // CCTAAGTCAA bias=1.0
    combined.bias_.at(0).at(766430) = 2000; // GTGTACTCTG bias=1.0
    combined.bias_.at(0).at(983295) = 192; // TTAAAATTTT bias=1.0
    combined.bias_.at(1).at(542591) = 1500; // GACACTCTTT bias=1.0
    combined.bias_.at(1).at(258513) = 192; // ATTTACTCAC bias=1.0
    //combined.bias_.at(1).at(787452) = 2*3200; // TAAAATTTTA bias=0.8; it is also a start surrounding
    combined.bias_.at(1).at(802803) = 2000; // TAATTTTTAT bias=1.0
    combined.bias_.at(1).at(254450) = 192; // ATTGACTTAG bias=1.0
    combined.bias_.at(2).at(1044692) = 1500; // TTTTAATCCA bias=1.0
    combined.bias_.at(2).at(674497) = 192; // GGCAGGTAAC bias=0.6
    combined.bias_.at(2).at(258513) = 3200; // ATTTACTCAC bias=1.0
    combined.bias_.at(2).at(505785) = 2000; // CTGTCTGTGC bias=1.0
    combined.bias_.at(2).at(739075) = 192; // GTCACTAAAT bias=0.8

    array<double, 4*Surrounding::Length()> separated;
    combined.SeparatePositions(separated);
    const size_t norm = Surrounding::Size()/4;
    double mean = (192 + 3200 + 192+3200+2000+1500+2000 + 1500+192+192)/4.0;
    EXPECT_DOUBLE_EQ((192-mean)/norm, separated.at(0));
    EXPECT_DOUBLE_EQ((3200-mean)/norm, separated.at(1));
    EXPECT_DOUBLE_EQ((192+3200+2000+1500+2000-mean)/norm, separated.at(2));
    EXPECT_DOUBLE_EQ((1500+192+192-mean)/norm, separated.at(3));
    mean = (1500+1500+192+3200 + 192+192 + 2000+2000 + 3200+192)/4.0;
    EXPECT_DOUBLE_EQ((1500+1500+192+3200-mean)/norm, separated.at(36));
    EXPECT_DOUBLE_EQ((192+192-mean)/norm, separated.at(37));
    EXPECT_DOUBLE_EQ((2000+2000-mean)/norm, separated.at(38));
    EXPECT_DOUBLE_EQ((3200+192-mean)/norm, separated.at(39));

    mean = (1500+2000+192+192 + 192 + 192+1500 + 2*3200+2000)/4.0;
    EXPECT_EQ((1500+2000+192+192-mean)/norm, separated.at(40));
    EXPECT_EQ((192-mean)/norm, separated.at(41));
    EXPECT_EQ((192+1500-mean)/norm, separated.at(42));
    EXPECT_EQ((2*3200+2000-mean)/norm, separated.at(43));
    mean = (2*3200+2000 + 1500+192 + 192 + 192+192+1500+2000)/4.0;
    EXPECT_EQ((2*3200+2000-mean)/norm, separated.at(76));
    EXPECT_EQ((1500+192-mean)/norm, separated.at(77));
    EXPECT_EQ((192-mean)/norm, separated.at(78));
    EXPECT_EQ((192+192+1500+2000-mean)/norm, separated.at(79));

    mean = (192+3200 + 2000+2000 + 192+192 + 1500+192+3200+1500)/4.0;
    EXPECT_EQ((192+3200-mean)/norm, separated.at(80));
    EXPECT_EQ((2000+2000-mean)/norm, separated.at(81));
    EXPECT_EQ((192+192-mean)/norm, separated.at(82));
    EXPECT_EQ((1500+192+3200+1500-mean)/norm, separated.at(83));
    mean = (192+192+1500 + 1500+2000+192+3200+2000 + 3200 + 192)/4.0;
    EXPECT_EQ((192+192+1500-mean)/norm, separated.at(116));
    EXPECT_EQ((1500+2000+192+3200+2000-mean)/norm, separated.at(117));
    EXPECT_EQ((3200-mean)/norm, separated.at(118));
    EXPECT_EQ((192-mean)/norm, separated.at(119));
}

namespace reseq{
	TEST_F(SurroundingTest, Basics){
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );

		TestBasics();
		TestSettersAndUpdaters(test_dir);
		TestSettersAndUpdatersWithN(test_dir);
	}

	TEST_F(SurroundingTest, Modifiers){
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		LoadReference(test_dir+"ecoli-GCF_000005845.2_ASM584v2_genomic.fa");

		TestModifiers();
		TestModifiersExtremCases();
	}

	TEST_F(SurroundingTest, Bias){
		TestCombiningBias();
		TestSeparatingBias();
	}
}
