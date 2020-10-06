#include "ReferenceTest.h"
using reseq::ReferenceTest;

#include <array>
using std::array;
#include <set>
using std::set;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include <seqan/seq_io.h>
using seqan::DnaString;

#include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::VectorAtomic;
using reseq::utilities::Complement;
using reseq::utilities::IsN;
#include "Surrounding.h"
using reseq::Surrounding;

#include "CMakeConfig.h"

void ReferenceTest::Register(){
	// Guarantees that library is included
}

void ReferenceTest::TestVariantClass(){
	Reference::Variant test(0, "", {7});
	EXPECT_TRUE(test.InAllele(2));
	EXPECT_TRUE(test.InAllele(1));
	EXPECT_TRUE(test.InAllele(0));
	EXPECT_EQ(0, test.FirstAllele());

	test.allele_.at(0) = 6;
	EXPECT_TRUE(test.InAllele(2));
	EXPECT_TRUE(test.InAllele(1));
	EXPECT_FALSE(test.InAllele(0));
	EXPECT_EQ(1, test.FirstAllele());

	test.allele_.at(0) = 5;
	EXPECT_TRUE(test.InAllele(2));
	EXPECT_FALSE(test.InAllele(1));
	EXPECT_TRUE(test.InAllele(0));
	EXPECT_EQ(0, test.FirstAllele());

	test.allele_.at(0) = 4;
	EXPECT_TRUE(test.InAllele(2));
	EXPECT_FALSE(test.InAllele(1));
	EXPECT_FALSE(test.InAllele(0));
	EXPECT_EQ(2, test.FirstAllele());

	test.allele_.at(0) = 3;
	EXPECT_FALSE(test.InAllele(2));
	EXPECT_TRUE(test.InAllele(1));
	EXPECT_TRUE(test.InAllele(0));
	EXPECT_EQ(0, test.FirstAllele());

	test.allele_.at(0) = 2;
	EXPECT_FALSE(test.InAllele(2));
	EXPECT_TRUE(test.InAllele(1));
	EXPECT_FALSE(test.InAllele(0));
	EXPECT_EQ(1, test.FirstAllele());

	test.allele_.at(0) = 1;
	EXPECT_FALSE(test.InAllele(2));
	EXPECT_FALSE(test.InAllele(1));
	EXPECT_TRUE(test.InAllele(0));
	EXPECT_EQ(0, test.FirstAllele());

	test.allele_.at(1) = 1;
	EXPECT_FALSE(test.InAllele(2));
	EXPECT_FALSE(test.InAllele(1));
	EXPECT_TRUE(test.InAllele(0));
	EXPECT_TRUE(test.InAllele(64));
	EXPECT_EQ(0, test.FirstAllele());

	test.allele_.at(0) = 0;
	EXPECT_FALSE(test.InAllele(0));
	EXPECT_TRUE(test.InAllele(64));
	EXPECT_EQ(64, test.FirstAllele());
}

void ReferenceTest::TestInsertVariant(){
	ref_.variants_.resize(1);
	ref_.InsertVariant(0, 0, DnaString("A"), {1});
	ref_.InsertVariant(0, 1, DnaString("A"), {2});
	ref_.InsertVariant(0, 1, DnaString("ACT"), {1});
	ref_.InsertVariant(0, 1, DnaString("C"), {4});
	ref_.InsertVariant(0, 1, DnaString(""), {8});
	ref_.InsertVariant(0, 1, DnaString("C"), {16});
	ref_.InsertVariant(0, 1, DnaString(""), {32});
	ref_.InsertVariant(0, 1, DnaString("ACT"), {64});
	ref_.InsertVariant(0, 1, DnaString("TG"), {128});

	const vector<Reference::Variant> &vars(ref_.variants_.at(0));
	EXPECT_EQ(6, vars.size());
	for(intVariantId nvar=0; nvar < vars.size(); ++nvar){
		switch(nvar){
		case 0:
			EXPECT_EQ(0, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_EQ(1, vars.at(nvar).allele_.at(0));
			break;
		case 1:
			EXPECT_EQ(1, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "" );
			EXPECT_EQ(40, vars.at(nvar).allele_.at(0));
			break;
		case 2:
			EXPECT_EQ(1, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_EQ(2, vars.at(nvar).allele_.at(0));
			break;
		case 3:
			EXPECT_EQ(1, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "C" );
			EXPECT_EQ(20, vars.at(nvar).allele_.at(0));
			break;
		case 4:
			EXPECT_EQ(1, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "TG" );
			EXPECT_EQ(128, vars.at(nvar).allele_.at(0));
			break;
		case 5:
			EXPECT_EQ(1, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "ACT" );
			EXPECT_EQ(65, vars.at(nvar).allele_.at(0));
			break;
		}
	}

	ref_.variants_.clear();
}

void ReferenceTest::TestVariationLoading(){
	// freebayes -f GCF_000005845.2_ASM584v2_genomic.fna SRR490124.bam | awk 'BEGIN{pos[1]=1;pos[17]=1;pos[21]=1;pos[11369]=1;pos[953166]=1;pos[3192438]=1;pos[3424236]=1;pos[3424237]=1}("#" == substr($0,1,1) || $2 in pos)' > test-var.vcf
	// Manually changed alternative for pos 3192438 from GT to G to create deletion
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.PrepareVariantFile(test_dir+"test-var.vcf") );
	ASSERT_TRUE( ref_.ReadFirstVariants() );

	EXPECT_EQ(2, ref_.num_alleles_);
	ASSERT_EQ(1, ref_.variants_.size());
	EXPECT_EQ(13, ref_.variants_.at(0).size());

	const vector<Reference::Variant> &vars(ref_.variants_.at(0));
	for(intVariantId nvar=0; nvar < vars.size(); ++nvar){
		switch(nvar){
		case 0:
			EXPECT_EQ(2, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "T" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 1:
			EXPECT_EQ(3, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 2:
			EXPECT_EQ(4, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 3:
			EXPECT_EQ(5, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "G" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 4:
			EXPECT_EQ(7, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 5:
			EXPECT_EQ(8, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "TTTTTCAGCTTTTCA" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 6:
			EXPECT_EQ(11368, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "T" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 7:
			EXPECT_EQ(11370, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "T" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 8:
			EXPECT_EQ(953165, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "G" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 9:
			EXPECT_EQ(3192437, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "G" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 10:
			EXPECT_EQ(3192438, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "" );
			EXPECT_FALSE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 11:
			EXPECT_EQ(3424235, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "C" );
			EXPECT_TRUE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		case 12:
			EXPECT_EQ(3424236, vars.at(nvar).position_);
			EXPECT_TRUE( vars.at(nvar).var_seq_ == "A" );
			EXPECT_TRUE( vars.at(nvar).InAllele(0) );
			EXPECT_TRUE( vars.at(nvar).InAllele(1) );
			break;
		}
	}

	ref_.ClearVariants(1);
	EXPECT_EQ(0, ref_.variants_.at(0).size());

	ref_.ClearAllVariants();
	EXPECT_EQ(0, ref_.variants_.size());
}

void ReferenceTest::TestVariationPositionLoading(){
	// freebayes -f GCF_000005845.2_ASM584v2_genomic.fna SRR490124.bam | awk 'BEGIN{pos[1]=1;pos[17]=1;pos[21]=1;pos[11369]=1;pos[953166]=1;pos[3192438]=1;pos[3424236]=1;pos[3424237]=1}("#" == substr($0,1,1) || $2 in pos)' > test-var.vcf
	// Manually changed alternative for pos 3192438 from GT to G to create deletion
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.PrepareVariantFile(test_dir+"test-var.vcf") );
	ASSERT_TRUE( ref_.ReadFirstVariantPositions() );

	ASSERT_EQ(1, ref_.variant_positions_.size());

	TestVectorEquality({0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 20, 11368, 11369, 11370, 953165, 3192437, 3192438, 3424235, 3424236}, ref_.variant_positions_.at(0), "ref_.variant_positions_.at(0)", "is wrong", " for VariationPositionLoading.");

	ref_.ClearVariantPositions(1);
	EXPECT_EQ(0, ref_.variant_positions_.at(0).size());

	ref_.ClearAllVariantPositions();
	EXPECT_EQ(0, ref_.variant_positions_.size());
}

void ReferenceTest::TestLoadingAndAccess(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-special-chars.fa").c_str() ) );
	for( auto base : ref_[0] ){
		EXPECT_TRUE(base == 'N') << "Special characters not handled properly\n";
	}

	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	ASSERT_EQ(2, ref_.NumberSequences()) << "The loaded test reference does not have the right number of sequences\n";
	EXPECT_TRUE("NC_000913.3_1-500 bla" == ref_.ReferenceId(0)) << ref_.ReferenceId(0) << "\nThe first loaded test reference sequence does not have the correct name.\n";
	EXPECT_TRUE("NC_000913.3_10000-10500 blub" == ref_.ReferenceId(1)) << ref_.ReferenceId(1) << "\nThe second loaded test reference sequence does not have the correct name\n";
	EXPECT_TRUE("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCAT" == ref_.ReferenceSequence(0)) << "The first loaded test reference sequence is not correct\n";
	EXPECT_TRUE("GATTGCGCTGGCACCGCAGATCAGCCCAATCCAGCCGGCAAAGTGGATGATTGCGGCGTTACCGGCAATGTTACCGATCGCCAGCAGGGCAAACAGCACGGTCAGGCTAAAGAAAACGAATTGCAGAACGCGTGCGCCTTTCAGCGTGCCGAAGAACATAAACAGCGTAAATACGCCCCACAGACCCAGGTAGACACCAAGGAACTGTGCATTTGGCGCATCGGTCAGACCCAGTTTCGGCATCAGCAGAATCGCAACCAGCGTCAGCCAGAAAGAACCGTAAGAGGTGAATGCGGTTAAACCGAAAGTGTTGCCTTTTTTGTACTCCAGCAGACCAGCAAAAATTTGCGCGATGCCGCCGTAGAAAATGCCCATGGCAAGAATAATACCGTCCAGAGCGAAATAACCCACGTTGTGCAGGTTAAGCAGAATGGTGGTCATGCCGAAGCCCATCAGGCCCAGCGGTGCCGGATTAGCCAACTTAGTGTTGCCCATAATTCC" == ref_[1]) << "The second loaded test reference sequence is not correct\n";
	EXPECT_EQ(500, ref_.SequenceLength(0)) << "SequenceLength function does not work\n";
	EXPECT_EQ(501, ref_.SequenceLength(1)) << "SequenceLength function does not work\n";

	DnaString test_string;
	ref_.ReferenceSequence(test_string, 0, 0, 10);
	EXPECT_TRUE("AGCTTTTCAT" == test_string) << test_string << "\nReferenceSequence function does return correct result\n";
	EXPECT_EQ(10, length(test_string)) << "ReferenceSequence function does return correct size\n";
	ref_.ReferenceSequence(test_string, 0, 500, 10, true);
	EXPECT_TRUE("ATGGTTTTTT" == test_string) << test_string << "\nReferenceSequence function does return correct reversed result\n";
	EXPECT_EQ(10, length(test_string)) << "ReferenceSequence function does return correct size if reversed\n";

	vector<Reference::Variant> test_variants;
	test_variants.push_back({2, "", {2}}); // Allele 1
	test_variants.push_back({4, "TAG", {3}}); // Allele 0 + 1
	test_variants.push_back({9, "C", {1}}); // Allele 0
	ref_.ReferenceSequence(test_string, 0, 0, 12, false, test_variants, {0, 0}, 0);
	EXPECT_TRUE("AGCTTAGTTCAC" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 0, 11, false, test_variants, {0, 0}, 1);
	EXPECT_TRUE("AGTTAGTTCAT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 4, 8, false, test_variants, {1, 0}, 1);
	EXPECT_TRUE("TAGTTCAT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 4, 7, false, test_variants, {1, 1}, 1);
	EXPECT_TRUE("AGTTCAT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 4, 6, false, test_variants, {1, 2}, 1);
	EXPECT_TRUE("GTTCAT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 10, 12, true, test_variants, {2, 0}, 0);
	EXPECT_TRUE("GTGAACTAAGCT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 10, 11, true, test_variants, {2, 0}, 1);
	EXPECT_TRUE("ATGAACTAACT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 7, true, test_variants, {1, 0}, 0);
	EXPECT_TRUE("CTAAGCT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 6, true, test_variants, {1, 2}, 0);
	EXPECT_TRUE("TAAGCT" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 5, true, test_variants, {1, 1}, 0);
	EXPECT_TRUE("AAGCT" == test_string) << test_string << "\n";

	test_variants.push_back({499, "TAG", {3}}); // Allele 0 + 1
	ref_.ReferenceSequence(test_string, 0, 500, 12, true, test_variants, {3, 0}, 0);
	EXPECT_TRUE("CTATGGTTTTTT" == test_string) << test_string << "\n";

	// Test what happens if given fragment length is shorter than variant length
	ref_.ReferenceSequence(test_string, 0, 4, 1, false, test_variants, {1, 0}, 1);
	EXPECT_TRUE("T" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 4, 1, false, test_variants, {1, 1}, 1);
	EXPECT_TRUE("A" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 4, 1, false, test_variants, {1, 2}, 1);
	EXPECT_TRUE("G" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 1, true, test_variants, {1, 0}, 0);
	EXPECT_TRUE("C" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 1, true, test_variants, {1, 2}, 0);
	EXPECT_TRUE("T" == test_string) << test_string << "\n";
	ref_.ReferenceSequence(test_string, 0, 5, 1, true, test_variants, {1, 1}, 0);
	EXPECT_TRUE("A" == test_string) << test_string << "\n";

	string error_msg = "\nThe extraction of only the first part by cutting at the first space does not work\n";
	EXPECT_TRUE(ref_.ReferenceIdFirstPart(0) == "NC_000913.3_1-500") << ref_.ReferenceIdFirstPart(0) << error_msg;
	EXPECT_TRUE(ref_.ReferenceIdFirstPart(1) == "NC_000913.3_10000-10500") << ref_.ReferenceIdFirstPart(1) << error_msg;
}

void ReferenceTest::TestGC(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	string error_msg = "The function GCContentAbsolut returns wrong results\n";
	EXPECT_EQ(2, ref_.GCContentAbsolut(0,0,7)) << error_msg;
	EXPECT_EQ(4, ref_.GCContentAbsolut(1,22,27)) << error_msg;

	error_msg = "The function GCContent returns wrong results\n";
	EXPECT_EQ(29, ref_.GCContent(0,0,7)) << error_msg;
	EXPECT_EQ(80, ref_.GCContent(1,22,27)) << error_msg;

	error_msg = "The function UpdateGC returns wrong results\n";
	uintSeqLen n_count(0);
	uintSeqLen gc = ref_.GCContentAbsolut(n_count,0,0,7);
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(0), 0, 7);
	EXPECT_EQ(3, gc) << error_msg;
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(0), 1, 8);
	EXPECT_EQ(2, gc) << error_msg;
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(0), 2, 9);
	EXPECT_EQ(1, gc) << error_msg;
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(0), 3, 10);
	EXPECT_EQ(1, gc) << error_msg;
	gc = ref_.GCContentAbsolut(n_count, 1,4,10);
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(1), 4, 10);
	EXPECT_EQ(5, gc) << error_msg;
	gc = ref_.GCContentAbsolut(n_count, 1,22,27);
	ref_.UpdateGC(gc, n_count, ref_.ReferenceSequence(1), 22, 27);
	EXPECT_EQ(4, gc) << error_msg;
	EXPECT_EQ(0, n_count) << "In one of the GC tests above somehow a non-existing N was detected.";
}

void ReferenceTest::TestSumBias(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	string error_msg = "The function SumBias returns wrong results\n";
	// cat <(seqtk seq reference-test.fa | awk '(2==NR)') <(seqtk seq -r reference-test.fa | awk '(2==NR)') | awk '{for(i=21;i<=length($0)-50; i+=1){print substr($0,i-10,30), NR, i-1}}' | awk '{print ">"$1, $2, $3; print $1}' | seqtk seq -r | awk '(1==NR%2){full=substr($1,2,length($1)-1); strand=$2-1; pos=$3}(0==NR%2){print strand, pos, substr(full,11,10), full, $0}' | awk '{print $1, $2, gsub(/[GC]/,"",$3), substr($4,1,10), substr($4, 11, 10), substr($4, 21, 10) , substr($5,1,10), substr($5, 11, 10), substr( $5, 21, 10)}' | awk '(0==$1 && (27 == $2 || 50 == $2 || 90 == $2 || 100 == $2 || 162 == $2 || 242 == $2) || 1==$1 && 381==$2)' | awk '{print $1, $2, $3; print $4; print $5; print $6; print $7; print $8; print $9}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{if(NR%7==1){print $0}else{mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print $0, sur}}'
	// cat <(seqtk seq reference-test.fa | awk '(2==NR)') <(seqtk seq -r reference-test.fa | awk '(2==NR)') | awk '{print substr($0,length($0)-9,length($0)) $0 substr($0,1,10)}' | awk '{for(i=11;i<=length($0)-20; i+=1){print substr($0,i-10,30), NR, i-11}}' | awk '{print ">"$1, $2, $3; print $1}' | seqtk seq -r | awk '(1==NR%2){full=substr($1,2,length($1)-1); strand=$2-1; pos=$3}(0==NR%2){print strand, pos, substr(full,11,10), full, $0}' | awk '{print $1, $2, gsub(/[GC]/,"",$3), substr($4,1,10), substr($4, 11, 10), substr($4, 21, 10) , substr($5,1,10), substr($5, 11, 10), substr( $5, 21, 10)}' | awk 'BEGIN{gc[0] = 1.0; gc[3] = 0.1; gc[4] = 0.3; gc[9] = 0.2; base=10; high=base+0.5; normal=base+0.0; low=base-0.5; very_low=base-1; sur["GCAACGGGCA"]=normal; sur["TGGATTAAAA"]=normal; sur["GTTACCTGCC"]=normal; sur["GTGAGTAAAT"]=normal; sur["GCACAGACAG"]=very_low; sur["CCACAGGTAA"]=normal; sur["ATTTAGTGAC"]=high; sur["ATATGTCTCT"]=normal; sur["AAAGAGTGTC"]=normal; sur["GTGAGTAAAT"]=normal; sur["TAAAATTTTA"]=high; sur["ATAAAAATTA"]=normal; sur["CGGTGCGGGC"]=normal; sur["CTAAGTCAAT"]=normal; sur["GTGTGGATTA"]=normal; sur["TGATAGCAGC"]=normal; sur["TAAAATTTTA"]=high; sur["TTGACTTAGG"]=normal; sur["CAGAGTACAC"]=normal; sur["TGACGCGTAC"]=normal; sur["AAAATTTTAA"]=low; sur["TAATCCACAC"]=normal; sur["GCTGCTATCA"]=normal; sur["TAAAATTTTA"]=high; sur["CCTAAGTCAA"]=normal; sur["GTGTACTCTG"]=normal; sur["GTACGCGTCA"]=normal; sur["TTAAAATTTT"]=normal; sur["AGAGACATAT"]=normal; sur["GACACTCTTT"]=normal; sur["ATTTACTCAC"]=normal; sur["TAATTTTTAT"]=normal; sur["GCCCGCACCG"]=normal; sur["ATTGACTTAG"]=normal; sur["TGCCCGTTGC"]=normal; sur["TTTTAATCCA"]=normal; sur["GGCAGGTAAC"]=low; sur["ATTTACTCAC"]=normal; sur["CTGTCTGTGC"]=normal; sur["TTACCTGTGG"]=normal; sur["GTCACTAAAT"]=high}{start=0; for(i=4; i <= 6; i++){if(sur[$i] > 5){start+=sur[$i]-base}else{start-=1000}}; end=0; for(i=7; i <= 9; i++){if(sur[$i] > 5){end+=sur[$i]-base}else{end-=1000}}; print $0, 0.0+gc[$3], start, end}' | awk '{print $0, 2/(1+exp(-$11)), 2/(1+exp(-$12))}' 2>/dev/null | awk '{print $0, $10*$13*$14}' | awk '{print $0, 0.5*$15}' | awk 'BEGIN{sum=0; max=0}{sum += $16; if($16 > max){max=$16}}END{print "sum:", sum, " max: ", max}'
	Vect<double> gc_bias;
	gc_bias[30] = 0.1; // Filtered because it is at ref seq ends for stats creation: So the SumBias function for simulations tested here counts it, GetFragmentSites completely excludes it
	gc_bias[40] = 0.3; // Right at the beginning of the accepted sequence, so should be counted in SumBias and the first one entered in GetFragmentSites
	// gc_bias[30] = 0.1;
	gc_bias[0] = 1.0;
	// gc_bias[0] = 1.0;
	gc_bias[90] = 0.2; // N's will block this one for GetFragmentSites(So placeholder expected)
	// gc_bias[30] = 0.1; // Reverse strand

	SurroundingBias sur_bias;
	for( auto &block : sur_bias.bias_){
		block.clear();
		block.resize(Surrounding::Size(), -1000.0); // Set to big negative number, so that bias will be very close to zero
	}

	double val_high = 0.5;
	double val_normal = 0.0;
	double val_low = -0.5;
	double val_very_low = -1.0;
	// Fragment starts
	sur_bias.bias_.at(0).at(591524) = val_normal; // GCAACGGGCA
	sur_bias.bias_.at(0).at(954112) = val_normal; // TGGATTAAAA
	sur_bias.bias_.at(0).at(771557) = val_normal; // GTTACCTGCC
	sur_bias.bias_.at(0).at(756483) = val_normal; // GTGAGTAAAT
	sur_bias.bias_.at(0).at(594450) = val_very_low; // GCACAGACAG
	sur_bias.bias_.at(0).at(332464) = val_normal; // CCACAGGTAA
	sur_bias.bias_.at(0).at(258785) = val_high; // ATTTAGTGAC
	sur_bias.bias_.at(1).at(211831) = val_normal; // ATATGTCTCT
	sur_bias.bias_.at(1).at(8941) = val_normal; // AAAGAGTGTC
	sur_bias.bias_.at(1).at(756483) = val_normal; // GTGAGTAAAT
	sur_bias.bias_.at(1).at(787452) = val_high; // TAAAATTTTA
	sur_bias.bias_.at(1).at(196668) = val_normal; // ATAAAAATTA
	sur_bias.bias_.at(1).at(440745) = val_normal; // CGGTGCGGGC
	sur_bias.bias_.at(1).at(461635) = val_normal; // CTAAGTCAAT
	sur_bias.bias_.at(2).at(768572) = val_normal; // GTGTGGATTA
	sur_bias.bias_.at(2).at(930377) = val_normal; // TGATAGCAGC
	sur_bias.bias_.at(2).at(787452) = val_high; // TAAAATTTTA
	sur_bias.bias_.at(2).at(1017802) = val_normal; // TTGACTTAGG
	sur_bias.bias_.at(2).at(297745) = val_normal; // CAGAGTACAC
	sur_bias.bias_.at(2).at(924081) = val_normal; // TGACGCGTAC
	sur_bias.bias_.at(2).at(4080) = val_low; // AAAATTTTAA

	// Fragment ends
	sur_bias.bias_.at(0).at(800017) = val_normal; // TAATCCACAC
	sur_bias.bias_.at(0).at(649012) = val_normal; // GCTGCTATCA
	sur_bias.bias_.at(0).at(787452) = val_high; // TAAAATTTTA
	sur_bias.bias_.at(0).at(377552) = val_normal; // CCTAAGTCAA
	sur_bias.bias_.at(0).at(766430) = val_normal; // GTGTACTCTG
	sur_bias.bias_.at(0).at(727476) = val_normal; // GTACGCGTCA
	sur_bias.bias_.at(0).at(983295) = val_normal; // TTAAAATTTT
	sur_bias.bias_.at(1).at(139571) = val_normal; // AGAGACATAT
	sur_bias.bias_.at(1).at(542591) = val_normal; // GACACTCTTT
	sur_bias.bias_.at(1).at(258513) = val_normal; // ATTTACTCAC
	//sur_bias.bias_.at(1).at(787452) = val_high; // TAAAATTTTA it is also a start surrounding
	sur_bias.bias_.at(1).at(802803) = val_normal; // TAATTTTTAT
	sur_bias.bias_.at(1).at(612630) = val_normal; // GCCCGCACCG
	sur_bias.bias_.at(1).at(254450) = val_normal; // ATTGACTTAG
	sur_bias.bias_.at(2).at(939769) = val_normal; // TGCCCGTTGC
	sur_bias.bias_.at(2).at(1044692) = val_normal; // TTTTAATCCA
	sur_bias.bias_.at(2).at(674497) = val_low; // GGCAGGTAAC
	sur_bias.bias_.at(2).at(258513) = val_normal; // ATTTACTCAC
	sur_bias.bias_.at(2).at(505785) = val_normal; // CTGTCTGTGC
	sur_bias.bias_.at(2).at(989114) = val_normal; // TTACCTGTGG
	sur_bias.bias_.at(2).at(739075) = val_high; // GTCACTAAAT

	double max_bias(0.0);
	EXPECT_NEAR(2.9367, 2*ref_.SumBias(max_bias, 0, 10, 0.5, gc_bias, sur_bias), 0.0001 ) << error_msg; // The factor of 2 is because SumBias only calculates for one strand
	EXPECT_NEAR(0.774915, max_bias, 0.00001); // Rounding errors during summation using floating point precision in awk gives lower precision than the value stated suggests
}

void ReferenceTest::TestGetFragmentSites(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );
	at(at(ref_.reference_sequences_, 0), 253) = 'N';
	at(at(ref_.reference_sequences_, 0), 256) = 'N';
	at(at(ref_.reference_sequences_, 0), 263) = 'N';
	ref_.PrepareExclusionRegions();
	ref_.ObtainExclusionRegions(ref_.NumberSequences(), 100);

	string error_msg = "The function GetFragmentSites returns wrong results\n";
	std::vector<FragmentSite> frag_sites;
	ref_.GetFragmentSites(frag_sites, 0, 10, 0, 1000);

	// seqtk seq reference-test.fa | awk '(2==NR){for(i=21;i<=length($0)-50; i+=1){print substr($0,i-10,30), i-1}}' | awk '{print ">"$1, $2; print $1}' | seqtk seq -r | awk '(1==NR%2){full=substr($1,2,length($1)-1); pos=$2}(0==NR%2){print pos, substr(full,11,10), full, $0}' | awk '{print $1, gsub(/[GC]/,"",$2), substr($3,1,10), substr($3, 11, 10), substr($3, 21, 10) , substr($4,1,10), substr($4, 11, 10), substr( $4, 21, 10)}' | awk '(50 == $1 || 100 == $1 || 242 == $1 || 400 == $1)' | awk '{print $1, $2; print $3; print $4; print $5; print $6; print $7; print $8}' | awk 'BEGIN{d["A"]=0;d["C"]=1;d["G"]=2;d["T"]=3}{if(NR%7==1){print $0}else{mult=1;sur=0;for(i=length($0);i>0;i-=1){sur+=mult*d[substr($0,i,1)];mult*=4}; print $0, sur}}'
	EXPECT_EQ(391, frag_sites.size() ) << error_msg; // 500(SequenceLength) -10(fragmentLength)+1 -2*50(min_dist_to_ref_seq_ends_)
	EXPECT_EQ(40, frag_sites.at(0).gc_) << error_msg;
	EXPECT_EQ(954112, frag_sites.at(0).start_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(8941, frag_sites.at(0).start_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(930377, frag_sites.at(0).start_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(649012, frag_sites.at(0).end_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(542591, frag_sites.at(0).end_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(1044692, frag_sites.at(0).end_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(0, frag_sites.at(0).count_forward_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(0).count_reverse_) << error_msg;
	EXPECT_EQ(1.0, frag_sites.at(0).bias_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(50).gc_) << error_msg;
	EXPECT_EQ(756483, frag_sites.at(50).start_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(787452, frag_sites.at(50).start_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(1017802, frag_sites.at(50).start_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(377552, frag_sites.at(50).end_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(787452, frag_sites.at(50).end_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(258513, frag_sites.at(50).end_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(0, frag_sites.at(50).count_forward_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(50).count_reverse_) << error_msg;
	EXPECT_EQ(1.0, frag_sites.at(50).bias_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(192).gc_) << error_msg; // N's are blocking this one (so we just have a placeholder here)
	EXPECT_EQ(-1, frag_sites.at(192).start_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(-1, frag_sites.at(192).start_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(-1, frag_sites.at(192).start_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(-1, frag_sites.at(192).end_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(-1, frag_sites.at(192).end_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(-1, frag_sites.at(192).end_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(0, frag_sites.at(192).count_forward_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(192).count_reverse_) << error_msg;
	EXPECT_EQ(0.0, frag_sites.at(192).bias_) << error_msg;
	EXPECT_EQ(30, frag_sites.at(350).gc_) << error_msg;
	EXPECT_EQ(454550, frag_sites.at(350).start_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(212456, frag_sites.at(350).start_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(37093, frag_sites.at(350).start_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(675743, frag_sites.at(350).end_surrounding_.sur_.at(0)) << error_msg;
	EXPECT_EQ(870451, frag_sites.at(350).end_surrounding_.sur_.at(1)) << error_msg;
	EXPECT_EQ(430150, frag_sites.at(350).end_surrounding_.sur_.at(2)) << error_msg;
	EXPECT_EQ(0, frag_sites.at(350).count_forward_) << error_msg;
	EXPECT_EQ(0, frag_sites.at(350).count_reverse_) << error_msg;
	EXPECT_EQ(1.0, frag_sites.at(350).bias_) << error_msg;

	ref_.GetFragmentSites(frag_sites, 0, 10, 50+50, 193+50); // +50 from kMinNToSplitContigs
	EXPECT_EQ(143, frag_sites.size() );
	EXPECT_EQ(0, frag_sites.at(0).gc_);
	EXPECT_EQ(756483, frag_sites.at(0).start_surrounding_.sur_.at(0));
	EXPECT_EQ(787452, frag_sites.at(0).start_surrounding_.sur_.at(1));
	EXPECT_EQ(1017802, frag_sites.at(0).start_surrounding_.sur_.at(2));
	EXPECT_EQ(377552, frag_sites.at(0).end_surrounding_.sur_.at(0));
	EXPECT_EQ(787452, frag_sites.at(0).end_surrounding_.sur_.at(1));
	EXPECT_EQ(258513, frag_sites.at(0).end_surrounding_.sur_.at(2));
	EXPECT_EQ(0, frag_sites.at(0).count_forward_);
	EXPECT_EQ(0, frag_sites.at(0).count_reverse_);
	EXPECT_EQ(1.0, frag_sites.at(0).bias_);
	EXPECT_EQ(0, frag_sites.at(142).gc_); // N's are blocking this one (so we just have a placeholder here)
	EXPECT_EQ(-1, frag_sites.at(142).start_surrounding_.sur_.at(0));
	EXPECT_EQ(-1, frag_sites.at(142).start_surrounding_.sur_.at(1));
	EXPECT_EQ(-1, frag_sites.at(142).start_surrounding_.sur_.at(2));
	EXPECT_EQ(-1, frag_sites.at(142).end_surrounding_.sur_.at(0));
	EXPECT_EQ(-1, frag_sites.at(142).end_surrounding_.sur_.at(1));
	EXPECT_EQ(-1, frag_sites.at(142).end_surrounding_.sur_.at(2));
	EXPECT_EQ(0, frag_sites.at(142).count_forward_);
	EXPECT_EQ(0, frag_sites.at(142).count_reverse_);
	EXPECT_EQ(0.0, frag_sites.at(142).bias_);
}

void ReferenceTest::TestReplaceN(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	for(uintSeqLen i=50; i<150; ++i){
		at(at(ref_.reference_sequences_, 0), i) = 'N';
		at(at(ref_.reference_sequences_, 1), i+300) = 'N';
	}
	for(uintSeqLen i=54; i<104; ++i){
		at(at(ref_.reference_sequences_, 1), i) = 'N';
	}
	at(at(ref_.reference_sequences_, 0), 253) = 'N';
	at(at(ref_.reference_sequences_, 0), 256) = 'N';
	at(at(ref_.reference_sequences_, 0), 263) = 'N';

	uintSeqLen n_in_reference(253);
	array< uintSeqLen, 4 > old_values = {{0,0,0,0}};
	for( const auto &seq : ref_.reference_sequences_ ){
		for( auto base : seq ){
			if( !IsN(base) ){
				old_values.at( static_cast<uintBaseCall>(base) ) += 1;
			}
		}
	}

	ref_.ReplaceN(317);

	array< uintSeqLen, 4 > new_values = {{0,0,0,0}};
	for( const auto &seq : ref_.reference_sequences_ ){
		for( auto base : seq ){
			if( !IsN(base) ){
				new_values.at( static_cast<uintBaseCall>(base) ) += 1;
			}
		}
	}

	EXPECT_EQ(old_values.at(0)+old_values.at(1)+old_values.at(2)+old_values.at(3)+n_in_reference, new_values.at(0)+new_values.at(1)+new_values.at(2)+new_values.at(3)) << "Not all N's properly replaced\n";
	for(uintBaseCall i=4; i--; ){
		EXPECT_LT(old_values.at(i), new_values.at(i)) << "No N's converted to " << i << '\n';
	}

	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 50)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 51)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 52)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 53)));
	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 54)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 55)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 56)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 57)));
	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 146)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 147)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 148)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 0), 149)));

	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 350)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 351)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 352)));
	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 353)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 354)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 355)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 356)));
	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 357)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 446)));
	EXPECT_EQ(0, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 447)));
	EXPECT_EQ(1, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 448)));
	EXPECT_EQ(2, static_cast<uintBaseCall>(at(at(ref_.reference_sequences_, 1), 449)));
}

void ReferenceTest::TestExclusionRegions(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa").c_str() ) );

	// Shorten second sequence so it is invalid
	resize(at(ref_.reference_sequences_, 1), 140);

	// Add repeats and N's to first sequence
	auto &ref_seq = at(ref_.reference_sequences_, 0);
	for(auto pos=0; pos<3; ++pos){
		at( ref_seq, pos ) = 'N'; // Add N's at the start that should be added to kMinDistToContigEnds
	}

	for(auto pos=14; pos<23; ++pos){
		at( ref_seq, pos ) = 'N'; // Add N's that are below kMinNToSplitContigs so should be ignored
	}

	for(auto pos=95; pos<104; ++pos){
		at( ref_seq, pos ) = 'N'; // Add N's that are below kMinNToSplitContigs so should be ignored
	}
	for(auto pos=104; pos<125; ++pos){
		at( ref_seq, pos ) = 'C'; // But add repeats at the end that result in a single position marked as exclusion region
	}

	for(auto pos=210; pos<220; ++pos){
		at( ref_seq, pos ) = 'N'; // Add N's that hit kMinNToSplitContigs so should be excluded
	}

	auto copy_pos = 349;
	for(auto pos=355; pos<425; ++pos){
		at( ref_seq, pos ) = Complement::Dna5(at( ref_seq, --copy_pos )); // Add inverse repeat that is not joined, but first part is joined with N exclusion region before
	}

	for(auto pos=length(ref_seq)-2; pos<length(ref_seq); ++pos){
		at( ref_seq, pos ) = 'N'; // Add N's at the end that should be added to kMinDistToContigEnds
	}

	// Get exclusion regions
	ref_.PrepareExclusionRegions();
	ref_.ObtainExclusionRegions(5, 50); // Check if setting a too high value for ref seq id is an issue

	// Test results
	EXPECT_EQ(1, ref_.excluded_regions_.at(1).size());
	auto excluded_regions( ref_.excluded_regions_.at(0) );
	EXPECT_EQ(3, excluded_regions.size()); //5
	for( auto id=0; id < excluded_regions.size(); ++id ){
		// Make sure we do not check a non existing region if number of regions is less than expected
		switch(id){
		case 0:
			EXPECT_EQ(0, excluded_regions.at(id).first );
			EXPECT_EQ(53, excluded_regions.at(id).second );
			break;
		//case 1:
		//	EXPECT_EQ(105, excluded_regions.at(id).first );
		//	EXPECT_EQ(106, excluded_regions.at(id).second );
		//	break;
		case 1:
			EXPECT_EQ(160, excluded_regions.at(id).first );
			EXPECT_EQ(270, excluded_regions.at(id).second ); //319
			break;
		//case 3:
		//	EXPECT_EQ(385, excluded_regions.at(id).first );
		//	EXPECT_EQ(395, excluded_regions.at(id).second );
		//	break;
		case 4:
			EXPECT_EQ(length(ref_seq)-52, excluded_regions.at(id).first );
			EXPECT_EQ(length(ref_seq), excluded_regions.at(id).second );
			break;
		}
	}


	// Test checking functions
	EXPECT_TRUE( ref_.ObtainedExclusionRegionsForSequence(0) );
	EXPECT_TRUE( ref_.ObtainedExclusionRegionsForSequence(1) );
	EXPECT_TRUE( ref_.ExclusionRegionsCompletelyObtained() );

	uintSeqLen last_region_id(0), last_ref_seq(0);
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, 140, 161) );
	EXPECT_FALSE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, 140, 160) );
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, length(ref_seq)-2, length(ref_seq)) );
	EXPECT_EQ(2, last_region_id); //4
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, 1, 5) );
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, 52, 70) );
	EXPECT_EQ(0, last_region_id);
	EXPECT_FALSE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, 53, 70) );
	last_region_id = 100;
	last_ref_seq = 1;
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 0, length(ref_seq)-2, length(ref_seq)) );
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 1, 10, 20) );
	EXPECT_EQ(1, last_ref_seq);
	EXPECT_EQ(0, last_region_id);
	EXPECT_TRUE( ref_.FragmentExcluded(last_region_id, last_ref_seq, 1, 60, 70) );

	EXPECT_FALSE( ref_.ReferenceSequenceExcluded(0) );
	EXPECT_TRUE( ref_.ReferenceSequenceExcluded(1) );
	EXPECT_EQ(215, ref_.SumExcludedBases(0)); //275
	EXPECT_EQ(140, ref_.SumExcludedBases(1));
	EXPECT_EQ(3, ref_.NumExcludedRegions(0)); //5
	EXPECT_EQ(1, ref_.NumExcludedRegions(1));


	// Test clearing exclusion regions
	ref_.ClearExclusionRegions(2);
	EXPECT_EQ(0, ref_.excluded_regions_.at(0).size());
	EXPECT_EQ(0, ref_.excluded_regions_.at(1).size());

	ref_.ClearAllExclusionRegions();
	EXPECT_EQ(0, ref_.excluded_regions_.size());
}

void ReferenceTest::TestMethylationLoading(){
	string test_dir;
	ASSERT_TRUE( GetTestDir(test_dir) );
	ASSERT_TRUE( ref_.ReadFasta( (test_dir+"drosophila-GCF_000001215.4_cut.fna").c_str() ) );
	ref_.num_alleles_ = 2;

	EXPECT_FALSE( ref_.MethylationLoaded() );

	ASSERT_TRUE( ref_.PrepareMethylationFile(test_dir+"drosophila-methylation.bed") );
	ASSERT_TRUE( ref_.ReadMethylation(2) ); // Not actually reading anything, because second sequence is the first with entries
	EXPECT_TRUE( ref_.MethylationLoaded() );
	EXPECT_TRUE( ref_.MethylationLoadedForSequence(1) );
	EXPECT_FALSE( ref_.MethylationLoadedForSequence(2) );
	EXPECT_FALSE( ref_.MethylationCompletelyLoaded() );

	ASSERT_TRUE( ref_.ReadMethylation(ref_.NumberSequences()+2) );
	EXPECT_TRUE( ref_.MethylationLoadedForSequence(2) );
	EXPECT_TRUE( ref_.MethylationLoadedForSequence(ref_.NumberSequences()-1) );
	EXPECT_TRUE( ref_.MethylationCompletelyLoaded() );

	EXPECT_EQ(0, ref_.UnmethylatedRegions(1).size() );
	EXPECT_EQ(0, ref_.Unmethylation(1, 0).size() );
	EXPECT_EQ(0, ref_.Unmethylation(1, 1).size() );

	EXPECT_EQ(3, ref_.UnmethylatedRegions(2).size() );
	EXPECT_EQ(0, ref_.UnmethylatedRegions(2).at(0).first );
	EXPECT_EQ(100, ref_.UnmethylatedRegions(2).at(0).second );
	EXPECT_EQ(100, ref_.UnmethylatedRegions(2).at(1).first );
	EXPECT_EQ(101, ref_.UnmethylatedRegions(2).at(1).second );
	EXPECT_EQ(34520, ref_.UnmethylatedRegions(2).at(2).first );
	EXPECT_EQ(34521, ref_.UnmethylatedRegions(2).at(2).second );
	EXPECT_EQ(3, ref_.Unmethylation(2, 0).size() );
	EXPECT_DOUBLE_EQ(0.75, ref_.Unmethylation(2, 0).at(0) );
	EXPECT_DOUBLE_EQ(0.7, ref_.Unmethylation(2, 0).at(1) );
	EXPECT_DOUBLE_EQ(0.6, ref_.Unmethylation(2, 0).at(2) );
	EXPECT_EQ(3, ref_.Unmethylation(2, 1).size() );
	EXPECT_DOUBLE_EQ(0.74, ref_.Unmethylation(2, 1).at(0) );
	EXPECT_DOUBLE_EQ(0.69, ref_.Unmethylation(2, 1).at(1) );
	EXPECT_DOUBLE_EQ(0.59, ref_.Unmethylation(2, 1).at(2) );

	EXPECT_EQ(1, ref_.UnmethylatedRegions(9).size() );
	EXPECT_EQ(5000, ref_.UnmethylatedRegions(9).at(0).first );
	EXPECT_EQ(35000, ref_.UnmethylatedRegions(9).at(0).second );
	EXPECT_EQ(1, ref_.Unmethylation(9, 0).size() );
	EXPECT_DOUBLE_EQ(0.9, ref_.Unmethylation(9, 0).at(0) );
	EXPECT_EQ(1, ref_.Unmethylation(9, 1).size() );
	EXPECT_DOUBLE_EQ(0.9, ref_.Unmethylation(9, 1).at(0) );

	EXPECT_EQ(0, ref_.UnmethylatedRegions(10).size() );
	EXPECT_EQ(0, ref_.Unmethylation(10, 0).size() );
	EXPECT_EQ(0, ref_.Unmethylation(10, 1).size() );

	EXPECT_EQ(0, ref_.UnmethylatedRegions(20).size() );
	EXPECT_EQ(0, ref_.Unmethylation(20, 0).size() );
	EXPECT_EQ(0, ref_.Unmethylation(20, 1).size() );

	EXPECT_TRUE( ref_.ReadMethylation(ref_.NumberSequences()+2) ); // Check whether calling the read after reading everything does not crash
	ref_.ClearMethylation(2);
	ref_.ClearMethylation(ref_.NumberSequences());
	EXPECT_TRUE( ref_.ReadMethylation(ref_.NumberSequences()+2) ); // Check whether reading after clearing works
	ref_.ClearAllMethylation();
}

namespace reseq{
	TEST_F(ReferenceTest, Variants){
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		ASSERT_TRUE( ref_.ReadFasta( (test_dir+"ecoli-GCF_000005845.2_ASM584v2_genomic.fa").c_str() ) );

		TestVariantClass();
		TestInsertVariant();
		TestVariationPositionLoading();
		TestVariationLoading();
	}

	TEST_F(ReferenceTest, Functionality){
		TestLoadingAndAccess();
		TestGC();
		TestSumBias();
		TestGetFragmentSites();
		TestReplaceN();
	}

	TEST_F(ReferenceTest, GZip){
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa.gz").c_str() ) );

		ASSERT_EQ(2, ref_.NumberSequences()) << "The gz test file did not load properly\n";
	}

	TEST_F(ReferenceTest, BZip2){
		string test_dir;
		ASSERT_TRUE( GetTestDir(test_dir) );
		ASSERT_TRUE( ref_.ReadFasta( (test_dir+"reference-test.fa.bz2").c_str() ) );

		ASSERT_EQ(2, ref_.NumberSequences()) << "The bz2 test file did not load properly\n";
	}

	TEST_F(ReferenceTest, ExclusionRegions){
		TestExclusionRegions();
	}

	TEST_F(ReferenceTest, Methylation){
		TestMethylationLoading();
	}
}
