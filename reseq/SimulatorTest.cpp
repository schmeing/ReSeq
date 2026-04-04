#include "SimulatorTest.h"
using reseq::SimulatorTest;

#include <algorithm>
using std::max;
// include <array>
using std::array;
#include <atomic>
using std::atomic;
#include <bitset>
using std::bitset;
#include <chrono>
#include <condition_variable>
using std::condition_variable;
#include <mutex>
using std::mutex;
using std::unique_lock;
#include <string>
using std::string;
#include <thread>
using std::thread;
// include <vector>
using std::vector;

// include <seqan/seq_io.h>
using seqan::DnaString;
using seqan::length;

// include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::Percent;
using reseq::utilities::ReverseComplementorDna;

void SimulatorTest::Register() {
    // Guarantees that library is included
}

void SimulatorTest::CreateTestObject() {
    test_ = std::make_unique<Simulator>();
    ASSERT_TRUE(test_) << "Could not allocate memory for Simulator object\n";
}

void SimulatorTest::DeleteTestObject() {
    test_.reset();
}

void SimulatorTest::TearDown() {
    BasicTestClass::TearDown();
    DeleteTestObject();
}

void SimulatorTest::ChooseAlleles(vector<uintAlleleId>& chosen_allele_ids, vector<bool>& reverse_selection,
                                  uintAlleleId non_zero_strands, uintAlleleId possible_strands) const {
    EXPECT_FALSE(non_zero_strands <= possible_strands / 2)
        << "Not all alleles were select in the test. Something went wrong";

    // More than half the possible alleles need to be drawn: Draw inverse selection
    chosen_allele_ids.clear();
    reverse_selection.clear();
    reverse_selection.resize(possible_strands, true);
    while (chosen_allele_ids.size() < possible_strands - non_zero_strands) {
        double random_value = 0.5;
        test_->SelectAllele(chosen_allele_ids, reverse_selection, possible_strands, random_value);
    }

    EXPECT_EQ(possible_strands - non_zero_strands, chosen_allele_ids.size());

    test_->ReverseSelection(chosen_allele_ids, reverse_selection, possible_strands);

    EXPECT_EQ(possible_strands, chosen_allele_ids.size());
}

void SimulatorTest::TestCoverageConversion() {
    // CoveragePropLostFromAdapters
    DataStats stats(nullptr);

    stats.read_lengths_by_fragment_length_.at(1)[200][150] = 10;
    stats.read_lengths_by_fragment_length_.at(0)[150][150] = 10;
    stats.read_lengths_by_fragment_length_.at(1)[100][150] = 5;           // 5*50
    stats.read_lengths_by_fragment_length_.at(1)[50][150] = 10;           // (10-5)*100 [-5 due to next line]
    stats.non_mapped_read_lengths_by_fragment_length_.at(1)[50][150] = 5; // 5*150
    stats.read_lengths_by_fragment_length_.at(0)[0][150] = 10;            // 10*150

    double adapter_part = test_->CoveragePropLostFromAdapters(stats);
    EXPECT_DOUBLE_EQ(3000.0 / (45 * 150), adapter_part);

    // CoverageToNumberPairs
    double coverage = 100;
    uintRefLenCalc total_ref_size = 50000;
    double average_read_length = 150;

    uintFragCount total_pairs =
        test_->CoverageToNumberPairs(coverage, total_ref_size, average_read_length, adapter_part);
    EXPECT_EQ(30000, total_pairs);

    // NumberPairsToCoverage
    EXPECT_NEAR(coverage, test_->NumberPairsToCoverage(total_pairs, total_ref_size, average_read_length, adapter_part),
                0.01);
}

void SimulatorTest::TestSelectAllele() {
    // Test 2 alleles
    vector<uintAlleleId> chosen_allele_ids;
    vector<bool> reverse_selection;

    reverse_selection.resize(2, true);
    while (chosen_allele_ids.size() < 2) {
        test_->SelectAllele(chosen_allele_ids, reverse_selection, 2, 0.5);
    }

    EXPECT_EQ(1, chosen_allele_ids.at(0));
    EXPECT_EQ(0, chosen_allele_ids.at(1));

    // Test 4 alleles
    chosen_allele_ids.clear();
    reverse_selection.clear();
    reverse_selection.resize(4, true);
    while (chosen_allele_ids.size() < 4) {
        test_->SelectAllele(chosen_allele_ids, reverse_selection, 4, 0.5);
    }

    EXPECT_EQ(2, chosen_allele_ids.at(0));
    EXPECT_EQ(1, chosen_allele_ids.at(1));
    EXPECT_EQ(3, chosen_allele_ids.at(2));
    EXPECT_EQ(0, chosen_allele_ids.at(3));
}

void SimulatorTest::TestVariationInInnerLoopOfSimulateFromGivenBlock(
    Simulator::VariantBiasVarModifiers& bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position,
    uintSeqLen frag_length_from, uintSeqLen frag_length_to, uintAlleleId num_valid_alleles,
    array<vector<intVariantId>, 2> unhandled_variant_id, array<vector<uintSeqLen>, 2> unhandled_bases_in_variant,
    array<vector<intSeqShift>, 2> gc_mod, array<vector<intSeqShift>, 2> end_pos_shift,
    array<uintSeqLen, 2> modified_start_pos, array<Reference*, 2> comp_ref) {
    // Preparation
    DataStats stats(nullptr);
    stats.read_lengths_.at(0)[100];
    stats.read_lengths_.at(1)[100];
    stats.errors_.PrepareSimulation();
    Simulator::SimPair sim_reads;
    Surrounding surrounding_end, comp_surrounding;

    auto num_tests = 0;

    uintPercent gc_perc;
    uintSeqLen cur_end_position;

    vector<uintAlleleId> possible_alleles, chosen_allele_ids;
    vector<bool> reverse_selection;

    double probability_chosen = 1.0;

    // Inner loop
    auto frag_len_start = frag_length_from;
    test_->GetPossibleAlleles(possible_alleles, species_reference_, bias_mod, cur_start_position, ref_seq_id);

    for (auto fragment_length = frag_len_start; fragment_length < frag_length_to; ++fragment_length) {
        if (test_->ProbabilityAboveThreshold(probability_chosen, ref_seq_id, fragment_length)) {
            auto non_zero_strands = stats.FragmentDistribution().DrawNumberNonZeroStrands(
                possible_alleles.size(), test_->NonZeroThreshold(ref_seq_id, fragment_length), probability_chosen);
            EXPECT_EQ(2 * possible_alleles.size(), non_zero_strands);

            if (non_zero_strands) {
                ChooseAlleles(chosen_allele_ids, reverse_selection, non_zero_strands, 2 * possible_alleles.size());

                for (auto chosen_id : chosen_allele_ids) {
                    auto allele = possible_alleles.at(chosen_id / 2);
                    bool strand = chosen_id % 2;

                    test_->PrepareBiasModForCurrentFragmentLength(bias_mod, ref_seq_id, species_reference_,
                                                                  cur_start_position, fragment_length, allele);

                    comp_ref.at(allele)->ReverseSurrounding(comp_surrounding, ref_seq_id,
                                                            modified_start_pos.at(allele) + fragment_length - 1);
                    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_end_.at(allele).sur_.at(0))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Allele: " << allele << std::endl;
                    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_end_.at(allele).sur_.at(1))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Allele: " << allele << std::endl;
                    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_end_.at(allele).sur_.at(2))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Allele: " << allele << std::endl;

                    auto result = fragment_length - frag_length_from;
                    EXPECT_EQ(unhandled_variant_id.at(allele).at(result), bias_mod.unhandled_variant_id_.at(allele))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
                    EXPECT_EQ(unhandled_bases_in_variant.at(allele).at(result),
                              bias_mod.unhandled_bases_in_variant_.at(allele))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
                    EXPECT_EQ(gc_mod.at(allele).at(result), bias_mod.gc_mod_.at(allele))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;
                    EXPECT_EQ(end_pos_shift.at(allele).at(result), bias_mod.end_pos_shift_.at(allele))
                        << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                        << " Variant position: " << bias_mod.start_variant_pos_ << " Allele: " << allele << std::endl;

                    cur_end_position = cur_start_position + fragment_length + bias_mod.end_pos_shift_.at(allele);
                    if (cur_end_position <= species_reference_.SequenceLength(ref_seq_id)) {
                        // Determine how many read pairs are generated for this strand and allele at this position with
                        // this fragment_length
                        gc_perc = test_->GetGCPercent(bias_mod, ref_seq_id, species_reference_, cur_end_position,
                                                      fragment_length, allele);
                        EXPECT_EQ(Percent(comp_ref.at(allele)->GCContentAbsolut(
                                              ref_seq_id, modified_start_pos.at(allele),
                                              modified_start_pos.at(allele) + fragment_length),
                                          fragment_length),
                                  gc_perc);

                        test_->GetOrgSeq(sim_reads, strand, allele, fragment_length, cur_start_position,
                                         cur_end_position, ref_seq_id, species_reference_, stats, bias_mod);

                        EXPECT_TRUE(infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id),
                                          modified_start_pos.at(allele),
                                          modified_start_pos.at(allele) + fragment_length) ==
                                    prefix(sim_reads.at(strand).org_seq_, fragment_length))
                            << prefix(sim_reads.at(strand).org_seq_, fragment_length) << std::endl
                            << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                            << " Allele: " << allele << std::endl;
                        EXPECT_TRUE(ReverseComplementorDna(infix(comp_ref.at(allele)->ReferenceSequence(ref_seq_id),
                                                                 modified_start_pos.at(allele),
                                                                 modified_start_pos.at(allele) + fragment_length)) ==
                                    prefix(sim_reads.at(!strand).org_seq_, fragment_length))
                            << prefix(sim_reads.at(!allele).org_seq_, fragment_length) << std::endl
                            << "Start position: " << cur_start_position << " Fragment length: " << fragment_length
                            << " Allele: " << allele << std::endl;

                        ++num_tests;
                    }
                }
            }
        }
    }

    EXPECT_EQ(num_valid_alleles * 2 * max(unhandled_variant_id.at(0).size(), unhandled_variant_id.at(1).size()),
              num_tests);
}

void SimulatorTest::TestVariationInSimulateFromGivenBlock() {
    species_reference_.num_alleles_ = 2;

    // samtools faidx ../test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa NC_000913.3:1001-1020
    // GTTGCGAGATTTGGACGGAC
    species_reference_.variants_.clear();
    species_reference_.variants_.resize(1);
    species_reference_.variants_.at(0).clear();
    std::array<uintAlleleBitArray, reseq::Reference::Variant::kMaxAlleles / 64> present_in_alleles = {
        2}; // bitwise [0:no, 1:yes]
    species_reference_.variants_.at(0).emplace_back(455, "", present_in_alleles);
    species_reference_.variants_.at(0).emplace_back(1002, "TAC", present_in_alleles);
    species_reference_.variants_.at(0).emplace_back(1004, "TGA", present_in_alleles);
    species_reference_.variants_.at(0).emplace_back(1008, "", present_in_alleles);
    species_reference_.variants_.at(0).emplace_back(1011, "C", present_in_alleles);
    species_reference_.variants_.at(0).emplace_back(1012, "", present_in_alleles);
    // GTT--GC--GAGATTTGGACGGAC
    // GTTACGCGAGAG-TTC-GACGGAC

    Reference test_ref;
    resize(test_ref.reference_sequences_, 1);
    at(test_ref.reference_sequences_, 0) = infix(species_reference_.ReferenceSequence(0), 0, 455);
    at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 456, 1003);
    at(test_ref.reference_sequences_, 0) += "AC";
    at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1003, 1004);
    at(test_ref.reference_sequences_, 0) += "TGA";
    at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1005, 1008);
    at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1009, 1012);
    at(at(test_ref.reference_sequences_, 0), 1013) = 'C';
    at(test_ref.reference_sequences_, 0) += infix(species_reference_.ReferenceSequence(0), 1013, 2000);
    Surrounding comp_surrounding;

    test_->coverage_groups_.resize(species_reference_.NumberSequences(), 0);
    test_->non_zero_thresholds_.resize(1);
    test_->non_zero_thresholds_.at(0).resize(100, {0.75, 0.31640625});

    Simulator::VariantBiasVarModifiers bias_mod(1, 2);
    uintRefSeqId ref_seq_id = 0;
    Surrounding surrounding_start;

    // Position 1003
    uintSeqLen cur_start_position = 1003;
    bias_mod.first_variant_id_ = 2; // Insertion at position 2 has already been processed

    species_reference_.ForwardSurrounding(surrounding_start, ref_seq_id, cur_start_position);
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
    EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
    EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

    test_ref.ForwardSurrounding(comp_surrounding, ref_seq_id, 1004);
    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

    TestVariationInInnerLoopOfSimulateFromGivenBlock(
        bias_mod, ref_seq_id, cur_start_position, 1, 13, 2,
        {{{2, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6}, {2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 6, 6}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -1, -2, -2, -2, -2, -1, -1, -1, 0, 0}}},
        {cur_start_position, 1004}, {&species_reference_, &test_ref});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(2, bias_mod.first_variant_id_);
    EXPECT_EQ(0, bias_mod.start_variant_pos_);

    // Position 1004 + 0
    species_reference_.ForwardSurrounding(surrounding_start, ref_seq_id, ++cur_start_position);
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(0)) << " (Result)";
    EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(1)) << " (Result)";
    EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(0).sur_.at(2)) << " (Result)";

    test_ref.ForwardSurrounding(comp_surrounding, ref_seq_id, 1005);
    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

    TestVariationInInnerLoopOfSimulateFromGivenBlock(
        bias_mod, ref_seq_id, cur_start_position, 1, 12, 2,
        {{{3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6}, {2, 2, 3, 3, 3, 3, 4, 4, 5, 6, 6}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}}},
        {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, -2, -2, -2, -2, -1, -1, -1, 0, 0}}}, {cur_start_position, 1005},
        {&species_reference_, &test_ref});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(2, bias_mod.first_variant_id_);
    EXPECT_EQ(1, bias_mod.start_variant_pos_);

    // Position 1004 + 1
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    test_ref.ForwardSurrounding(comp_surrounding, ref_seq_id, 1006);
    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

    TestVariationInInnerLoopOfSimulateFromGivenBlock(
        bias_mod, ref_seq_id, cur_start_position, 1, 11, 1, {{{}, {3, 3, 3, 3, 3, 4, 4, 5, 6, 6}}},
        {{{}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}}, {{{}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}}},
        {{{}, {0, -1, -1, -1, -1, 0, 0, 0, 1, 1}}}, {0, 1006}, {nullptr, &test_ref});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(2, bias_mod.first_variant_id_);
    EXPECT_EQ(2, bias_mod.start_variant_pos_);

    // Position 1004 + 2
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    test_ref.ForwardSurrounding(comp_surrounding, ref_seq_id, 1007);
    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

    TestVariationInInnerLoopOfSimulateFromGivenBlock(
        bias_mod, ref_seq_id, cur_start_position, 1, 10, 1, {{{}, {3, 3, 3, 3, 4, 4, 5, 6, 6}}},
        {{{}, {0, 0, 0, 0, 0, 0, 0, 0, 0}}}, {{{}, {-1, -1, -1, -1, -1, -1, 0, -1, -1}}},
        {{{}, {0, 0, 0, 0, 1, 1, 1, 2, 2}}}, {0, 1007}, {nullptr, &test_ref});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(3, bias_mod.first_variant_id_);
    EXPECT_EQ(0, bias_mod.start_variant_pos_);

    // Position 1008
    cur_start_position = 1008;
    species_reference_.ForwardSurrounding(surrounding_start, ref_seq_id, cur_start_position);
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
    EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
    EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

    TestVariationInInnerLoopOfSimulateFromGivenBlock(bias_mod, ref_seq_id, cur_start_position, 1, 8, 1,
                                                     {{{4, 4, 4, 5, 6, 6, 6}, {}}}, {{{0, 0, 0, 0, 0, 0, 0}, {}}},
                                                     {{{0, 0, 0, 0, 0, 0, 0}, {}}}, {{{0, 0, 0, 0, 0, 0, 0}, {}}},
                                                     {cur_start_position, 0}, {&species_reference_, nullptr});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(4, bias_mod.first_variant_id_);
    EXPECT_EQ(0, bias_mod.start_variant_pos_);

    // Position 1011
    cur_start_position = 1011;
    species_reference_.ForwardSurrounding(surrounding_start, ref_seq_id, cur_start_position);
    test_->PrepareBiasModForCurrentStartPos(
        bias_mod, ref_seq_id, species_reference_, cur_start_position, 1,
        surrounding_start); // cur_end_position = cur_start_position -> Fragment length = 1
    EXPECT_EQ(surrounding_start.sur_.at(0), bias_mod.surrounding_start_.at(0).sur_.at(0));
    EXPECT_EQ(surrounding_start.sur_.at(1), bias_mod.surrounding_start_.at(0).sur_.at(1));
    EXPECT_EQ(surrounding_start.sur_.at(2), bias_mod.surrounding_start_.at(0).sur_.at(2));

    test_ref.ForwardSurrounding(comp_surrounding, ref_seq_id, 1013);
    EXPECT_EQ(comp_surrounding.sur_.at(0), bias_mod.surrounding_start_.at(1).sur_.at(0))
        << bitset<20>(surrounding_start.sur_.at(0)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(0)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(0)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(1), bias_mod.surrounding_start_.at(1).sur_.at(1))
        << bitset<20>(surrounding_start.sur_.at(1)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(1)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(1)) << " (Result)";
    EXPECT_EQ(comp_surrounding.sur_.at(2), bias_mod.surrounding_start_.at(1).sur_.at(2))
        << bitset<20>(surrounding_start.sur_.at(2)) << " (Original)" << std::endl
        << bitset<20>(comp_surrounding.sur_.at(2)) << " (Goal)" << std::endl
        << bitset<20>(bias_mod.surrounding_start_.at(1).sur_.at(2)) << " (Result)";

    TestVariationInInnerLoopOfSimulateFromGivenBlock(bias_mod, ref_seq_id, cur_start_position, 1, 5, 2,
                                                     {{{5, 6, 6, 6}, {5, 6, 6, 6}}}, {{{0, 0, 0, 0}, {0, 0, 0, 0}}},
                                                     {{{0, 0, 0, 0}, {1, 0, 0, 0}}}, {{{0, 0, 0, 0}, {0, 1, 1, 1}}},
                                                     {cur_start_position, 1013}, {&species_reference_, &test_ref});

    test_->CheckForInsertedBasesToStartFrom(bias_mod, 0, cur_start_position, species_reference_);
    EXPECT_EQ(5, bias_mod.first_variant_id_);
    EXPECT_EQ(0, bias_mod.start_variant_pos_);
}

void SimulatorTest::TestWrittenBlocksSynchronization() {
    // Test that written_blocks_ is properly incremented so that threads
    // waiting to write later blocks are unblocked. This is the core fix
    // for the 10k read hang (GitHub issue #24).
    //
    // Without the ++written_blocks_ after successful flush, the condition
    // variable wait in WriteSingleReads would never be satisfied for
    // block 1+, causing a deadlock after the first block (10,000 reads).

    const uintFragCount num_blocks = 5; // Enough to test multi-block ordering
    atomic<uintFragCount> completed_blocks(0);
    atomic<bool> deadlock_detected(false);

    // Reset state
    test_->written_blocks_ = 0;
    test_->simulation_error_ = false;
    test_->written_records_ = 0;

    // Open a temp output file so FlushWriteValues has somewhere to write
    string tmp_file = "/tmp/reseq_sync_test.fq";
    seqan::open(test_->dest_.at(0), tmp_file.c_str());

    // Launch threads that call WriteSingleReads with increasing block numbers.
    // Each thread writes a small batch. If written_blocks_ is not incremented,
    // threads for block >= 1 will wait forever (deadlock).
    vector<thread> threads;
    for (uintFragCount block = 0; block < num_blocks; ++block) {
        threads.emplace_back([this, block, &completed_blocks, &deadlock_detected]() {
            seqan::StringSet<seqan::CharString> ids;
            seqan::StringSet<seqan::Dna5String> seqs;
            seqan::StringSet<seqan::CharString> quals;

            // Add one read per block
            seqan::appendValue(ids, "test_read");
            seqan::Dna5String seq = "ACGT";
            seqan::appendValue(seqs, seq);
            seqan::appendValue(quals, "IIII");

            bool success = test_->WriteSingleReads(block, ids, seqs, quals);
            if (success) {
                ++completed_blocks;
            }
        });
    }

    // Wait with a timeout to detect deadlocks
    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    for (auto& t : threads) {
        if (t.joinable()) {
            // We can't do a timed join in standard C++, so we join and rely
            // on the test framework's timeout. But we track completion count.
            t.join();
        }
    }

    seqan::close(test_->dest_.at(0));
    std::remove(tmp_file.c_str());

    // All blocks must have completed
    EXPECT_EQ(num_blocks, completed_blocks)
        << "Not all blocks completed - written_blocks_ synchronization is broken (10k hang bug)";
    EXPECT_EQ(num_blocks, test_->written_blocks_) << "written_blocks_ counter does not match expected count";
}

void SimulatorTest::TestBlockLifecycle() {
    // Verify that AllocBlock assigns correct indices and that FreeBlock + re-alloc reuses the slot.

    // Initially no blocks exist
    EXPECT_EQ(0u, test_->blocks_.size());
    EXPECT_TRUE(test_->free_block_indices_.empty());

    // Allocate first block
    size_t idx0 = test_->AllocBlock(0, 100, SIZE_MAX, 42);
    EXPECT_EQ(0u, idx0);
    EXPECT_EQ(1u, test_->blocks_.size());
    ASSERT_NE(nullptr, test_->blocks_[idx0].get());
    EXPECT_EQ(idx0, test_->blocks_[idx0]->block_idx_);

    // Allocate second block
    size_t idx1 = test_->AllocBlock(0, 200, SIZE_MAX, 43);
    EXPECT_EQ(1u, idx1);
    EXPECT_EQ(2u, test_->blocks_.size());
    ASSERT_NE(nullptr, test_->blocks_[idx1].get());
    EXPECT_EQ(idx1, test_->blocks_[idx1]->block_idx_);

    // Free the first block — its slot should enter the free list
    test_->FreeBlock(idx0);
    EXPECT_EQ(nullptr, test_->blocks_[idx0].get());
    ASSERT_EQ(1u, test_->free_block_indices_.size());
    EXPECT_EQ(idx0, test_->free_block_indices_.back());

    // Next allocation must reuse the freed index
    size_t idx_reused = test_->AllocBlock(0, 300, SIZE_MAX, 44);
    EXPECT_EQ(idx0, idx_reused) << "Freed slot should be reused";
    EXPECT_TRUE(test_->free_block_indices_.empty()) << "Free list should be empty after reuse";
    EXPECT_EQ(2u, test_->blocks_.size()) << "Deque size must not grow when reusing a slot";
    ASSERT_NE(nullptr, test_->blocks_[idx_reused].get());
    EXPECT_EQ(idx_reused, test_->blocks_[idx_reused]->block_idx_);
}

void SimulatorTest::TestPartnerLinkage() {
    // Allocate a forward block (no partner yet) and then a reverse block pointing to the forward block.
    size_t fwd_idx = test_->AllocBlock(0, 1000, SIZE_MAX, 10);
    size_t rev_idx = test_->AllocBlock(0, 1000, fwd_idx, 11);

    ASSERT_NE(nullptr, test_->blocks_[fwd_idx].get());
    ASSERT_NE(nullptr, test_->blocks_[rev_idx].get());

    // Forward block has no partner (was created with SIZE_MAX)
    EXPECT_EQ(SIZE_MAX, test_->blocks_[fwd_idx]->partner_block_idx_);

    // Reverse block's partner points back to the forward block index
    EXPECT_EQ(fwd_idx, test_->blocks_[rev_idx]->partner_block_idx_);

    // They must be at distinct indices
    EXPECT_NE(fwd_idx, rev_idx);
}

void SimulatorTest::TestCleanupFreesList() {
    // Allocate several blocks, free a subset, verify the free list captures exactly those indices.

    size_t idx0 = test_->AllocBlock(0, 0, SIZE_MAX, 1);
    size_t idx1 = test_->AllocBlock(0, 10, SIZE_MAX, 2);
    size_t idx2 = test_->AllocBlock(0, 20, SIZE_MAX, 3);
    size_t idx3 = test_->AllocBlock(0, 30, SIZE_MAX, 4);

    EXPECT_EQ(4u, test_->blocks_.size());
    EXPECT_TRUE(test_->free_block_indices_.empty());

    // Free idx1 and idx3 (every other block)
    test_->FreeBlock(idx1);
    test_->FreeBlock(idx3);

    EXPECT_EQ(2u, test_->free_block_indices_.size());

    // Both freed indices must appear in the free list (order: LIFO)
    EXPECT_EQ(idx3, test_->free_block_indices_.at(1));
    EXPECT_EQ(idx1, test_->free_block_indices_.at(0));

    // Freed slots must be null; live slots must still be valid
    EXPECT_EQ(nullptr, test_->blocks_[idx1].get());
    EXPECT_EQ(nullptr, test_->blocks_[idx3].get());
    ASSERT_NE(nullptr, test_->blocks_[idx0].get());
    ASSERT_NE(nullptr, test_->blocks_[idx2].get());
}

namespace reseq {
TEST_F(SimulatorTest, BasicFunctonality) {
    CreateTestObject();

    TestCoverageConversion();
    TestSelectAllele();
}

TEST_F(SimulatorTest, WrittenBlocksSynchronization) {
    CreateTestObject();
    TestWrittenBlocksSynchronization();
}

TEST_F(SimulatorTest, Variants) {
    CreateTestObject();
    string test_dir;
    ASSERT_TRUE(GetTestDir(test_dir));
    LoadReference(test_dir + "ecoli-GCF_000005845.2_ASM584v2_genomic.fa");

    TestVariationInSimulateFromGivenBlock();
}

TEST_F(SimulatorTest, BlockLifecycle) {
    CreateTestObject();
    TestBlockLifecycle();
}

TEST_F(SimulatorTest, PartnerLinkage) {
    CreateTestObject();
    TestPartnerLinkage();
}

TEST_F(SimulatorTest, CleanupFreesList) {
    CreateTestObject();
    TestCleanupFreesList();
}
} // namespace reseq
