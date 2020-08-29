#ifndef SIMULATORTEST_H
#define SIMULATORTEST_H
#include "Simulator.h"

#include <array>
#include <stdint.h>
#include <vector>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class SimulatorTest : public BasicTestClassWithReference{
	public: 
		static void Register();

	protected:
		Simulator *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

		void ChooseAlleles(std::vector<uintAlleleId> &chosen_allele_ids, std::vector<bool> &reverse_selection, uintAlleleId non_zero_strands, uintAlleleId possible_strands) const;

		void TestCoverageConversion();
		void TestSelectAllele();

		void TestVariationInInnerLoopOfSimulateFromGivenBlock( Simulator::VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position, uintSeqLen frag_length_from, uintSeqLen frag_length_to, uintAlleleId num_valid_alleles, std::array<std::vector<intVariantId>, 2> unhandled_variant_id, std::array<std::vector<uintSeqLen>, 2> unhandled_bases_in_variant, std::array<std::vector<intSeqShift>, 2> gc_mod, std::array<std::vector<intSeqShift>, 2> end_pos_shift, std::array<uintSeqLen, 2> modified_start_pos, std::array<Reference *, 2> comp_ref );
		void TestVariationInSimulateFromGivenBlock();



	public:
		SimulatorTest():
			test_(NULL)
			{
		}
	};
}

#endif // SIMULATORTEST_H
