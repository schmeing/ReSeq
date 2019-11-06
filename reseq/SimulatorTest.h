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

		void TestCoverageConversion();

		void TestVariationInInnerLoopOfSimulateFromGivenBlock(Simulator::VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position, uintSeqLen frag_length_from, uintSeqLen frag_length_to, const std::vector<std::array<intVariantId, 2>> unhandled_variant_id, const std::vector<std::array<uintSeqLen, 2>> unhandled_bases_in_variant, const std::vector<std::array<intSeqShift, 2>> gc_mod, const std::vector<std::array<intSeqShift, 2>> end_pos_shift, std::array<uintSeqLen, 2> modified_start_pos, std::array<const Reference *, 2> comp_ref);
		void TestVariationInInnerLoopOfSimulateFromGivenBlock(Simulator::VariantBiasVarModifiers &bias_mod, uintRefSeqId ref_seq_id, uintSeqLen cur_start_position, uintSeqLen frag_length_from, uintSeqLen frag_length_to, uintAlleleId allele, const std::vector<intVariantId> unhandled_variant_id, const std::vector<uintSeqLen> unhandled_bases_in_variant, const std::vector<intSeqShift> gc_mod, const std::vector<intSeqShift> end_pos_shift, uintSeqLen modified_start_pos, const Reference &comp_ref);
		void TestVariationInSimulateFromGivenBlock();

	public:
		SimulatorTest():
			test_(NULL)
			{
		}
	};
}

#endif // SIMULATORTEST_H
