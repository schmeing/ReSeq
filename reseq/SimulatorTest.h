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
	protected:
		Simulator *test_;

		void CreateTestObject();
		void DeleteTestObject();

		virtual void TearDown();

		void TestCoverageConversion();
		void TestSurroundingModifiers();

		void TestVariationInInnerLoopOfSimulateFromGivenBlock(Simulator::VariantBiasVarModifiers &bias_mod, uint32_t ref_seq_id, uint32_t cur_start_position, uint32_t frag_length_from, uint32_t frag_length_to, const std::vector<std::array<uint32_t, 2>> unhandled_variant_id, const std::vector<std::array<uint32_t, 2>> unhandled_bases_in_variant, const std::vector<std::array<int32_t, 2>> gc_mod, const std::vector<std::array<int32_t, 2>> end_pos_shift, std::array<uint32_t, 2> modified_start_pos, std::array<const Reference *, 2> comp_ref);
		void TestVariationInInnerLoopOfSimulateFromGivenBlock(Simulator::VariantBiasVarModifiers &bias_mod, uint32_t ref_seq_id, uint32_t cur_start_position, uint32_t frag_length_from, uint32_t frag_length_to, uint16_t allele, const std::vector<uint32_t> unhandled_variant_id, const std::vector<uint32_t> unhandled_bases_in_variant, const std::vector<int32_t> gc_mod, const std::vector<int32_t> end_pos_shift, uint32_t modified_start_pos, const Reference &comp_ref);
		void TestVariationInSimulateFromGivenBlock();

	public:
		SimulatorTest():
			test_(NULL)
			{
		}
	};
}

#endif // SIMULATORTEST_H
