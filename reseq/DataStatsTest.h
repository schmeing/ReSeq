#ifndef DATASTATSTEST_H
#define DATASTATSTEST_H
#include "DataStats.h"

#include <stdint.h>
#include <string>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class DataStatsTest : public BasicTestClassWithReference{
	protected:
		DataStats *test_;

		void CreateTestObject(Reference *ref);
		void DeleteTestObject();
		void LoadStats(const std::string &stats_file, bool ignore_tiles=false, bool calculate_biases=false, const std::string adapter="TruSeq_v2", const std::string variants="");

		void TestSequenceContent(uint16_t template_segment, uint32_t at_pos, uint64_t cont_a, uint64_t cont_c, uint64_t cont_g, uint64_t cont_t, uint64_t cont_n, const char * context );
		void TestSequenceContentReference(uint16_t template_segment, uint16_t strand, uint32_t at_pos, uint64_t cont_a, uint64_t cont_c, uint64_t cont_g, uint64_t cont_t, const char * context );

		void TestSrr490124Equality(const char *context, bool test_tile_information=true);
		void TestTiles();
		void TestDuplicates();
		void TestVariants();

		void TestCrossDuplicates();
		void TestCoverage();

		void TestAdapters();
		void TestNexteraAdapters();

		virtual void TearDown();

	public:
		DataStatsTest():
			test_(NULL)
			{
		}
	};
}

#endif // DATASTATSTEST_H
