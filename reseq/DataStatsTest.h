#ifndef DATASTATSTEST_H
#define DATASTATSTEST_H
#include "DataStats.h"

#include <stdint.h>
#include <string>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class DataStatsTest : public BasicTestClassWithReference{
	public:
		static void Register();

	protected:
		DataStats *test_;

		void CreateTestObject(Reference *ref);
		void DeleteTestObject();
		void LoadStats(const std::string &stats_file, bool ignore_tiles=false, bool calculate_biases=false, const std::string adapter="TruSeq_single", const std::string variants="");

		void TestSequenceContent(uintTempSeq template_segment, uintReadLen at_pos, uintNucCount cont_a, uintNucCount cont_c, uintNucCount cont_g, uintNucCount cont_t, uintNucCount cont_n, const char * context );
		void TestSequenceContentReference(uintTempSeq template_segment, uintTempSeq strand, uintReadLen at_pos, uintNucCount cont_a, uintNucCount cont_c, uintNucCount cont_g, uintNucCount cont_t, const char * context );

		void TestSrr490124Equality(const char *context, bool test_tile_information=true, bool bwa=false);
		void TestTiles();
		void TestDuplicates();
		void TestVariants();

		void TestCrossDuplicates();
		void TestCoverage();

		void TestAdapters(const char *context, bool bwa=false);
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
