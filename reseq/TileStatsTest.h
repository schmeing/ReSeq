#ifndef TILESTATSTEST_H
#define TILESTATSTEST_H
#include "TileStats.h"

#include <stdint.h>
#include <string>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"

namespace reseq{
	class TileStatsTest : public BasicTestClass{
	protected:
		TileStats *test_;

		void CreateTestObject();
		void DeleteTestObject();

		void TestProperTileFormat(const char *read_id_string, const std::string &format, uint16_t expected_id, uint16_t expected_tile, uint16_t expected_colon_number);
		void TestProperTileFormatWithGivenColonNumber(const char * read_id_string, const std::string &format, uint16_t expected_id);
		void TestErroneousTileFormat(const char * read_id_string, bool reset_accession_info, uint16_t expected_id, const std::string& comment, const std::string& expected_error);

		virtual void TearDown();

	public:
		TileStatsTest():
			test_(NULL)
			{
		}

		static void TestSrr490124Equality(const TileStats &test, const char *context, bool test_tile_information=true);
		static void TestTiles(const TileStats &test);
	};
}

#endif // TILESTATSTEST_H
