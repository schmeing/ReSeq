#include "TileStatsTest.h"
using reseq::TileStatsTest;

#include <sstream>
using std::stringstream;
//include <string>
using std::string;
using std::to_string;

#include <seqan/bam_io.h>
using seqan::CharString;
using seqan::toCString;

void TileStatsTest::Register(){
	// Guarantees that library is included
}

void TileStatsTest::CreateTestObject(){
	ASSERT_TRUE( test_ = new TileStats() ) << "Could not allocate memory for DataStats object\n";
}

void TileStatsTest::DeleteTestObject(){
	if( test_ ){
		delete test_;
		test_ = NULL;
	}
}

void TileStatsTest::TestProperTileFormat(
		const char *read_id_string,
		const string &format,
		uint16_t expected_id,
		uint16_t expected_tile,
		uint16_t expected_colon_number){
	test_->ResetTileIdAccessionInfo();

	string error_msg = "Could not determine tile id for " + format + " format\n";

	stringstream returned_error_message("");
	uint16_t tile_id;

	EXPECT_TRUE( test_->EnterTile(read_id_string, &returned_error_message) ) << error_msg;
	EXPECT_EQ( "", returned_error_message.str() ) << error_msg;
	EXPECT_TRUE( test_->tile_accessible_ ) << error_msg;
	EXPECT_EQ( expected_colon_number, test_->tile_colon_number_ ) << "tile_colon_number_ wrongly set for " + format + " format\n";

	EXPECT_TRUE( test_->GetTileId(tile_id, read_id_string) ) << error_msg;
	EXPECT_EQ( expected_id, tile_id ) << error_msg;
	EXPECT_EQ( expected_tile, test_->tiles_[tile_id] ) << "Tile wrongly determined for " + format + " format\n";
}

void TileStatsTest::TestProperTileFormatWithGivenColonNumber(const char * read_id_string, const string &format, uint16_t expected_id){
	string error_msg = "Could not determine tile id for " + format + " format with given tile_colon_number_\n";
	stringstream returned_error_message("");
	uint16_t tile_id;
	EXPECT_TRUE( test_->EnterTile(read_id_string, &returned_error_message) ) << error_msg;
	EXPECT_EQ( "", returned_error_message.str() ) << error_msg;

	EXPECT_TRUE( test_->GetTileId(tile_id, read_id_string) ) << error_msg;
	EXPECT_EQ( expected_id, tile_id ) << error_msg;
}

void TileStatsTest::TestErroneousTileFormat(
		const char * read_id_string,
		bool reset_accession_info,
		uint16_t expected_id,
		const string& comment,
		const string& expected_error){
	if( reset_accession_info ){
		test_->ResetTileIdAccessionInfo();
	}

	string error_msg = "'GetTileId' did not set tile to default 0 in case of an error:'";
	error_msg += read_id_string;
	error_msg += "'\n";

	stringstream returned_error_message("");
	uint16_t tile_id;

	EXPECT_TRUE( !test_->EnterTile(read_id_string, &returned_error_message) || !test_->tile_accessible_ ) << error_msg;
	EXPECT_EQ( expected_error, returned_error_message.str() ) << comment;

	if( !test_->tile_accessible_ ){
		EXPECT_TRUE( test_->GetTileId(tile_id, read_id_string) ) << error_msg;
		EXPECT_EQ( 0, tile_id ) << "tile_accessible_ checking doesn't work properly\n"; // EnterTile might have set tile_accessible_ to false and then we set tile_id simply to 0 in GetTileId
		test_->tile_accessible_ = true; // Now we set it manually to true to test the other error cases
	}
	EXPECT_TRUE( test_->GetTileId(tile_id, read_id_string) ) << error_msg;
	EXPECT_EQ( expected_id, tile_id ) << error_msg;
	EXPECT_EQ( 0, test_->tiles_[tile_id] ) << error_msg;
}

void TileStatsTest::TearDown(){
	BasicTestClass::TearDown();
	DeleteTestObject();
}

void TileStatsTest::TestSrr490124Equality(const TileStats &test, const char *context, bool test_tile_information){
	EXPECT_EQ(1, test.tiles_.size()) << "SRR490124-4pairs tiles_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.tiles_[0]) << "SRR490124-4pairs tiles_ wrong for " << context << '\n';
	EXPECT_EQ(1, test.tile_ids_.size()) << "SRR490124-4pairs tile_ids_ wrong for " << context << '\n';
	EXPECT_EQ(0, test.tile_ids_.at(1)) << "SRR490124-4pairs tile_ids_ wrong for " << context << '\n';
	if(test_tile_information){
		EXPECT_EQ(2, test.tile_colon_number_) << "SRR490124-4pairs tile_colon_number_ wrong for " << context << '\n';
		EXPECT_TRUE(test.tile_accessible_) << "SRR490124-4pairs tile_accessible_ wrong for " << context << '\n';
	}

	EXPECT_EQ(1, test.abundance_.size()) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';
	EXPECT_EQ(4, test.abundance_[0]) << "SRR490124-4pairs abundance_ wrong for " << context << '\n';
}

void TileStatsTest::TestTiles(const TileStats &test){
	EXPECT_EQ(2, test.abundance_.size()) << "abundance_ wrong in tile test\n";
	EXPECT_EQ(2, test.abundance_[0]) << "abundance_ wrong in tile test\n";
	EXPECT_EQ(2, test.abundance_[1]) << "abundance_ wrong in tile test\n";
}

namespace reseq{
	TEST_F(TileStatsTest, TileIdentification){
		CreateTestObject();

		// GetTileId: old Illumina format
		TestProperTileFormat("HWUSI-EAS100R:6:73:941:1973#0/1", "old Illumina", 0, 73, 2);
		TestProperTileFormatWithGivenColonNumber("HWUSI-EAS100R:6:73:941:1973#0/1", "old Illumina", 0);
		TestErroneousTileFormat(
				"EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG",
				false,
				1,
				"'GetTileId' did not notice tile id format change\n",
				"The read id encoding changed and the tile could not be found at colon number 2 anymore. Tile will be treated as 0: EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG\nException thrown: Unable to convert 'FC706VJ' into unsigned short.\n"
				);
		TestProperTileFormat("HWUSI-EAS100R:6:73:941:1973", "shortened old Illumina", 0, 73, 2);

		// GetTileId: new Illumina format
		TestProperTileFormat("EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG", "new Illumina", 2, 2104, 4);
		TestProperTileFormatWithGivenColonNumber("EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG", "new Illumina", 2);

		// GetTileId: NCBI format
		TestProperTileFormat("SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36", "NCBI", 3, 1, 2);
		TestProperTileFormatWithGivenColonNumber("SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36", "NCBI", 3);
		TestProperTileFormat("SSRR490124.258144_1706:3:1:20328:165629", "bam converted NCBI", 3, 1, 2);

		// Error detection
		TestErroneousTileFormat(
				"SomeReadIdWithoutAColon",
				false,
				1,
				"Failed for read id without a colon with given colon number\n",
				"The read id encoding changed and the id ended before colon number 3. Tile will be treated as 0: SomeReadIdWithoutAColon\n"
				);
		TestErroneousTileFormat(
				"Some read id without a colon",
				true,
				1,
				"Failed for read id without a colon without given colon number\n",
				"Exception: Premature end of read after first field."
				);
		TestErroneousTileFormat(
				"EAS139:i36:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG",
				true,
				1,
				"Failed to notice non-integer second field\n",
				"Exception: Failed to convert second field into an integer."
				);
		TestErroneousTileFormat(
				"EAS139:136;FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG",
				true,
				1,
				"Failed to notice wrong second delimiter\n",
				"Exception: Wrong delimiter: ';' found instead of ':' at second occurrence."
				);
		TestErroneousTileFormat(
				"SRR001666.1 071112_SLXA-EAS1_s_7:5:1J:817:345 length=36",
				true,
				1,
				"Failed to notice non-integer third field\n",
				"Exception: Could not convert third field to int: 1J"
				);
		TestErroneousTileFormat(
				"EAS139:136:FC706VJ;2:2104:15343:197393 1:Y:18:ATCACG",
				true,
				1,
				"Failed to notice wrong third delimiter\n",
				"Exception: Wrong delimiter: ' ' found instead of ':' at sixth occurrence."
				);
		TestErroneousTileFormat(
				"EAS139:136:FC706VJ:2:2104:153k43:197393 1:Y:18:ATCACG",
				true,
				1,
				"Failed to notice non-integer sixth field\n",
				"Exception: Wrong delimiter: 'k' found instead of ':' at sixth occurrence."
				);
		TestErroneousTileFormat(
				"EAS139:136:FC706VJ:2:2104:15343.197393 1:Y:18:ATCACG",
				true,
				1,
				"Failed to notice wrong sixth delimiter\n",
				"Exception: Wrong delimiter: '.' found instead of ':' at sixth occurrence."
				);
		TestErroneousTileFormat(
				"SRR001666.1",
				true,
				1,
				"Failed to notice premature end of read\n",
				"Exception: Premature end of read after first field."
				);
		TestErroneousTileFormat(
				"SRR001666.1:",
				true,
				1,
				"Failed to notice premature end of read id with ':'-delimiter\n",
				"Exception: Premature end of read after second field."
				);
	}
}
