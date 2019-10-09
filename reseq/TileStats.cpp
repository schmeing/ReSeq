#include "TileStats.h"
using reseq::TileStats;

#include <exception>
using std::exception;
//include <sstream>
using std::stringstream;
#include <stdexcept>
using std::runtime_error;
#include <stdio.h>
#include <stdlib.h>
//include <string>
using std::char_traits;
using std::string;
//include <vector>
using std::vector;

#include "reportingUtils.hpp"

//include <seqan/stream.h>
using seqan::CharString;
using seqan::Exception;
using seqan::lexicalCast;

inline void TileStats::ReadUntil(stringstream &id_stream, char *buffer, uint16_t buffer_size, char until, const string &field) const{
	id_stream.getline( buffer, buffer_size, until);
	if( id_stream.eof() ){
		throw runtime_error(string("Premature end of read after ") + field + " field.");
	}
}

inline void TileStats::CheckDelimiter(char delimiter, char expected_delimiter, const string &occurence) const{
	if( expected_delimiter != delimiter ){
		throw runtime_error( string("Wrong delimiter: '") + delimiter + "' found instead of '" + expected_delimiter + "' at " + occurence + " occurrence." );
	}
}

inline void TileStats::CheckDelimiter(stringstream &id_stream, char expected_delimiter, const string &occurence) const{
	char delimiter;
	id_stream.get(delimiter); // Need to use get, because << automatically skips ' '
	CheckDelimiter(delimiter, expected_delimiter, occurence);
}

inline uint32_t TileStats::ReadInt(stringstream &id_stream, const string &field, bool read_must_continue) const{
	uint32_t int_value;
	id_stream >> int_value;
	if( read_must_continue && id_stream.eof() ){
		throw runtime_error(string("Premature end of read after ") + field + " field.");
	}
	if( id_stream.fail() ){
		throw runtime_error(string("Failed to convert ") + field + " field into an integer.");
	}

	return int_value;
}

reseq::uintTile TileStats::GetKnownTile(const CharString &read_id, uint16_t tile_colon_number, bool print_warnings, stringstream *error_message) const{
	uint16_t pos(0);

	for(uint16_t num_colons = 0; num_colons != tile_colon_number && pos < length(read_id); ){ // Stop at colon tile_colon_number
		if( ':' == read_id[pos++] ) ++num_colons; // Count number of colons
	}

	uint16_t start_pos(pos);

	while( ':' != read_id[++pos] && pos < length(read_id)); // Stop at next colon

	if( pos < length(read_id) ){
		try{
			return lexicalCast<uintTile>(infix(read_id, start_pos, pos));
		}
		catch(const Exception &e){
			if(print_warnings){
				(error_message ? *error_message : printWarn) <<  "The read id encoding changed and the tile could not be found at colon number " << tile_colon_number << " anymore. Tile will be treated as 0: " << read_id << "\nException thrown: " << e.what() << std::endl;
			}
			return 0;
		}
	}
	else{
		if(print_warnings){
			(error_message ? *error_message : printWarn) <<  "The read id encoding changed and the id ended before colon number " << tile_colon_number+1 << ". Tile will be treated as 0: " << read_id << std::endl;
		}
		return 0;
	}
}

TileStats::TileStats():
	tile_colon_number_(0), // 0 indicates that the correct tile colon number has not been found yet
	tile_accessible_(true) // Tile number is assumed to be accessible before the first read has been read
	{
}

bool TileStats::TileId(uintTileId &id, uintTile tile) const{
	// Determine tile id
	auto it = tile_ids_.find(tile);
	if(tile_ids_.end() == it){ // Tile hasn't been found
		return false;
	}
	else{ // Tile was found already in a previous read
		id = it->second;
		return true;
	}
}

reseq::uintTile TileStats::GetTile(const CharString &read_id, uint16_t &tile_colon_number, bool &tile_accessible, stringstream *error_message) const{
	uint32_t tile(0); // Make it a bigger int, because it also has to accept a position info to check for integer convertibility of that string

	if(tile_colon_number){ // The colon number where to find the tile is know from previous reads
		tile = GetKnownTile(read_id, tile_colon_number, true, error_message);
	}
	else if(tile_accessible){ // First read and tile colon number has to be determined
		// Formats to check for:
		// Old Illumina format:						@HWUSI-EAS100R:6:73:941:1973#0/1							@(unique instrument name):(flowcell lane):(tile number):(x-coordinate):(y-coordinate)#(index number)/(member of pair)
		// New Illumina format since Casava 1.8:	@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG		@(unique instrument name):(run id):(flowcell id):(flowcell lane):(tile number):(x-coordinate):(y-coordinate) (member of pair):(filtered?):(control bits):(index)
		// NCBI Sequence Read Archive:				@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36		@(NCBI identifier) (original format identifier as above) length=(read length)

		stringstream id_stream(toCString(read_id));

		try{
			char buffer[kReadIdBufferSize];

			 // read until the first colon (unique instrument name) + potentially (NCBI identifier)
			ReadUntil(id_stream, buffer, kReadIdBufferSize, ':', "first");

			// flowcell lane or run id: Has to be int
			ReadInt(id_stream, "second");
			CheckDelimiter(id_stream, ':', "second");

			// tile number or flowcell id: Might contain letters
			ReadUntil(id_stream, buffer, kReadIdBufferSize, ':', "third");

			// x-coordinate or flowcell lane: Has to be int
			ReadInt(id_stream, "fourth");
			CheckDelimiter(id_stream, ':', "fourth");

			// y-coordinate or tile number: Has to be int
			tile = ReadInt(id_stream, "fifth", false);

			char delimiter = id_stream.get();
			if( '#' == delimiter || ' ' == delimiter || '_' == delimiter || char_traits<char>::eof() == delimiter ){ // Old Illumina format
				// Member of pair and index number are not checked for since they might not exist
				tile = atoi(buffer); // Tile number have been stored in buffer before
				string tileString(buffer);
				if( tile && strlen(buffer) == snprintf( buffer, kReadIdBufferSize, "%u", tile ) ){ // Make sure atoi found digits to convert and all characters it found where digits
					tile_colon_number = 2;
				}
				else{
					throw runtime_error(string("Could not convert third field to int: ")+tileString);
				}
			}
			else{ // New Illumina format
				CheckDelimiter(delimiter, ':', "fifth");

				// x-coordinate
				ReadInt(id_stream, "sixth", false);
				CheckDelimiter(id_stream, ':', "sixth");

				// y-coordinate
				ReadInt(id_stream, "seventh", false);

				delimiter = id_stream.get();
				if( ' ' != delimiter &&  '_' != delimiter && char_traits<char>::eof() != delimiter ){
					throw runtime_error(string("Delimiter '") + delimiter + "' is non of the possible delimiters at seventh occurence (' ','_',eof)");
				}

				// Not checking for the part after the space delimiter, because it is frequently changed by data processing software like error correction
				// Tile is already set correctly
				tile_colon_number = 4;
			}
		}
		catch(const exception & e){
			if( error_message ){
				*error_message << "Exception: " << e.what();
			}
			else{
				printWarn << "The read enconding is unknown. All tiles will be treated as 0: " << read_id << "\nException thrown: " << e.what() << std::endl;
			}
			tile=0; // Not exiting because tile=0 has still to be entered into the dictionary.
			tile_accessible=false; // Don't try further to read tile from read id
		}
	}

	return tile;
}

bool TileStats::EnterTile(const CharString &read_id, stringstream *error_message){
	uintTile tile = GetTile(read_id, tile_colon_number_, tile_accessible_, error_message);

	// Determine tile id
	auto it = tile_ids_.find(tile);
	if(tile_ids_.end() == it){ // Tile hasn't been found
		// Enter tile + id into map and lookup table
		tile_ids_.emplace(tile, tile_ids_.size());
		tiles_.push_back(tile);
		abundance_.push_back(1);
	}
	else{ // Tile was found already in a previous read
		++abundance_.at(it->second);
	}

	// If we have a tile everything is fine, if tiles have been set to being ignored everything is fine, if tiles are found to be inaccessible not everything is fine, but we still want to continue
	// Only if we already found accessible tiles and then find inaccessible tiles we want to quit
	return tile || !tile_accessible_;
}

bool TileStats::GetTileId(uintTileId &tile_id, const CharString &read_id) const{
	if(tile_accessible_){
		auto tile = GetKnownTile(read_id, tile_colon_number_, false);

		// Determine tile id
		auto it = tile_ids_.find(tile);
		if(tile_ids_.end() == it){ // Tile hasn't been found
			return false;
		}
		else{ // Tile was found already in a previous read
			tile_id = it->second;
			return true;
		}
	}
	else{
		tile_id = 0;
		return true;
	}
}

void TileStats::Shrink(){
	tiles_.shrink_to_fit();
	abundance_.shrink_to_fit();
}
