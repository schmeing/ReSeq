#ifndef TILESTATS_H
#define TILESTATS_H

#include <sstream>
#include <stdint.h>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/serialization/vector.hpp>

#include <seqan/stream.h>

#include "utilities.hpp"

namespace reseq{
	class TileStats{
	private:
		// Definitions
		const uint16_t kReadIdBufferSize = 200; // Size of the buffer used for parsing the read id; Has to be long enough to contain everything until the first colon in read id; Minimum 70 so some generated error messages fit

		// Variables for Simulation
		std::vector<uintTile> tiles_; // Array that contains the tile numbers at their corresponding ID
		std::unordered_map<uintTile,uintTileId> tile_ids_; // Map that contains the IDs corresponding to each tile

		std::vector<uintFragCount> abundance_; // abundance_[tile_id] = #pairs

		// Temporary variables for read in
		uint16_t tile_colon_number_; // The count at which colon beginning from the left the tile number is starting (and then goes to the next colon)
		bool tile_accessible_; // Whether the read id can be parsed to extract the tile
		
		// Helper functions
		inline void ReadUntil(std::stringstream &id_stream, char *buffer, uint16_t buffer_size, char until, const std::string &field) const;
		inline void CheckDelimiter(char delimiter, char expected_delimiter, const std::string &occurence) const;
		inline void CheckDelimiter(std::stringstream &id_stream, char expected_delimiter, const std::string &occurence) const;
		inline uint32_t ReadInt(std::stringstream &id_stream, const std::string &field, bool read_must_continue=true) const;
		
		uintTile GetKnownTile(const seqan::CharString &read_id, uint16_t tile_colon_number, bool print_warnings, std::stringstream *error_message=NULL) const;

		// Boost archive functions
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int UNUSED(version)){
			ar & tiles_;
			ar & abundance_;
			
			if(!tile_ids_.size()){
				// As tile_ids_ isn't set yet serialization is used for loading: Recreate tile_ids_ from tiles_
				for(auto id=tiles_.size(); id--; ){
					tile_ids_.emplace(tiles_.at(id), id);
				}
			}
		}

		// Google test
		friend class TileStatsTest;
	public:
		//  Constructor
		TileStats();
	
		// Getter functions
		bool TileId(uintTileId &id, uintTile tile) const;
		const std::vector<uintTileId> &Tiles() const{ return tiles_; }
		uintTileId NumTiles() const{ return tiles_.size(); }
		const std::vector<uintFragCount> &Abundance() const{ return abundance_; }
		
		// Setter functions
		void IgnoreTiles(){ tile_accessible_ = false; }
		inline void ResetTileIdAccessionInfo(){ // Resets tile_colon_number_ and tile_accessible_
			tile_colon_number_ = 0;
			tile_accessible_ = true;
		}
	
		// Main functions

		uintTile GetTile(const seqan::CharString &read_id, uint16_t &tile_colon_number, bool &tile_accessible, std::stringstream *error_message=NULL) const;
		bool EnterTile(const seqan::CharString &read_id, std::stringstream *error_message=NULL); // Reads the tile from the read id and returns its tile id. In case the tile is new it enters them into tiles_ and tile_ids_
		bool GetTileId(uintTileId &tile_id, const seqan::CharString &read_id) const;
		void Shrink();
	};
}

#endif // TILESTATS_H
