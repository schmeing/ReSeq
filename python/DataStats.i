/* File : example.i */
%module DataStats
#pragma SWIG nowarn=389

%{
#include <stdint.h>
namespace reseq{
	uint16_t kVerbosityLevel = 99;
	bool kNoDebugOutput = false;
}
#include "utilities.hpp"
#include "Vect.hpp"
#include "Surrounding.h"
#include "Reference.h"
#include "DataStatsInterface.h"
%}

%include <std_array.i>
%include <stdint.i>
%include <std_map.i>
%include <std_pair.i>
%include <std_vector.i>
typedef unsigned long int size_t;
%include "utilities.hpp"

%template(CharVector) std::vector<signed char>;
%template(CharHistVector) std::pair< size_t, std::vector<signed char> >;
%template(UCharVector) std::vector<unsigned char>;
%template(UCharHistVector) std::pair< size_t, std::vector<unsigned char> >;
%template(UInt16Vector) std::vector<uint16_t>;
%template(UInt64Vector) std::vector<uint64_t>;
%template(UInt64Array1024) std::array<uint64_t, 1024>;
%template(UInt16HistVector) std::pair< size_t, std::vector<uint16_t> >;
%template(UInt64HistVector) std::pair< size_t, std::vector<uint64_t> >;
%template(DoubleVector) std::vector<double>;
%template(DoubleHistVector) std::pair< size_t, std::vector<double> >;
%template(UInt32Pair) std::pair<uint32_t,uint32_t>;
%template(UInt64MapUInt32Pair) std::map<std::pair<uint32_t,uint32_t>, uint64_t>;

%include "SurroundingBase.hpp"
%template(SurroundingBaseSurrounding) reseq::SurroundingBase<3,10,10,int32_t>;
%include "Surrounding.h"
%include "Reference.h"
/* Let's just grab the original header file here */
%include "DataStatsInterface.h"
