#ifndef TYPES_H
#define TYPES_H

#include <atomic>
#include <cstdint>
#include <stdexcept>

#include <seqan/bam_io.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/version.h>

// SeqAn 2.5+ renamed namespace from seqan to seqan2
#if defined(SEQAN_VERSION_MINOR) && SEQAN_VERSION_MINOR >= 5
namespace seqan = seqan2;
#endif

// https://stackoverflow.com/questions/3599160/how-to-suppress-unused-parameter-warnings-in-c
#ifdef __GNUC__
#define UNUSED(x) UNUSED_##x __attribute__((__unused__))
#else
#define UNUSED(x) UNUSED_##x
#endif

namespace reseq {
// Type definitions
using intQualDiff = int8_t; // Difference of quality values

using uintInDelType = uint8_t; // InDel Type (0,1)
using uintPercent = uint8_t;   // Percent values
using uintQual = uint8_t;      // Quality values
using uintTempSeq = uint8_t;   // Template segment, strand (0,1)

using uintPercentShift = int16_t; // Percent values with direction information for shifts

using uintAdapterId = uint16_t;      // Adapter id
using uintAlleleId = uint16_t;       // Allele number, id, etc.
using uintBaseCall = uint16_t;       // Base call in case it is converted from seqan::Dna5(or other) to int
using uintDupCount = uint16_t;       // Count of fragments at a given site
using uintErrorCount = uint16_t;     // Count of errors that occurred during execution
using uintInDelTypePrint = uint16_t; // InDel Type (0,1) for printing (so they are not printed as characters)
using uintMarginId = uint16_t;       // Id or size for probability margin or dimension
using uintNumThreads = uint16_t;     // Id for probability margin
using uintPercentPrint = uint16_t;   // Percent values for printing (so they are not printed as characters)
using uintQualPrint = uint16_t;      // Quality values for printing (so they are not printed as characters)
using uintReadLen = uint16_t;        // Position on read, length of sequence
using uintSurBlockId = uint16_t;     // Surrounding block number, id, etc.
using uintSurPos = uint16_t;         // Position in surrounding
using uintTempSeqPrint =
    uint16_t;                // Template segment, strand (0,1) for printing (so they are not printed as characters)
using uintTile = uint16_t;   // Tile encoding (e.g. 2308)
using uintTileId = uint16_t; // Tile ids

using intSeqShift = int32_t;  // Sequence length with direction information for shifts
using intVariantId = int32_t; // Variant id, number, etc. (id can be invalid < 0)

using uintCovCount = uint32_t;    // Count for position coverage
using uintMatrixIndex = uint32_t; // Index or size for probability matrix
using uintNumFits = uint32_t;     // Number of fits, bias sum parameters(ref seq, insert length), iterations, parameters
using uintReadLenCalc = uint32_t; // Calculations on read length producing temporary larger values, like calculating the
                                  // mean of something over all positions
using uintRefSeqBin = uint32_t;   // Reference sequence bin/block id, number, etc.
using uintRefSeqId = uint32_t;    // Reference sequence id, number, etc.
using uintSeqLen = uint32_t;      // Position on reference sequence, length of sequence, distance

using intExtSurrounding = int64_t; // Extended surrounding value to combine blocks, etc.
using intFragCountShift = int64_t; // General count of fragments/reads with direction information for shifts
using intSeqPos = int64_t;         // Position on reference sequence that can be invalid (-1)

using uintAlleleBitArray = uint64_t; // bit array for storing yes/no for 64 alleles
using uintFragCount = uint64_t;      // General count of fragments/reads
using uintMatrixCount = uint64_t;    // Count in probability matrix independent of origin (at least as big as
                                     // uintFragCount, uintNucCount)
using uintNucCount = uint64_t;       // General count of fragments/reads
using uintRefLenCalc = uint64_t; // Calculations with reference sequences producing temporary large values, like summing
                                 // up length of reference sequences or getting a mean
using uintSeed = uint64_t;       // Seed value for random number generator

#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the seqan stuff
namespace utilities {
// Type definitions
using CigarString = seqan::String<seqan::CigarElement<>>;
using ModComplementIupac = seqan::ModView<seqan::FunctorComplement<seqan::Iupac>>;

using ConstDna5StringReverseComplement =
    const seqan::ModifiedString<seqan::ModifiedString<const seqan::Dna5String, seqan::ModComplementDna5>,
                                seqan::ModReverse>;
using ConstDnaStringReverseComplement =
    const seqan::ModifiedString<seqan::ModifiedString<const seqan::DnaString, seqan::ModComplementDna>,
                                seqan::ModReverse>;
using ConstIupacStringReverseComplement =
    const seqan::ModifiedString<seqan::ModifiedString<const seqan::IupacString, ModComplementIupac>, seqan::ModReverse>;

using ComplementedConstDna5String = const seqan::ModifiedString<const seqan::Dna5String, seqan::ModComplementDna5>;

using ReversedConstCharString = const seqan::ModifiedString<const seqan::CharString, seqan::ModReverse>;
using ReversedConstCigarString = const seqan::ModifiedString<const CigarString, seqan::ModReverse>;

template <typename T> struct VectorAtomic {
    typename std::atomic<T> value_;

    VectorAtomic() : value_(0) {}

    VectorAtomic(const VectorAtomic& UNUSED(right)) {
        throw std::runtime_error(
            "This function should never be called. Did you resize an object with a non-zero length?");
    }

    operator T() const { return value_; }

    template <typename U> inline VectorAtomic<T>& operator=(const U& value) {
        value_ = value;
        return *this;
    }

    template <typename U> bool operator==(const U& comp) const { return value_ == comp; }

    inline VectorAtomic<T>& operator++() {
        ++value_;
        return *this;
    }

    inline T operator++(int) { return value_++; }

    template <typename U> inline T operator+=(const U& rhs) { return value_ += rhs; }

    inline VectorAtomic<T>& operator--() {
        --value_;
        return *this;
    }

    inline T operator--(int) { return value_--; }

    template <typename U> inline T operator-=(const U& rhs) { return value_ -= rhs; }
};
// In C++20, the member operator== supports reversed arguments automatically
// (i.e. `5 == vectorAtomic` resolves to `vectorAtomic.operator==(5)`).
// A free-function reverse delegation (`return rhs == lhs`) causes infinite
// recursion because C++20 overload resolution can pick it as the reversed
// candidate for its own body.

} // namespace utilities
#endif // SWIG
} // namespace reseq

#endif // TYPES_H
