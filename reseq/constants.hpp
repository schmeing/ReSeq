#ifndef RESEQ_CONSTANTS_H
#define RESEQ_CONSTANTS_H

#include "utilities.hpp"

namespace reseq {

inline constexpr uintBaseCall kNumBases = 4;        // A, C, G, T
inline constexpr uintBaseCall kNumBasesN = 5;       // A, C, G, T, N
inline constexpr uintTempSeq kTemplateSegments = 2; // first, second read
inline constexpr uintTempSeq kStrands = 2;          // forward, reverse
inline constexpr uintQual kPhredSangerOffset = 33;
inline constexpr uintQual kPhredIlluminaOffset = 64;
inline constexpr uint16_t kGCBins = 101; // 0-100% inclusive

} // namespace reseq

#endif // RESEQ_CONSTANTS_H
