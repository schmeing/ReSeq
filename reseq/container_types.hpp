#ifndef RESEQ_CONTAINER_TYPES_H
#define RESEQ_CONTAINER_TYPES_H

#include <array>
#include <vector>

#include "constants.hpp"
#include "utilities.hpp"

namespace reseq {

// Atomic accumulator vector (thread-safe, filled during multi-threaded BAM pass)
template <typename T> using AtomicVec = std::vector<utilities::VectorAtomic<T>>;

// Domain-indexed array aliases
template <typename T> using PerSegment = std::array<T, kTemplateSegments>;

template <typename T> using PerStrand = std::array<T, kStrands>;

template <typename T> using PerBase = std::array<T, kNumBases>;

template <typename T> using PerBaseN = std::array<T, kNumBasesN>;

// Composite aliases
template <typename T> using PerSegmentPerBase = PerSegment<PerBase<T>>;

template <typename T> using PerSegmentPerBaseN = PerSegment<PerBaseN<T>>;

template <typename T> using PerSegmentPerStrand = PerSegment<PerStrand<T>>;

template <typename T> using PerSegmentPerStrandPerBase = PerSegment<PerStrand<PerBase<T>>>;

} // namespace reseq

#endif // RESEQ_CONTAINER_TYPES_H
