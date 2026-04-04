#ifndef ARCHIVE_FORMAT_H
#define ARCHIVE_FORMAT_H

#include <cstdint>
#include <istream>
#include <ostream>

namespace reseq::format {

// Magic bytes identifying binary format families
constexpr char kStatsMagic[3] = {'R', 'S', 'Q'};
constexpr char kIpfMagic[3] = {'I', 'P', 'F'};

// Format version: 1 = gzip-compressed Boost binary archive
constexpr uint8_t kFormatVersionCompressedBinary = 1;

// Total header size: 3-byte magic + 1-byte version
constexpr size_t kHeaderSize = 4;

enum class ArchiveFormat {
    kText,                    // Legacy Boost text archive (no header)
    kCompressedBinaryV1,      // Gzip-compressed Boost binary archive
    kUnsupportedBinaryVersion // Magic matches but version is unknown
};

/// Peek first kHeaderSize bytes from stream.
/// If magic matches and version is known, return kCompressedBinaryV1
///   with stream positioned after header.
/// If magic matches but version is unknown, return kUnsupportedBinaryVersion
///   (caller must emit error — do NOT fall through to text).
/// If magic does not match, rewind stream and return kText.
ArchiveFormat DetectFormat(std::istream& is, const char (&magic)[3]);

/// Write magic + version header to stream.
void WriteHeader(std::ostream& os, const char (&magic)[3]);

} // namespace reseq::format

#endif // ARCHIVE_FORMAT_H
