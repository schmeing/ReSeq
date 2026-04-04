#include "archive_format.h"

namespace reseq::format {

ArchiveFormat DetectFormat(std::istream& is, const char (&magic)[3]) {
    char header[kHeaderSize];
    is.read(header, kHeaderSize);

    if (is.gcount() < static_cast<std::streamsize>(kHeaderSize)) {
        // File too short for a binary header — rewind and treat as text
        is.clear();
        is.seekg(0);
        return ArchiveFormat::kText;
    }

    if (header[0] == magic[0] && header[1] == magic[1] && header[2] == magic[2]) {
        // Magic matches — check version
        auto version = static_cast<uint8_t>(header[3]);
        if (version == kFormatVersionCompressedBinary) {
            return ArchiveFormat::kCompressedBinaryV1;
        }
        // Known magic, unknown version — do NOT fall back to text
        return ArchiveFormat::kUnsupportedBinaryVersion;
    }

    // No magic — rewind for text archive parser
    is.clear();
    is.seekg(0);
    return ArchiveFormat::kText;
}

void WriteHeader(std::ostream& os, const char (&magic)[3]) {
    os.write(magic, 3);
    char version = static_cast<char>(kFormatVersionCompressedBinary);
    os.write(&version, 1);
}

} // namespace reseq::format
