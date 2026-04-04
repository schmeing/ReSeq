# Phase 6: Binary Serialization & Profile Optimization — Design Spec

**Date:** 2026-04-03
**Status:** Approved
**Prerequisite:** Phase 5a (main.cpp decomposition) complete

---

## Problem

ReSeq stores sequencing profiles (`.reseq`) and probability estimates (`.reseq.ipf`) using Boost text archives. These files are:

- **Large:** 85-98MB per `.reseq` profile, 14-55MB per `.reseq.ipf`
- **Slow to load:** ~660ms to deserialize a 98MB text profile (dominates `queryProfile` runtime)
- **Inefficient:** text format requires string↔number conversion for millions of values

Five profiles exist in the upstream [schmeing/ReSeq-profiles](https://github.com/schmeing/ReSeq-profiles) repository, all in text format. Only one (`Hs-Nova-TruSeq`) is currently on Zenodo.

## Goal

Switch serialization to gzip-compressed Boost binary archives. Maintain full backward compatibility with existing text-format files. Publish all 5 profiles in both formats on Zenodo.

---

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Format detection | Magic header + auto-detect on load | Zero friction — old files just work |
| New format | Compressed binary (gzip + binary archive) | Maximum size reduction in one step |
| Uncompressed binary | Skipped — never shipped | No reason for intermediate format |
| `--textFormat` flag | On all writing commands | Flexible, minimal wiring |
| `convertProfile` command | New top-level command | Follows Phase 5a pattern |
| File extensions | Same `.reseq` / `.reseq.ipf` | Auto-detect makes extension irrelevant |
| Platform dependence | Accept (x86_64), document | Text format is portable fallback |
| Magic bytes / versions | Named constants in `reseq/archive_format.h` | No hardcoded literals; core infra, not CLI |
| Zenodo release | New version of existing record | Single DOI, old version remains accessible |
| Compression library | boost::iostreams + gzip | zlib already a dependency, iostreams in same Boost distribution |

---

## Architecture

### Archive Format Constants

New header `reseq/archive_format.h` (core serialization infrastructure, not CLI-specific):

```cpp
#ifndef ARCHIVE_FORMAT_H
#define ARCHIVE_FORMAT_H

#include <cstdint>
#include <istream>
#include <ostream>

namespace reseq::format {

// Magic bytes identifying binary format families
constexpr char kStatsMagic[3] = {'R', 'S', 'Q'};
constexpr char kIpfMagic[3]   = {'I', 'P', 'F'};

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
```

**Detection logic:** Text archives always start with `"22 serialization::archive"` (ASCII digits). The magic bytes `RSQ` / `IPF` are not valid Boost text preamble, so detection is unambiguous. When magic matches but version is unrecognized, `DetectFormat` returns `kUnsupportedBinaryVersion` so the caller can emit a clear error ("unsupported profile format version N — please upgrade ReSeq") rather than silently misinterpreting the data.

### Implementation in `archive_format.cpp`

```cpp
ArchiveFormat DetectFormat(std::istream& is, const char (&magic)[3]) {
    char header[kHeaderSize];
    is.read(header, kHeaderSize);

    if (is.gcount() == kHeaderSize &&
        header[0] == magic[0] &&
        header[1] == magic[1] &&
        header[2] == magic[2]) {
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
```

---

## Load Path (Auto-Detect)

Both `DataStats::Load()` and `ProbabilityEstimates::Load()` follow this pattern:

```cpp
bool DataStats::Load(const char* archive_file) {
    ifstream ifs(archive_file, ios::binary);
    auto fmt = reseq::format::DetectFormat(ifs, reseq::format::kStatsMagic);

    switch (fmt) {
    case reseq::format::ArchiveFormat::kCompressedBinaryV1: {
        boost::iostreams::filtering_istream fis;
        fis.push(boost::iostreams::gzip_decompressor());
        fis.push(ifs);
        boost::archive::binary_iarchive ia(fis);
        ia >> *this;
        break;
    }
    case reseq::format::ArchiveFormat::kText: {
        boost::archive::text_iarchive ia(ifs);
        ia >> *this;
        break;
    }
    case reseq::format::ArchiveFormat::kUnsupportedBinaryVersion:
        printErr << "Unsupported profile format version in '"
                 << archive_file << "'. Please upgrade ReSeq." << std::endl;
        return false;
    }
}
```

`ProbabilityEstimates::Load()` is identical but uses `reseq::format::kIpfMagic`.

**No changes** to any `serialize()` template in any sub-object. Boost handles text↔binary transparently.

---

## Save Path (Binary Default, Text Optional)

```cpp
bool DataStats::Save(const char* archive_file, bool text_format = false) const {
    ofstream ofs(archive_file, ios::binary);
    if (text_format) {
        boost::archive::text_oarchive oa(ofs);
        oa << *this;
    } else {
        reseq::format::WriteHeader(ofs, reseq::format::kStatsMagic);
        boost::iostreams::filtering_ostream fos;
        fos.push(boost::iostreams::gzip_compressor());
        fos.push(ofs);
        boost::archive::binary_oarchive oa(fos);
        oa << *this;
    }
}
```

`ProbabilityEstimates::Save()` is identical but uses `reseq::format::kIpfMagic`.

`DataStatsInterface::Save()` passes `text_format` through.

The `text_format` parameter defaults to `false`, so all existing internal call sites produce compressed binary without code changes.

---

## CLI Changes

### `--textFormat` Flag

Added to the option descriptions of:

- **`illuminaPE`** (Stats group) — controls both `.reseq` and `.reseq.ipf` output
- **`seqToIllumina`** — controls `.reseq.ipf` output
- **`convertProfile`** — controls output format

For `illuminaPE`, the flag is a simple `bool text_format = opts_map.count("textFormat")` passed to `DataStats::Save()` and `ProbabilityEstimates::Save()` at the call site.

**Note on `ProbabilityEstimates::Estimate()`:** This method calls `this->Save(output)` internally (line 1219 in `ProbabilityEstimates.cpp`) after fitting completes. The `--textFormat` flag cannot be threaded through as simple CLI wiring. The design adds a `text_format` parameter to `Estimate()` itself, which it passes through to the internal `Save()` call. This is preferable to refactoring save responsibility out of `Estimate()`, since the save is tightly coupled to the precision-improvement check and error cleanup logic within that method.

### New Command: `convertProfile`

```
Usage:  reseq convertProfile -s <stats.reseq> [-p <probs.reseq.ipf>] [options]

Converts profile files between text and compressed binary formats.

Options:
  -s  Input stats file (.reseq)
  -o  Output stats file [overwrites input if omitted]
  -p  Input probabilities file (.reseq.ipf)       [optional]
  -P  Output probabilities file [overwrites input if omitted]
  --textFormat  Write legacy text format (default: compressed binary)
```

**Files:**
- `reseq/cli/convert_profile.h` + `reseq/cli/convert_profile.cpp`
- `reseq/main.cpp` — add dispatch entry
- `reseq/CMakeLists.txt` — add to `target_sources`

**Implementation:** Load file (auto-detects format) → Save in target format. When `-o` / `-P` are omitted, write to a temporary file in the same directory, then atomically rename over the input (safe in-place conversion).

---

## Build Changes

### CMakeLists.txt

Add `iostreams` to Boost components:

```cmake
find_package(Boost 1.48.0 REQUIRED
  filesystem iostreams math_c99 math_c99f math_c99l
  math_tr1 math_tr1f math_tr1l program_options serialization system)
```

### New Includes (in .cpp files only)

```cpp
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
```

These go in `DataStats.cpp`, `ProbabilityEstimates.cpp`, and `reseq/archive_format.cpp` — not in headers, to minimize compile-time impact.

### CI

No changes needed. Ubuntu 24.04 CI already installs `libboost-all-dev` which includes `libboost-iostreams-dev`. zlib is already required.

---

## Verification

### Unit Tests

1. **Semantic round-trip: text → binary → reload** — Load text profile, save as compressed binary, reload from binary. Verify semantic equivalence: successful reload plus equality of representative queried outputs (maxReadLength, maxLenDeletion, TotalNumberReads, etc.). Do NOT require bitwise file identity — archive and gzip implementation details may cause benign differences.

2. **Semantic round-trip: binary → reload** — Save as binary, reload. Verify same semantic equivalence as above.

3. **Query equivalence across formats** — Load same profile in text and binary formats, verify all `queryProfile` outputs match (`maxReadLength`, `maxLenDeletion`, `fragLenBias`).

4. **Format detection** — Construct files with text preamble, valid binary header, and magic-match-but-unknown-version header. Verify `DetectFormat` returns `kText`, `kCompressedBinaryV1`, and `kUnsupportedBinaryVersion` respectively.

5. **Unsupported version error path** — Attempt to Load a file with valid magic but unknown version byte. Verify Load returns false with an error message (not a misleading parse failure).

6. **Existing DataStats tests** — 2 existing Save/Load tests in `DataStatsTest.cpp` continue to pass (now produce binary by default). Add mirrored tests that force `text_format = true`.

7. **Existing ProbabilityEstimates tests** — 4 existing Save/Load round-trip tests in `ProbabilityEstimatesTest.cpp` continue to pass (now produce binary by default). Add mirrored tests that force `text_format = true`.

### Regression Tests (RegressionTest.cpp)

6. **convertProfile exit codes** — `convertProfile -s <file>` exits 0, converts successfully.

7. **`--textFormat` flag** — Verify `illuminaPE --statsOnly --textFormat` produces a file starting with `"22 serialization::archive"`.

### Benchmark

8. **All 5 profiles** — Load time and file size in both formats. Compare against text baseline.

### Profile Validation

9. **Each of 5 profiles** — Load text, convert to binary, load binary, compare all queryable outputs. Ensures no data loss across the full profile collection.

---

## Available Profiles

All from [schmeing/ReSeq-profiles](https://github.com/schmeing/ReSeq-profiles):

| Profile | Sequencer | Species | Read len | Text .reseq | Text .ipf |
|---------|-----------|---------|----------|-------------|-----------|
| `Ec-Hi2000-TruSeq` | HiSeq 2000 | *E. coli* | 100bp | 87MB | 31MB |
| `Ec-Hi4000-Nextera` | HiSeq 4000 | *E. coli* | 151bp | 85MB | 14MB |
| `Bc-Hi4000-Nextera` | HiSeq 4000 | *B. cereus* | 151bp | 88MB | 25MB |
| `Hs-Nova-TruSeq` | NovaSeq 6000 | *H. sapiens* | 151bp | 98MB | 53MB |
| `Ec-Mi-TruSeq` | MiSeq | *E. coli* | 251bp | 90MB | 55MB |
| **Total** | | | | **448MB** | **178MB** |

Expected compressed binary total: ~45-90MB (vs 626MB text). Exact sizes determined during implementation.

---

## Zenodo Release

Update existing Zenodo record (DOI: 10.5281/zenodo.19383555) with new version:

- All 5 profiles in both formats (20 files)
- Organized in `text/` and `binary/` subdirectories
- Updated description with profile metadata table

Update `test/download_test_data.sh`:
- Download binary format by default (smaller, faster)
- `--text` flag for text format
- Download all 5 profiles
- Verify checksums

---

## Platform Notes

Boost binary archives encode native endianness, type sizes, and ABI details. Compressed binary profiles are **only guaranteed to be portable within the project's supported native build environment / ABI family** (currently GCC and Clang on x86_64 Linux). Profiles generated on one such toolchain should load on another, but cross-architecture or cross-platform portability is not guaranteed. Text format remains the portable fallback for any environment where binary loading fails.

---

## Files Modified

| File | Change |
|------|--------|
| `reseq/archive_format.h` | **New** — format constants, enums, DetectFormat, WriteHeader |
| `reseq/archive_format.cpp` | **New** — implementation |
| `reseq/cli/convert_profile.h` | **New** — RunConvertProfile declaration |
| `reseq/cli/convert_profile.cpp` | **New** — convertProfile command |
| `reseq/DataStats.cpp` | Modify Load/Save for auto-detect + compressed binary |
| `reseq/DataStats.h` | Add `text_format` parameter to Save |
| `reseq/DataStatsInterface.h` | Pass through `text_format` |
| `reseq/DataStatsInterface.cpp` | Pass through `text_format` |
| `reseq/ProbabilityEstimates.cpp` | Modify Load/Save/Estimate for auto-detect + compressed binary |
| `reseq/ProbabilityEstimates.h` | Add `text_format` parameter to Save and Estimate |
| `reseq/cli/illumina_pe.cpp` | Add `--textFormat` option, pass to Save and Estimate |
| `reseq/cli/seq_to_illumina.cpp` | Add `--textFormat` option, pass to Estimate |
| `reseq/main.cpp` | Add convertProfile dispatch |
| `reseq/CMakeLists.txt` | Add new source files |
| `CMakeLists.txt` | Add `iostreams` to Boost components |
| `reseq/RegressionTest.cpp` | Add format conversion and detection tests |
| `reseq/DataStatsTest.cpp` | Add text-format mirrored Save/Load tests |
| `reseq/ProbabilityEstimatesTest.cpp` | Add text-format mirrored Save/Load tests |
| `test/download_test_data.sh` | Download all 5 profiles, binary default |
