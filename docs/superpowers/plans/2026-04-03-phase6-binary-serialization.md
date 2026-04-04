# Phase 6: Binary Serialization & Profile Optimization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Switch `.reseq` and `.reseq.ipf` serialization from Boost text archives to gzip-compressed binary archives with backward-compatible auto-detection, and add a `convertProfile` command.

**Architecture:** A new `reseq/archive_format.{h,cpp}` module provides format constants, magic header detection (`DetectFormat`), and header writing (`WriteHeader`). `DataStats` and `ProbabilityEstimates` Load/Save methods gain format-switching logic. A `--textFormat` CLI flag threads through to all write paths including `Estimate()`. A new `convertProfile` top-level command enables batch conversion.

**Tech Stack:** C++20, Boost.Serialization (binary_archive), Boost.Iostreams (gzip), CMake

**Design spec:** `docs/superpowers/specs/2026-04-03-phase6-binary-serialization-design.md`

---

## Task 1: Add Boost.Iostreams Dependency

**Files:**
- Modify: `CMakeLists.txt:44`

- [ ] **Step 1: Update Boost find_package**

In `CMakeLists.txt`, change line 44 from:

```cmake
find_package(Boost 1.48.0 REQUIRED filesystem math_c99 math_c99f math_c99l math_tr1 math_tr1f math_tr1l program_options serialization system)
```

to:

```cmake
find_package(Boost 1.48.0 REQUIRED filesystem iostreams math_c99 math_c99f math_c99l math_tr1 math_tr1f math_tr1l program_options serialization system)
```

- [ ] **Step 2: Verify the build still works**

```bash
make build && make test
```

Expected: all tests pass, no new warnings.

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt
git commit -m "build(6): add Boost.Iostreams dependency for gzip compression"
```

---

## Task 2: Create archive_format Module

**Files:**
- Create: `reseq/archive_format.h`
- Create: `reseq/archive_format.cpp`
- Modify: `reseq/CMakeLists.txt:10-24` (add to `reseq_lib`)

- [ ] **Step 1: Create reseq/archive_format.h**

```cpp
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
```

- [ ] **Step 2: Create reseq/archive_format.cpp**

```cpp
#include "archive_format.h"

namespace reseq::format {

ArchiveFormat DetectFormat(std::istream& is, const char (&magic)[3]) {
    char header[kHeaderSize];
    is.read(header, kHeaderSize);

    if (is.gcount() == static_cast<std::streamsize>(kHeaderSize) && header[0] == magic[0] &&
        header[1] == magic[1] && header[2] == magic[2]) {
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
```

- [ ] **Step 3: Add archive_format.cpp to reseq_lib in CMakeLists.txt**

In `reseq/CMakeLists.txt`, add `archive_format.cpp` to the `reseq_lib` sources (after `AdapterStats.cpp`):

```cmake
add_library(reseq_lib STATIC
  AdapterStats.cpp
  archive_format.cpp
  CoverageStats.cpp
  ...
```

- [ ] **Step 4: Build and test**

```bash
make build && make test
```

Expected: all tests pass. The module is compiled but not yet called.

- [ ] **Step 5: Commit**

```bash
git add reseq/archive_format.h reseq/archive_format.cpp reseq/CMakeLists.txt
git commit -m "refactor(6): add archive_format module with format detection and header writing"
```

---

## Task 3: Switch DataStats Load/Save to Binary

**Files:**
- Modify: `reseq/DataStats.h:335-336`
- Modify: `reseq/DataStats.cpp:1407-1445`
- Modify: `reseq/DataStatsInterface.h:195-196`
- Modify: `reseq/DataStatsInterface.cpp:467-478`

- [ ] **Step 1: Update DataStats.h Save signature**

In `reseq/DataStats.h`, change line 336 from:

```cpp
    bool Save(const char* archive_file) const;
```

to:

```cpp
    bool Save(const char* archive_file, bool text_format = false) const;
```

- [ ] **Step 2: Update DataStats.cpp Load**

Replace `DataStats::Load` (lines 1407-1426) with:

```cpp
bool DataStats::Load(const char* archive_file) {
    if (!FileExists(archive_file)) {
        printErr << "File '" << archive_file << "' does not exists or no read permission given." << std::endl;
        return false;
    }

    try {
        ifstream ifs(archive_file, std::ios::binary);
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
            printErr << "Unsupported profile format version in '" << archive_file
                     << "'. Please upgrade ReSeq." << std::endl;
            return false;
        }
    } catch (const exception& e) {
        printErr << "Could not load data statistics: " << e.what() << std::endl;
        return false;
    }

    return true;
}
```

- [ ] **Step 3: Update DataStats.cpp Save**

Replace `DataStats::Save` (lines 1428-1445) with:

```cpp
bool DataStats::Save(const char* archive_file, bool text_format) const {
    try {
        CreateDir(archive_file);

        ofstream ofs(archive_file, std::ios::binary);
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
    } catch (const exception& e) {
        printErr << "Could not save data statistics: " << e.what() << std::endl;
        return false;
    }

    return true;
}
```

- [ ] **Step 4: Add new includes to DataStats.cpp**

Add after the existing includes at the top of `DataStats.cpp`:

```cpp
#include "archive_format.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
```

- [ ] **Step 5: Update DataStatsInterface**

In `reseq/DataStatsInterface.h`, change line 196 from:

```cpp
    bool Save(const char* archive_file) const;
```

to:

```cpp
    bool Save(const char* archive_file, bool text_format = false) const;
```

In `reseq/DataStatsInterface.cpp`, change `DataStatsInterface::Save` (line 476-478) from:

```cpp
bool DataStatsInterface::Save(const char* archive_file) const {
    return stats_.Save(archive_file);
}
```

to:

```cpp
bool DataStatsInterface::Save(const char* archive_file, bool text_format) const {
    return stats_.Save(archive_file, text_format);
}
```

- [ ] **Step 6: Build and test**

```bash
make build && make test
```

Expected: all tests pass. Existing DataStats Save/Load tests now produce binary by default and reload successfully (Boost `serialize()` templates work identically across archive types).

- [ ] **Step 7: Commit**

```bash
git add reseq/DataStats.h reseq/DataStats.cpp reseq/DataStatsInterface.h reseq/DataStatsInterface.cpp
git commit -m "refactor(6): switch DataStats Load/Save to compressed binary with text fallback

Load auto-detects format via magic header. Save writes compressed binary
by default, text_format parameter for legacy output. Unsupported version
produces explicit error."
```

---

## Task 4: Switch ProbabilityEstimates Load/Save/Estimate to Binary

**Files:**
- Modify: `reseq/ProbabilityEstimates.h:1633-1637`
- Modify: `reseq/ProbabilityEstimates.cpp:1066-1107,1109,1219`

- [ ] **Step 1: Update ProbabilityEstimates.h signatures**

In `reseq/ProbabilityEstimates.h`, change lines 1633-1637 from:

```cpp
    bool Load(const char* archive_file);
    bool Save(const char* archive_file) const;

    bool Estimate(const DataStats& stats, uintNumFits max_iterations, double precision_aim, uintNumThreads num_threads,
                  const char* output, const char* input = "");
```

to:

```cpp
    bool Load(const char* archive_file);
    bool Save(const char* archive_file, bool text_format = false) const;

    bool Estimate(const DataStats& stats, uintNumFits max_iterations, double precision_aim, uintNumThreads num_threads,
                  const char* output, const char* input = "", bool text_format = false);
```

- [ ] **Step 2: Update ProbabilityEstimates.cpp Load**

Replace `ProbabilityEstimates::Load` (lines 1066-1088) with:

```cpp
bool ProbabilityEstimates::Load(const char* archive_file) {
    if (!FileExists(archive_file)) {
        printErr << "File '" << archive_file << "' does not exists or no read permission given." << std::endl;
        return false;
    }

    try {
        ifstream ifs(archive_file, std::ios::binary);
        auto fmt = reseq::format::DetectFormat(ifs, reseq::format::kIpfMagic);

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
            printErr << "Unsupported probability format version in '" << archive_file
                     << "'. Please upgrade ReSeq." << std::endl;
            return false;
        }
    } catch (const exception& e) {
        printErr << "Could not load probability estimates from '" << archive_file << "': " << e.what() << std::endl;
        return false;
    }

    error_during_fitting_ = false;
    precision_improved_ = false;

    return true;
}
```

- [ ] **Step 3: Update ProbabilityEstimates.cpp Save**

Replace `ProbabilityEstimates::Save` (lines 1090-1107) with:

```cpp
bool ProbabilityEstimates::Save(const char* archive_file, bool text_format) const {
    try {
        CreateDir(archive_file);

        ofstream ofs(archive_file, std::ios::binary);
        if (text_format) {
            boost::archive::text_oarchive oa(ofs);
            oa << *this;
        } else {
            reseq::format::WriteHeader(ofs, reseq::format::kIpfMagic);
            boost::iostreams::filtering_ostream fos;
            fos.push(boost::iostreams::gzip_compressor());
            fos.push(ofs);
            boost::archive::binary_oarchive oa(fos);
            oa << *this;
        }
    } catch (const exception& e) {
        printErr << "Could not save probability estimates to '" << archive_file << "': " << e.what() << std::endl;
        return false;
    }

    return true;
}
```

- [ ] **Step 4: Update Estimate signature and internal Save call**

Change the `Estimate` function signature (line 1109) from:

```cpp
bool ProbabilityEstimates::Estimate(const DataStats& stats, uintNumFits max_iterations, double precision_aim,
                                    uintNumThreads num_threads, const char* output, const char* input) {
```

to:

```cpp
bool ProbabilityEstimates::Estimate(const DataStats& stats, uintNumFits max_iterations, double precision_aim,
                                    uintNumThreads num_threads, const char* output, const char* input,
                                    bool text_format) {
```

Then change the internal Save call (line 1219) from:

```cpp
            if (!this->Save(output)) {
```

to:

```cpp
            if (!this->Save(output, text_format)) {
```

- [ ] **Step 5: Add new includes to ProbabilityEstimates.cpp**

Add after the existing includes at the top:

```cpp
#include "archive_format.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
```

- [ ] **Step 6: Build and test**

```bash
make build && make test
```

Expected: all tests pass. The 4 existing ProbabilityEstimates Save/Load round-trip tests now produce binary and reload successfully.

- [ ] **Step 7: Commit**

```bash
git add reseq/ProbabilityEstimates.h reseq/ProbabilityEstimates.cpp
git commit -m "refactor(6): switch ProbabilityEstimates Load/Save/Estimate to compressed binary

Same pattern as DataStats. Estimate() gains text_format parameter
threaded to internal Save() call."
```

---

## Task 5: Add --textFormat Flag to CLI Commands

**Files:**
- Modify: `reseq/cli/illumina_pe.cpp:242,487`
- Modify: `reseq/cli/seq_to_illumina.cpp:113`

- [ ] **Step 1: Add --textFormat to illuminaPE options**

In `reseq/cli/illumina_pe.cpp`, add `"textFormat"` to the Stats option group. Find the line:

```cpp
        ("vcfIn,v", value<string>(), "Ignore all positions with a listed variant for stats generation");
```

Add after it (inside the same `opt_desc.add_options()` chain):

```cpp
        ("textFormat", "Write profile files in legacy text format instead of compressed binary");
```

Then, after the `opts_map` is parsed, extract the flag. Find a suitable location after `notify(opts_map)` and add:

```cpp
            bool text_format = opts_map.count("textFormat");
```

Pass it to the DataStats Save call. Change line 242 from:

```cpp
                    real_data_stats.Save(stats_file.c_str());
```

to:

```cpp
                    real_data_stats.Save(stats_file.c_str(), text_format);
```

Pass it to the Estimate call. Change line 487 from:

```cpp
                    if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads,
                                                probs_out.c_str(), probs_in.c_str())) {
```

to:

```cpp
                    if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads,
                                                probs_out.c_str(), probs_in.c_str(), text_format)) {
```

- [ ] **Step 2: Add --textFormat to seqToIllumina options**

In `reseq/cli/seq_to_illumina.cpp`, add to the option description chain, after `"statsIn,s"`:

```cpp
        ("textFormat", "Write profile files in legacy text format instead of compressed binary");
```

Extract the flag after `notify(opts_map)`:

```cpp
            bool text_format = opts_map.count("textFormat");
```

Pass it to the Estimate call. Change line 113 from:

```cpp
            if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads, probs_out.c_str(),
                                        probs_in.c_str())) {
```

to:

```cpp
            if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads, probs_out.c_str(),
                                        probs_in.c_str(), text_format)) {
```

- [ ] **Step 3: Build and test**

```bash
make build && make test
```

- [ ] **Step 4: Commit**

```bash
git add reseq/cli/illumina_pe.cpp reseq/cli/seq_to_illumina.cpp
git commit -m "feat(6): add --textFormat flag to illuminaPE and seqToIllumina

Threads text_format through to DataStats::Save() and
ProbabilityEstimates::Estimate() for legacy output."
```

---

## Task 6: Add convertProfile Command

**Files:**
- Create: `reseq/cli/convert_profile.h`
- Create: `reseq/cli/convert_profile.cpp`
- Modify: `reseq/main.cpp`
- Modify: `reseq/CMakeLists.txt:38`

- [ ] **Step 1: Create reseq/cli/convert_profile.h**

```cpp
#ifndef CLI_CONVERT_PROFILE_H
#define CLI_CONVERT_PROFILE_H

#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace reseq::cli {

int RunConvertProfile(const std::vector<std::string>& args,
                      const boost::program_options::variables_map& general_opts,
                      boost::program_options::options_description& opt_desc_full);

} // namespace reseq::cli

#endif // CLI_CONVERT_PROFILE_H
```

- [ ] **Step 2: Create reseq/cli/convert_profile.cpp**

```cpp
#include "cli/convert_profile.h"

#include <cstdio>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "DataStats.h"
#include "logging.hpp"
#include "ProbabilityEstimates.h"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::cerr;
using std::exception;
using std::string;

namespace reseq::cli {

namespace {

/// Save to a temporary file then rename over the target for safe in-place conversion.
bool SafeRename(const string& tmp_path, const string& target_path) {
    std::error_code ec;
    std::filesystem::rename(tmp_path, target_path, ec);
    if (ec) {
        printErr << "Failed to rename '" << tmp_path << "' to '" << target_path << "': " << ec.message() << std::endl;
        std::remove(tmp_path.c_str());
        return false;
    }
    return true;
}

} // anonymous namespace

int RunConvertProfile(const std::vector<std::string>& args, const variables_map& general_opts,
                      options_description& opt_desc_full) {
    options_description opt_desc("convertProfile");
    opt_desc.add_options()("statsIn,s", value<string>(), "Input stats file (.reseq)")(
        "statsOut,o", value<string>(), "Output stats file [overwrites input if omitted]")(
        "probsIn,p", value<string>(), "Input probabilities file (.reseq.ipf)")(
        "probsOut,P", value<string>(), "Output probabilities file [overwrites input if omitted]")(
        "textFormat", "Write legacy text format (default: compressed binary)");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq convertProfile -s <stats.reseq> [-p <probs.reseq.ipf>] [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse convertProfile command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
        return 0;
    }

    bool text_format = opts_map.count("textFormat");
    bool has_stats = opts_map.count("statsIn");
    bool has_probs = opts_map.count("probsIn");

    if (!has_stats && !has_probs) {
        printErr << "At least one of statsIn (-s) or probsIn (-p) is required." << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (has_stats) {
        string stats_in = opts_map["statsIn"].as<string>();
        string stats_out;
        bool in_place = false;
        if (opts_map.count("statsOut")) {
            stats_out = opts_map["statsOut"].as<string>();
        } else {
            stats_out = stats_in + ".tmp";
            in_place = true;
        }

        printInfo << "Loading stats from " << stats_in << std::endl;
        DataStats data_stats(nullptr);
        if (!data_stats.Load(stats_in.c_str())) {
            return 1;
        }

        printInfo << "Saving stats to " << (in_place ? stats_in : stats_out)
                  << (text_format ? " (text format)" : " (compressed binary)") << std::endl;
        if (!data_stats.Save(stats_out.c_str(), text_format)) {
            return 1;
        }

        if (in_place && !SafeRename(stats_out, stats_in)) {
            return 1;
        }
    }

    if (has_probs) {
        string probs_in = opts_map["probsIn"].as<string>();
        string probs_out;
        bool in_place = false;
        if (opts_map.count("probsOut")) {
            probs_out = opts_map["probsOut"].as<string>();
        } else {
            probs_out = probs_in + ".tmp";
            in_place = true;
        }

        printInfo << "Loading probabilities from " << probs_in << std::endl;
        ProbabilityEstimates probs;
        if (!probs.Load(probs_in.c_str())) {
            return 1;
        }

        printInfo << "Saving probabilities to " << (in_place ? probs_in : probs_out)
                  << (text_format ? " (text format)" : " (compressed binary)") << std::endl;
        if (!probs.Save(probs_out.c_str(), text_format)) {
            return 1;
        }

        if (in_place && !SafeRename(probs_out, probs_in)) {
            return 1;
        }
    }

    printInfo << "Conversion complete." << std::endl;
    return 0;
}

} // namespace reseq::cli
```

- [ ] **Step 3: Add dispatch in main.cpp**

Add include after the existing cli includes (after line 33):

```cpp
#include "cli/convert_profile.h"
```

Add the convertProfile command to the general_usage string (after `seqToIllumina` line):

```cpp
        "  convertProfile\t" + "converts profiles between text and binary formats\n" +
```

Add dispatch block before the `else` (unknown command) branch:

```cpp
        } else if ("convertProfile" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in convertProfile mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunConvertProfile(unrecognized_opts, general_opts_map, opt_desc_full);
```

- [ ] **Step 4: Add to CMakeLists.txt**

In `reseq/CMakeLists.txt`, add `cli/convert_profile.cpp` to `target_sources` (line 38):

```cmake
target_sources(reseq PRIVATE cli/cli_common.cpp cli/convert_profile.cpp cli/illumina_pe.cpp cli/query_profile.cpp cli/replace_n.cpp cli/seq_to_illumina.cpp)
```

- [ ] **Step 5: Build and test**

```bash
make build && make test
```

- [ ] **Step 6: Commit**

```bash
git add reseq/cli/convert_profile.h reseq/cli/convert_profile.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "feat(6): add convertProfile command for text/binary profile conversion

Supports -s (stats), -p (probabilities), --textFormat flag.
In-place conversion via atomic rename when output path is omitted."
```

---

## Task 7: Add Unit Tests for Format Detection and Round-Trips

**Files:**
- Modify: `reseq/DataStatsTest.cpp:530-539,612-623`
- Modify: `reseq/ProbabilityEstimatesTest.cpp:1142,1269,1311,1348`
- Modify: `reseq/RegressionTest.cpp`

- [ ] **Step 1: Add format detection and unsupported version tests to RegressionTest.cpp**

Add before the closing `} // namespace reseq`:

```cpp
// --- Phase 6: binary serialization tests ---

TEST_F(RegressionTest, FormatDetectionText) {
    auto profile = GenerateEcoliProfile();
    // The profile generated by illuminaPE is now binary by default
    // Convert to text for testing text detection
    int rc = RunReseqExitOnly("convertProfile -s " + profile.string() + " --textFormat");
    ASSERT_EQ(0, rc);

    // Load and query — should auto-detect text format
    std::string output = RunReseqCapture("queryProfile -s " + profile.string() + " --maxReadLength");
    EXPECT_NE(std::string::npos, output.find("maxReadLength:"))
        << "Text format profile should load and query successfully";
}

TEST_F(RegressionTest, FormatDetectionBinary) {
    auto profile = GenerateEcoliProfile();
    // Profile is already binary by default — query directly
    std::string output = RunReseqCapture("queryProfile -s " + profile.string() + " --maxReadLength");
    EXPECT_NE(std::string::npos, output.find("maxReadLength:"))
        << "Binary format profile should load and query successfully";
}

TEST_F(RegressionTest, ConvertProfileRoundTrip) {
    auto profile = GenerateEcoliProfile();

    // Query original (binary)
    std::string binary_output = RunReseqCapture("queryProfile -s " + profile.string() + " --maxReadLength");

    // Convert to text
    auto text_profile = tmp_dir_ / "ecoli-text.reseq";
    int rc = RunReseqExitOnly("convertProfile -s " + profile.string() + " -o " + text_profile.string() + " --textFormat");
    ASSERT_EQ(0, rc);

    // Query text version
    std::string text_output = RunReseqCapture("queryProfile -s " + text_profile.string() + " --maxReadLength");
    EXPECT_EQ(binary_output, text_output) << "Binary and text profiles should produce identical query output";

    // Convert text back to binary
    auto roundtrip_profile = tmp_dir_ / "ecoli-roundtrip.reseq";
    rc = RunReseqExitOnly("convertProfile -s " + text_profile.string() + " -o " + roundtrip_profile.string());
    ASSERT_EQ(0, rc);

    // Query roundtrip version
    std::string roundtrip_output = RunReseqCapture("queryProfile -s " + roundtrip_profile.string() + " --maxReadLength");
    EXPECT_EQ(binary_output, roundtrip_output) << "Round-trip should preserve query output";
}

TEST_F(RegressionTest, ConvertProfileHelpExitCode) {
    int rc = RunReseqExitOnly("convertProfile --help");
    EXPECT_EQ(0, rc) << "convertProfile --help should exit 0";
}

TEST_F(RegressionTest, UnsupportedVersionError) {
    // Create a file with valid magic but unknown version byte
    auto bad_file = tmp_dir_ / "bad-version.reseq";
    {
        std::ofstream ofs(bad_file, std::ios::binary);
        ofs.write("RSQ", 3);
        char bad_version = 99;
        ofs.write(&bad_version, 1);
        ofs.write("garbage data", 12);
    }

    int rc = RunReseqExitOnly("queryProfile -s " + bad_file.string() + " --maxReadLength");
    EXPECT_NE(0, rc) << "Unsupported version should fail";

    std::string err = RunReseqCaptureStderr("queryProfile -s " + bad_file.string() + " --maxReadLength");
    EXPECT_NE(std::string::npos, err.find("Unsupported"))
        << "Error should mention unsupported version, got:\n" << err;
}
```

- [ ] **Step 2: Add text-format mirrored test to DataStatsTest.cpp**

After the existing Save/Load round-trip at line 539, add a text-format variant. Find:

```cpp
    EXPECT_EQ(0, remove(save_test_file.c_str())) << "Error deleting file: " << save_test_file << '\n';
```

(the first occurrence, around line 539) and add after it:

```cpp
    // Text format round-trip
    DeleteTestObject();
    CreateTestObject(&species_reference_);
    LoadStats(test_dir + "ecoli-SRR490124-4pairs.bam");
    string save_text_file = test_dir + "saveTestText.reseq";
    ASSERT_TRUE(test_->Save(save_text_file.c_str(), true)); // text_format = true
    TestSrr490124Equality("text save");

    DeleteTestObject();
    CreateTestObject(&species_reference_);
    ASSERT_TRUE(test_->Load(save_text_file.c_str()));
    test_->PrepareTesting();
    TestSrr490124Equality("text save and reload", false);
    EXPECT_EQ(0, remove(save_text_file.c_str())) << "Error deleting file: " << save_text_file << '\n';
```

Do the same for the second Save/Load block (adapters, around line 623).

- [ ] **Step 3: Add text-format mirrored tests to ProbabilityEstimatesTest.cpp**

For each of the 4 existing Save/Load blocks (lines 1142, 1269, 1311, 1348), add a text-format round-trip after the `EXPECT_EQ(0, remove(...))` line. Pattern for the first block (after line 1161):

```cpp
    // Text format round-trip
    ASSERT_TRUE(test_.Save(save_file.c_str(), true)); // text_format = true
    ProbabilityEstimates test_text;
    ASSERT_TRUE(test_text.Load(save_file.c_str()));
    auto text_iterations = GetIterationsQual(test_text, base);
    auto text_precision = GetPrecisionQual(test_text, base);
    EXPECT_EQ(iterations, text_iterations) << "Text round-trip iterations mismatch";
    EXPECT_EQ(precision, text_precision) << "Text round-trip precision mismatch";
    EXPECT_EQ(0, remove(save_file.c_str())) << "Error deleting file: " << save_file << '\n';
```

Repeat the same pattern for the other 3 blocks, adapted to each block's specific verification method.

- [ ] **Step 4: Build and test**

```bash
make build && make test
```

Expected: all existing + new tests pass.

- [ ] **Step 5: Commit**

```bash
git add reseq/RegressionTest.cpp reseq/DataStatsTest.cpp reseq/ProbabilityEstimatesTest.cpp
git commit -m "test(6): add format detection, round-trip, and text-format mirrored tests

Cover: text auto-detect, binary auto-detect, text→binary→text round-trip,
unsupported version error, and text-format variants of all existing
Save/Load tests for DataStats and ProbabilityEstimates."
```

---

## Task 8: Format Check, Smoke Tests, and Benchmark

**Files:**
- No new files

- [ ] **Step 1: Run format check**

```bash
make format-check
```

If formatting drift, run `make format`.

- [ ] **Step 2: Run full test suite**

```bash
make build && make test
```

- [ ] **Step 3: Manual smoke tests**

```bash
# convertProfile help
build/bin/reseq convertProfile --help 2>&1 | head -5

# Convert Hs-Nova-TruSeq text profile to binary
cp test/data/Hs-Nova-TruSeq.reseq /tmp/test-convert.reseq
build/bin/reseq convertProfile -s /tmp/test-convert.reseq -o /tmp/test-binary.reseq --verbosity 0 2>/dev/null
ls -lh /tmp/test-convert.reseq /tmp/test-binary.reseq

# Verify binary loads and queries correctly
build/bin/reseq queryProfile -s /tmp/test-binary.reseq --maxReadLength --verbosity 0 2>/dev/null
# Expected: maxReadLength: 151

# Convert back to text
build/bin/reseq convertProfile -s /tmp/test-binary.reseq -o /tmp/test-roundtrip.reseq --textFormat --verbosity 0 2>/dev/null

# Verify text round-trip
build/bin/reseq queryProfile -s /tmp/test-roundtrip.reseq --maxReadLength --verbosity 0 2>/dev/null
# Expected: maxReadLength: 151
```

- [ ] **Step 4: Benchmark all 5 profiles**

```bash
for f in test/data/*.reseq; do
    name=$(basename "$f")
    # Convert to binary
    build/bin/reseq convertProfile -s "$f" -o "/tmp/${name}.bin" --verbosity 0 2>/dev/null
    text_size=$(stat -c%s "$f")
    bin_size=$(stat -c%s "/tmp/${name}.bin")
    # Time text load
    text_ms=$( { time build/bin/reseq queryProfile -s "$f" --maxReadLength --verbosity 0 2>/dev/null; } 2>&1 | grep real | awk '{print $2}')
    # Time binary load
    bin_ms=$( { time build/bin/reseq queryProfile -s "/tmp/${name}.bin" --maxReadLength --verbosity 0 2>/dev/null; } 2>&1 | grep real | awk '{print $2}')
    printf "%-25s text: %6dKB %s  binary: %6dKB %s\n" "$name" $((text_size/1024)) "$text_ms" $((bin_size/1024)) "$bin_ms"
done
```

Record results.

- [ ] **Step 5: Commit format fixes if any**

```bash
git add -A && git commit -m "style(6): apply clang-format to Phase 6 files"
```

(Only if `make format` changed files.)

---

## Verification Protocol

**After every task:**

```bash
make build && make test
```

All unit tests + regression tests must pass.

**After Task 8 (final):**

- Normal build + test
- Format check clean
- Manual smoke tests for convertProfile, --textFormat, format auto-detection
- Benchmark all 5 profiles in text vs binary format
- Verify query output equivalence across formats for all profiles
