# Phase 2: Golden File Regression Tests — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Capture current CLI output as committed golden files and write GoogleTest fixtures that compare against them, creating a safety net for all subsequent refactoring phases. Uses a real Illumina error profile hosted on Zenodo for full pipeline regression tests.

**Architecture:** Subprocess-based regression tests in the existing `reseq_test` binary invoke the `reseq` CLI with fixed arguments (seed=42, threads=1), write output to per-test temp directories, and diff against committed golden files in `test/expected/`. A `RegressionTest` fixture provides helpers for running the binary, comparing files, and managing temp dirs. A real 98MB `Hs-Nova-TruSeq.reseq` profile is hosted on Zenodo (DOI: 10.5281/zenodo.19383555) and downloaded by a helper script to `test/data/` (gitignored). Tests that need the profile skip gracefully if it hasn't been downloaded. CI always downloads it.

**Tech Stack:** GoogleTest, `std::filesystem`, `popen()` for subprocess execution, Zenodo for large test data, curl for download

**Spec:** [docs/superpowers/specs/2026-04-01-comprehensive-refactoring-design.md](../specs/2026-04-01-comprehensive-refactoring-design.md) — Phase 2

---

## Zenodo Test Data

| File | Size | MD5 | URL |
|------|------|-----|-----|
| `Hs-Nova-TruSeq.reseq` | 98MB | `c374ef7198effa8ad6ed7fc9e6d9431e` | `https://zenodo.org/api/records/19383555/files/Hs-Nova-TruSeq.reseq/content` |
| `Hs-Nova-TruSeq.reseq.ipf` | 53MB | `fbea201bcae34ea025d6408e5e9e91d5` | `https://zenodo.org/api/records/19383555/files/Hs-Nova-TruSeq.reseq.ipf/content` |

---

## File Map

| Action | File | Purpose |
|--------|------|---------|
| Create | `test/download_test_data.sh` | Downloads Zenodo profiles to `test/data/`, verifies MD5 |
| Create | `reseq/RegressionTest.h` | Test fixture: subprocess helpers, file comparison, temp dir, profile discovery |
| Create | `reseq/RegressionTest.cpp` | All regression test cases |
| Create | `test/expected/replaceN.fa` | Golden file: `replaceN` output (seed=42) |
| Create | `test/expected/queryProfile_maxReadLength.txt` | Golden file: queryProfile stdout (small BAM profile) |
| Create | `test/expected/queryProfile_maxLenDeletion.txt` | Golden file: queryProfile stdout (small BAM profile) |
| Create | `test/expected/queryProfile_fragLenBias.tsv` | Golden file: queryProfile bias (small BAM profile, 2001 lines) |
| Create | `test/expected/queryProfile_real_maxReadLength.txt` | Golden file: queryProfile stdout (Zenodo Hs-Nova-TruSeq profile) |
| Create | `test/expected/queryProfile_real_maxLenDeletion.txt` | Golden file: queryProfile stdout (Zenodo Hs-Nova-TruSeq profile) |
| Create | `test/expected/queryProfile_real_fragLenBias.tsv` | Golden file: queryProfile bias (Zenodo profile, 2001 lines) |
| Modify | `reseq/test_main.cpp` | Add `#include "RegressionTest.h"` and `Register()` call |
| Modify | `reseq/CMakeLists.txt` | Add `RegressionTest.cpp` to `reseq_test` sources |
| Modify | `.gitignore` | Add `test/data/` for downloaded Zenodo files |
| Modify | `.github/workflows/ci.yml` | Add test data download step before test execution |

---

### Task 1: Create Test Data Download Script

**Files:**
- Create: `test/download_test_data.sh`
- Modify: `.gitignore`

- [ ] **Step 1: Add test/data/ to .gitignore**

In `.gitignore`, add:

```
test/data/
```

- [ ] **Step 2: Create the download script**

Create `test/download_test_data.sh`:

```bash
#!/usr/bin/env bash
# Downloads large test data files from Zenodo for regression testing.
# Files are cached in test/data/ and verified by MD5 checksum.
# Usage: ./test/download_test_data.sh
#
# Zenodo record: https://zenodo.org/records/19383555
# DOI: 10.5281/zenodo.19383555

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
ZENODO_RECORD="https://zenodo.org/api/records/19383555/files"

mkdir -p "$DATA_DIR"

download_and_verify() {
    local filename="$1"
    local expected_md5="$2"
    local dest="$DATA_DIR/$filename"

    if [ -f "$dest" ]; then
        # Verify existing file
        local actual_md5
        actual_md5=$(md5sum "$dest" | awk '{print $1}')
        if [ "$actual_md5" = "$expected_md5" ]; then
            echo "OK: $filename already exists and checksum matches"
            return 0
        else
            echo "WARN: $filename exists but checksum mismatch, re-downloading"
            rm -f "$dest"
        fi
    fi

    echo "Downloading $filename..."
    curl -L -o "$dest" "$ZENODO_RECORD/$filename/content"

    local actual_md5
    actual_md5=$(md5sum "$dest" | awk '{print $1}')
    if [ "$actual_md5" != "$expected_md5" ]; then
        echo "ERROR: MD5 mismatch for $filename"
        echo "  Expected: $expected_md5"
        echo "  Actual:   $actual_md5"
        rm -f "$dest"
        return 1
    fi
    echo "OK: $filename downloaded and verified"
}

download_and_verify "Hs-Nova-TruSeq.reseq" "c374ef7198effa8ad6ed7fc9e6d9431e"
download_and_verify "Hs-Nova-TruSeq.reseq.ipf" "fbea201bcae34ea025d6408e5e9e91d5"

echo ""
echo "All test data files ready in $DATA_DIR"
```

- [ ] **Step 3: Make it executable**

```bash
chmod +x test/download_test_data.sh
```

- [ ] **Step 4: Test the download script**

```bash
./test/download_test_data.sh
ls -lh test/data/
```

Expected:
```
OK: Hs-Nova-TruSeq.reseq downloaded and verified
OK: Hs-Nova-TruSeq.reseq.ipf downloaded and verified

All test data files ready in .../test/data
```

- [ ] **Step 5: Run it again to verify caching works**

```bash
./test/download_test_data.sh
```

Expected:
```
OK: Hs-Nova-TruSeq.reseq already exists and checksum matches
OK: Hs-Nova-TruSeq.reseq.ipf already exists and checksum matches
```

- [ ] **Step 6: Commit**

```bash
git add test/download_test_data.sh .gitignore
git commit -m "test: add Zenodo test data download script

Downloads Hs-Nova-TruSeq.reseq profile (98MB) and .ipf (53MB) from
Zenodo (DOI: 10.5281/zenodo.19383555) with MD5 verification and caching.
Files stored in test/data/ (gitignored)."
```

---

### Task 2: Generate and Commit Golden Files

**Files:**
- Create: `test/expected/replaceN.fa`
- Create: `test/expected/queryProfile_maxReadLength.txt`
- Create: `test/expected/queryProfile_maxLenDeletion.txt`
- Create: `test/expected/queryProfile_fragLenBias.tsv`
- Create: `test/expected/queryProfile_real_maxReadLength.txt`
- Create: `test/expected/queryProfile_real_maxLenDeletion.txt`
- Create: `test/expected/queryProfile_real_fragLenBias.tsv`

- [ ] **Step 1: Build the current tree**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
```

- [ ] **Step 2: Create the expected output directory**

```bash
mkdir -p test/expected
```

- [ ] **Step 3: Generate replaceN golden file**

```bash
build/bin/reseq replaceN \
  -r test/reference-test.fa \
  -R test/expected/replaceN.fa \
  --seed 42 --verbosity 0
```

Verify:
```bash
wc -l test/expected/replaceN.fa
# Expected: 18 lines
```

- [ ] **Step 4: Generate queryProfile golden files from small BAM profile**

Generate a profile at runtime from the small test BAM (not committed — 79MB):

```bash
TMPD=$(mktemp -d)
build/bin/reseq illuminaPE \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -b test/ecoli-SRR490124-4pairs.bam \
  --adapterFile adapters/TruSeq_single.fa \
  --adapterMatrix adapters/TruSeq_single.mat \
  --statsOnly --noBias \
  -S "$TMPD/ecoli.reseq" \
  -j 1 --verbosity 1 2>&1 | tail -3

build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s "$TMPD/ecoli.reseq" \
  --maxReadLength --verbosity 0 2>/dev/null \
  > test/expected/queryProfile_maxReadLength.txt

build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s "$TMPD/ecoli.reseq" \
  --maxLenDeletion --verbosity 0 2>/dev/null \
  > test/expected/queryProfile_maxLenDeletion.txt

build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s "$TMPD/ecoli.reseq" \
  --fragLenBias test/expected/queryProfile_fragLenBias.tsv \
  --verbosity 0 2>/dev/null

rm -rf "$TMPD"
```

Verify:
```bash
cat test/expected/queryProfile_maxReadLength.txt
# Expected: "maxReadLength: 100"
cat test/expected/queryProfile_maxLenDeletion.txt
# Expected: "maxLenDeletion: 0"
wc -l test/expected/queryProfile_fragLenBias.tsv
# Expected: 2001 lines
```

- [ ] **Step 5: Generate queryProfile golden files from Zenodo profile**

Ensure Zenodo data is downloaded:

```bash
./test/download_test_data.sh
```

Then generate golden files using the real Hs-Nova-TruSeq profile:

```bash
build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s test/data/Hs-Nova-TruSeq.reseq \
  --maxReadLength --verbosity 0 2>/dev/null \
  > test/expected/queryProfile_real_maxReadLength.txt

build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s test/data/Hs-Nova-TruSeq.reseq \
  --maxLenDeletion --verbosity 0 2>/dev/null \
  > test/expected/queryProfile_real_maxLenDeletion.txt

build/bin/reseq queryProfile \
  -r test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
  -s test/data/Hs-Nova-TruSeq.reseq \
  --fragLenBias test/expected/queryProfile_real_fragLenBias.tsv \
  --verbosity 0 2>/dev/null
```

Verify:
```bash
cat test/expected/queryProfile_real_maxReadLength.txt
# Expected: "maxReadLength: 151"
cat test/expected/queryProfile_real_maxLenDeletion.txt
# Expected: "maxLenDeletion: 29"
wc -l test/expected/queryProfile_real_fragLenBias.tsv
# Expected: 2001 lines
```

- [ ] **Step 6: Commit golden files**

```bash
git add test/expected/
git commit -m "test: add golden files for regression tests

- replaceN output (seed=42, reference-test.fa)
- queryProfile outputs from small BAM profile (maxReadLength=100,
  maxLenDeletion=0, fragLenBias 2001 lines)
- queryProfile outputs from real Hs-Nova-TruSeq Zenodo profile
  (maxReadLength=151, maxLenDeletion=29, fragLenBias 2001 lines)"
```

---

### Task 3: Create RegressionTest Header

**Files:**
- Create: `reseq/RegressionTest.h`

- [ ] **Step 1: Create RegressionTest.h**

```cpp
#ifndef REGRESSIONTEST_H
#define REGRESSIONTEST_H

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "gtest/gtest.h"

#include "BasicTestClass.hpp"
#include "CMakeConfig.h"

namespace reseq {

class RegressionTest : public BasicTestClass {
  protected:
    std::filesystem::path tmp_dir_;
    std::filesystem::path reseq_bin_;
    std::filesystem::path test_dir_;
    std::filesystem::path expected_dir_;
    std::filesystem::path adapter_dir_;
    std::filesystem::path data_dir_;

    void SetUp() override {
        BasicTestClass::SetUp();

        reseq_bin_ = std::filesystem::path(RESEQ_BINARY_DIR) / "reseq";
        ASSERT_TRUE(std::filesystem::exists(reseq_bin_))
            << "reseq binary not found at: " << reseq_bin_;

        std::string td;
        ASSERT_TRUE(GetTestDir(td));
        test_dir_ = std::filesystem::path(td);
        expected_dir_ = test_dir_ / "expected";
        ASSERT_TRUE(std::filesystem::exists(expected_dir_))
            << "Expected output directory not found: " << expected_dir_;

        adapter_dir_ = test_dir_.parent_path() / "adapters";
        data_dir_ = test_dir_ / "data";

        auto tmpl = std::filesystem::temp_directory_path() / "reseq_regtest_XXXXXX";
        std::string tmpl_str = tmpl.string();
        char* result = mkdtemp(tmpl_str.data());
        ASSERT_NE(result, nullptr) << "Failed to create temp directory";
        tmp_dir_ = std::filesystem::path(tmpl_str);
    }

    void TearDown() override {
        if (!::testing::Test::HasFailure() && std::filesystem::exists(tmp_dir_)) {
            std::filesystem::remove_all(tmp_dir_);
        }
        BasicTestClass::TearDown();
    }

    int RunReseq(const std::string& args) const {
        std::string cmd = reseq_bin_.string() + " " + args + " --verbosity 0 2>/dev/null";
        return std::system(cmd.c_str());
    }

    std::string RunReseqCapture(const std::string& args) const {
        std::string cmd = reseq_bin_.string() + " " + args + " --verbosity 0 2>/dev/null";
        std::string output;
        FILE* pipe = popen(cmd.c_str(), "r");
        if (pipe) {
            char buffer[256];
            while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                output += buffer;
            }
            pclose(pipe);
        }
        return output;
    }

    // Generate a .reseq profile from the small ecoli BAM with explicit adapters.
    std::filesystem::path GenerateEcoliProfile() {
        auto profile = tmp_dir_ / "ecoli.reseq";
        int rc = RunReseq(
            "illuminaPE"
            " -r " + (test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa").string() +
            " -b " + (test_dir_ / "ecoli-SRR490124-4pairs.bam").string() +
            " --adapterFile " + (adapter_dir_ / "TruSeq_single.fa").string() +
            " --adapterMatrix " + (adapter_dir_ / "TruSeq_single.mat").string() +
            " --statsOnly --noBias"
            " -S " + profile.string() +
            " -j 1");
        EXPECT_EQ(rc, 0) << "Failed to generate .reseq profile";
        return profile;
    }

    // Returns the path to the Zenodo Hs-Nova-TruSeq profile, or empty if not downloaded.
    std::filesystem::path ZenodoProfile() const {
        auto profile = data_dir_ / "Hs-Nova-TruSeq.reseq";
        if (std::filesystem::exists(profile)) {
            return profile;
        }
        return {};
    }

    static void ExpectTextFilesEqual(const std::filesystem::path& expected,
                                     const std::filesystem::path& actual) {
        ASSERT_TRUE(std::filesystem::exists(expected))
            << "Expected file not found: " << expected;
        ASSERT_TRUE(std::filesystem::exists(actual))
            << "Actual file not found: " << actual;

        std::ifstream exp_stream(expected);
        std::ifstream act_stream(actual);
        std::string exp_line, act_line;
        int line_num = 0;

        while (std::getline(exp_stream, exp_line)) {
            line_num++;
            ASSERT_TRUE(std::getline(act_stream, act_line))
                << "Actual file has fewer lines at line " << line_num
                << "\n  expected: " << expected << "\n  actual: " << actual;
            EXPECT_EQ(exp_line, act_line)
                << "Files differ at line " << line_num
                << "\n  expected: " << expected << "\n  actual: " << actual;
        }
        EXPECT_FALSE(std::getline(act_stream, act_line))
            << "Actual file has more lines than expected"
            << "\n  expected: " << expected << "\n  actual: " << actual;
    }

    static void ExpectFilesEqual(const std::filesystem::path& expected,
                                 const std::filesystem::path& actual) {
        ASSERT_TRUE(std::filesystem::exists(expected))
            << "Expected file not found: " << expected;
        ASSERT_TRUE(std::filesystem::exists(actual))
            << "Actual file not found: " << actual;

        auto exp_size = std::filesystem::file_size(expected);
        auto act_size = std::filesystem::file_size(actual);
        ASSERT_EQ(exp_size, act_size)
            << "File sizes differ: expected=" << exp_size << " actual=" << act_size
            << "\n  expected: " << expected << "\n  actual: " << actual;

        std::ifstream exp_stream(expected, std::ios::binary);
        std::ifstream act_stream(actual, std::ios::binary);
        std::string exp_content((std::istreambuf_iterator<char>(exp_stream)),
                                std::istreambuf_iterator<char>());
        std::string act_content((std::istreambuf_iterator<char>(act_stream)),
                                std::istreambuf_iterator<char>());
        EXPECT_EQ(exp_content, act_content)
            << "File contents differ:\n  expected: " << expected
            << "\n  actual: " << actual;
    }

  public:
    static void Register() {}
};

} // namespace reseq

#endif // REGRESSIONTEST_H
```

- [ ] **Step 2: Commit**

```bash
git add reseq/RegressionTest.h
git commit -m "test: add RegressionTest base class with subprocess and file helpers"
```

---

### Task 4: Wire RegressionTest into Build and Add All Test Cases

**Files:**
- Create: `reseq/RegressionTest.cpp`
- Modify: `reseq/test_main.cpp`
- Modify: `reseq/CMakeLists.txt`

- [ ] **Step 1: Create RegressionTest.cpp with all test cases**

```cpp
#include "RegressionTest.h"

namespace reseq {

// ── replaceN ────────────────────────────────────────────────────────────────

TEST_F(RegressionTest, ReplaceN) {
    auto output = tmp_dir_ / "replaced.fa";
    auto expected = expected_dir_ / "replaceN.fa";

    int rc = RunReseq(
        "replaceN"
        " -r " + (test_dir_ / "reference-test.fa").string() +
        " -R " + output.string() +
        " --seed 42");
    ASSERT_EQ(rc, 0) << "reseq replaceN exited with non-zero status";
    ExpectTextFilesEqual(expected, output);
}

// ── queryProfile from small BAM profile (generated at test time) ────────────

TEST_F(RegressionTest, QueryProfileMaxReadLength) {
    auto profile = GenerateEcoliProfile();
    ASSERT_TRUE(std::filesystem::exists(profile)) << "Profile generation failed";

    auto expected = expected_dir_ / "queryProfile_maxReadLength.txt";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";

    std::string output = RunReseqCapture(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --maxReadLength");

    std::ifstream exp_file(expected);
    std::string exp_content((std::istreambuf_iterator<char>(exp_file)),
                            std::istreambuf_iterator<char>());
    EXPECT_EQ(output, exp_content)
        << "queryProfile --maxReadLength output differs from golden file";
}

TEST_F(RegressionTest, QueryProfileMaxLenDeletion) {
    auto profile = GenerateEcoliProfile();
    ASSERT_TRUE(std::filesystem::exists(profile)) << "Profile generation failed";

    auto expected = expected_dir_ / "queryProfile_maxLenDeletion.txt";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";

    std::string output = RunReseqCapture(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --maxLenDeletion");

    std::ifstream exp_file(expected);
    std::string exp_content((std::istreambuf_iterator<char>(exp_file)),
                            std::istreambuf_iterator<char>());
    EXPECT_EQ(output, exp_content)
        << "queryProfile --maxLenDeletion output differs from golden file";
}

TEST_F(RegressionTest, QueryProfileFragLenBias) {
    auto profile = GenerateEcoliProfile();
    ASSERT_TRUE(std::filesystem::exists(profile)) << "Profile generation failed";

    auto expected = expected_dir_ / "queryProfile_fragLenBias.tsv";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";
    auto output = tmp_dir_ / "fragLenBias.tsv";

    int rc = RunReseq(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --fragLenBias " + output.string());
    ASSERT_EQ(rc, 0) << "reseq queryProfile --fragLenBias exited with non-zero status";
    ExpectTextFilesEqual(expected, output);
}

// ── queryProfile from real Zenodo Hs-Nova-TruSeq profile ────────────────────

TEST_F(RegressionTest, RealProfileMaxReadLength) {
    auto profile = ZenodoProfile();
    if (profile.empty()) {
        GTEST_SKIP() << "Zenodo test data not downloaded (run test/download_test_data.sh)";
    }

    auto expected = expected_dir_ / "queryProfile_real_maxReadLength.txt";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";

    std::string output = RunReseqCapture(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --maxReadLength");

    std::ifstream exp_file(expected);
    std::string exp_content((std::istreambuf_iterator<char>(exp_file)),
                            std::istreambuf_iterator<char>());
    EXPECT_EQ(output, exp_content)
        << "queryProfile --maxReadLength (real profile) differs from golden file";
}

TEST_F(RegressionTest, RealProfileMaxLenDeletion) {
    auto profile = ZenodoProfile();
    if (profile.empty()) {
        GTEST_SKIP() << "Zenodo test data not downloaded (run test/download_test_data.sh)";
    }

    auto expected = expected_dir_ / "queryProfile_real_maxLenDeletion.txt";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";

    std::string output = RunReseqCapture(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --maxLenDeletion");

    std::ifstream exp_file(expected);
    std::string exp_content((std::istreambuf_iterator<char>(exp_file)),
                            std::istreambuf_iterator<char>());
    EXPECT_EQ(output, exp_content)
        << "queryProfile --maxLenDeletion (real profile) differs from golden file";
}

TEST_F(RegressionTest, RealProfileFragLenBias) {
    auto profile = ZenodoProfile();
    if (profile.empty()) {
        GTEST_SKIP() << "Zenodo test data not downloaded (run test/download_test_data.sh)";
    }

    auto expected = expected_dir_ / "queryProfile_real_fragLenBias.tsv";
    auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";
    auto output = tmp_dir_ / "fragLenBias.tsv";

    int rc = RunReseq(
        "queryProfile"
        " -r " + ref.string() +
        " -s " + profile.string() +
        " --fragLenBias " + output.string());
    ASSERT_EQ(rc, 0) << "queryProfile --fragLenBias (real profile) exited with non-zero status";
    ExpectTextFilesEqual(expected, output);
}

// ── Error handling ──────────────────────────────────────────────────────────

TEST_F(RegressionTest, ErrorBadCommand) {
    int rc = RunReseq("nonexistent_command");
    EXPECT_NE(rc, 0) << "reseq should exit non-zero on unknown command";
}

TEST_F(RegressionTest, ErrorMissingRef) {
    int rc = RunReseq(
        "replaceN"
        " -r /nonexistent/file.fa"
        " -R " + (tmp_dir_ / "out.fa").string() +
        " --seed 42");
    EXPECT_NE(rc, 0) << "reseq replaceN should fail with missing reference";
}

TEST_F(RegressionTest, VersionOutput) {
    std::string output = RunReseqCapture("--version");
    EXPECT_FALSE(output.empty()) << "reseq --version should produce output";
    EXPECT_NE(output.find("ReSeq"), std::string::npos)
        << "reseq --version should contain 'ReSeq'";
}

} // namespace reseq
```

- [ ] **Step 2: Add RegressionTest.cpp to CMakeLists.txt**

In `reseq/CMakeLists.txt`, add `RegressionTest.cpp` to the `reseq_test` source list:

```cmake
  add_executable(reseq_test
    test_main.cpp
    AdapterStatsTest.cpp CoverageStatsTest.cpp DataStatsTest.cpp
    ErrorStatsTest.cpp FragmentDistributionStatsTest.cpp
    FragmentDuplicationStatsTest.cpp ProbabilityEstimatesTest.cpp
    QualityStatsTest.cpp ReferenceTest.cpp SeqQualityStatsTest.cpp
    SimulatorTest.cpp SurroundingTest.cpp TileStatsTest.cpp
    utilitiesTest.cpp VectTest.cpp
    RegressionTest.cpp
  )
```

- [ ] **Step 3: Add include and Register to test_main.cpp**

In `reseq/test_main.cpp`, add after the existing test headers:

```cpp
#include "RegressionTest.h"
```

And add before `::testing::InitGoogleTest`:

```cpp
    reseq::RegressionTest::Register();
```

- [ ] **Step 4: Build and run all tests**

```bash
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure && cd ..
```

Expected: All tests pass. Zenodo profile tests run if downloaded, skip if not.

- [ ] **Step 5: Run only regression tests**

```bash
build/bin/reseq_test --gtest_filter="RegressionTest.*"
```

Expected: 10 regression tests (4 basic + 3 real profile + 3 error). Real profile tests pass if Zenodo data is present, skip otherwise.

- [ ] **Step 6: Commit**

```bash
git add reseq/RegressionTest.cpp reseq/CMakeLists.txt reseq/test_main.cpp
git commit -m "test: add regression tests for replaceN, queryProfile, and error handling

- ReplaceN: golden file comparison (seed=42)
- QueryProfile (small BAM): generates profile at test time, compares
  maxReadLength/maxLenDeletion/fragLenBias against golden files
- QueryProfile (real Zenodo profile): uses Hs-Nova-TruSeq.reseq from
  Zenodo (DOI: 10.5281/zenodo.19383555), skips if not downloaded
- Error handling: bad command, missing reference, version output"
```

---

### Task 5: Add Zenodo Data Download to CI

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Add download step to build-and-test job**

In `.github/workflows/ci.yml`, add a step after "Install dependencies" and before "Configure" in the `build-and-test` job:

```yaml
      - name: Download test data
        run: ./test/download_test_data.sh
```

Also add it to the `sanitizers` and `coverage` jobs (same position — after "Install dependencies", before "Configure"):

```yaml
      - name: Download test data
        run: ./test/download_test_data.sh
```

- [ ] **Step 2: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: download Zenodo test data before running tests

Ensures regression tests using the real Hs-Nova-TruSeq profile
run in CI. Data is cached via MD5 verification."
```

---

### Task 6: Add Makefile test-data Target

**Files:**
- Modify: `Makefile`

- [ ] **Step 1: Add test-data target**

In `Makefile`, add after the `test` target:

```makefile
test-data:
	./test/download_test_data.sh
```

Update `test` to depend on `test-data`:

```makefile
test: build test-data
	cd $(BUILD_DIR) && ctest --output-on-failure
```

Add `test-data` to `.PHONY`:

```makefile
.PHONY: all configure build test test-data coverage format format-check lint clean install changelog help
```

Add to `help`:

```makefile
	@echo "  test-data    Download large test data from Zenodo"
```

- [ ] **Step 2: Test**

```bash
make test
```

Expected: Downloads data (or uses cache), builds, runs tests — all pass.

- [ ] **Step 3: Commit**

```bash
git add Makefile
git commit -m "build: add test-data Makefile target, make test depend on it"
```

---

### Task 7: Final Verification

**Files:**
- None modified — verification only

- [ ] **Step 1: Clean build from scratch**

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure && cd ..
```

Expected: All tests pass.

- [ ] **Step 2: List regression tests**

```bash
build/bin/reseq_test --gtest_filter="RegressionTest.*" --gtest_list_tests
```

Expected:
```
RegressionTest.
  ReplaceN
  QueryProfileMaxReadLength
  QueryProfileMaxLenDeletion
  QueryProfileFragLenBias
  RealProfileMaxReadLength
  RealProfileMaxLenDeletion
  RealProfileFragLenBias
  ErrorBadCommand
  ErrorMissingRef
  VersionOutput
```

- [ ] **Step 3: Verify format check passes**

```bash
make format-check
```

- [ ] **Step 4: Verify golden files are committed**

```bash
git ls-files test/expected/
```

Expected:
```
test/expected/queryProfile_fragLenBias.tsv
test/expected/queryProfile_maxLenDeletion.txt
test/expected/queryProfile_maxReadLength.txt
test/expected/queryProfile_real_fragLenBias.tsv
test/expected/queryProfile_real_maxLenDeletion.txt
test/expected/queryProfile_real_maxReadLength.txt
test/expected/replaceN.fa
```

- [ ] **Step 5: Verify Zenodo data is NOT tracked**

```bash
git ls-files test/data/
# Expected: empty (gitignored)
```

---

## Summary

| Aspect | Before | After |
|--------|--------|-------|
| Golden files | None | 7 files in `test/expected/` |
| Regression tests | None | 10 test cases in `RegressionTest` fixture |
| Test types | Unit tests only | Unit tests + subprocess regression tests |
| Test data hosting | N/A | Zenodo DOI: 10.5281/zenodo.19383555 |
| Profile handling | N/A | Small: generated at test time. Real: downloaded from Zenodo |
| Safety net | Internal logic only | CLI output verified against committed baselines |
| CI integration | N/A | Zenodo data auto-downloaded in all CI jobs |

### Test Matrix

| Test | Profile Source | Skippable? | Golden File |
|------|---------------|------------|-------------|
| ReplaceN | None needed | No | `replaceN.fa` |
| QueryProfileMaxReadLength | Generated from ecoli-4pairs BAM | No | `queryProfile_maxReadLength.txt` |
| QueryProfileMaxLenDeletion | Generated from ecoli-4pairs BAM | No | `queryProfile_maxLenDeletion.txt` |
| QueryProfileFragLenBias | Generated from ecoli-4pairs BAM | No | `queryProfile_fragLenBias.tsv` |
| RealProfileMaxReadLength | Zenodo Hs-Nova-TruSeq | Yes (skip if not downloaded) | `queryProfile_real_maxReadLength.txt` |
| RealProfileMaxLenDeletion | Zenodo Hs-Nova-TruSeq | Yes (skip if not downloaded) | `queryProfile_real_maxLenDeletion.txt` |
| RealProfileFragLenBias | Zenodo Hs-Nova-TruSeq | Yes (skip if not downloaded) | `queryProfile_real_fragLenBias.tsv` |
| ErrorBadCommand | None | No | N/A |
| ErrorMissingRef | None | No | N/A |
| VersionOutput | None | No | N/A |
