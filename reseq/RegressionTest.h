#ifndef REGRESSIONTEST_H
#define REGRESSIONTEST_H

#include <array>
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
    std::filesystem::path reseq_bin_;
    std::filesystem::path test_dir_;
    std::filesystem::path expected_dir_;
    std::filesystem::path adapter_dir_;
    std::filesystem::path data_dir_;
    std::filesystem::path tmp_dir_;

    void SetUp() override {
        BasicTestClass::SetUp();

        // Locate the reseq binary
        reseq_bin_ = std::filesystem::path(RESEQ_BINARY_DIR) / "reseq";
        ASSERT_TRUE(std::filesystem::exists(reseq_bin_)) << "reseq binary not found at: " << reseq_bin_;

        // Locate test directory via BasicTestClass::GetTestDir
        std::string test_dir_str;
        ASSERT_TRUE(GetTestDir(test_dir_str)) << "Could not find test data directory";
        test_dir_ = std::filesystem::path(test_dir_str);

        // Derived directories
        expected_dir_ = test_dir_ / "expected";
        ASSERT_TRUE(std::filesystem::is_directory(expected_dir_)) << "Expected dir not found: " << expected_dir_;

        // GetTestDir returns path with trailing slash; use canonical parent
        auto project_root = std::filesystem::path(PROJECT_SOURCE_DIR);
        adapter_dir_ = project_root / "adapters";
        ASSERT_TRUE(std::filesystem::is_directory(adapter_dir_)) << "Adapter dir not found: " << adapter_dir_;

        data_dir_ = test_dir_ / "data";
        // data_dir_ may not exist if Zenodo data was not downloaded -- that is OK

        // Create unique temp directory
        std::string tmp_template = (std::filesystem::temp_directory_path() / "reseq_test_XXXXXX").string();
        std::vector<char> tmp_buf(tmp_template.begin(), tmp_template.end());
        tmp_buf.push_back('\0');
        char* result = mkdtemp(tmp_buf.data());
        ASSERT_NE(result, nullptr) << "Failed to create temp directory";
        tmp_dir_ = std::filesystem::path(result);
    }

    void TearDown() override {
        // Preserve temp dir on failure for debugging
        if (tmp_dir_.empty()) {
            BasicTestClass::TearDown();
            return;
        }

        const auto* test_info = ::testing::UnitTest::GetInstance()->current_test_info();
        if (test_info && test_info->result() && test_info->result()->Passed()) {
            std::filesystem::remove_all(tmp_dir_);
        } else {
            std::cerr << "Test failed -- preserving temp dir: " << tmp_dir_ << std::endl;
        }

        BasicTestClass::TearDown();
    }

    /// Run reseq with the given arguments, suppressing stderr and verbosity.
    /// Returns the process exit code.
    int RunReseq(const std::string& args) {
        std::string cmd = reseq_bin_.string() + " " + args + " --verbosity 0 2>/dev/null";
        int status = std::system(cmd.c_str());
        return WEXITSTATUS(status);
    }

    /// Run reseq and capture stdout. Stderr is suppressed.
    /// Returns the captured stdout string.
    std::string RunReseqCapture(const std::string& args) {
        std::string cmd = reseq_bin_.string() + " " + args + " --verbosity 0 2>/dev/null";
        std::string output;
        std::array<char, 4096> buffer;

        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            ADD_FAILURE() << "popen() failed for command: " << cmd;
            return "";
        }

        while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
            output += buffer.data();
        }

        pclose(pipe);
        return output;
    }

    /// Run reseq and capture stderr. Stdout is suppressed.
    /// Returns the captured stderr string.
    std::string RunReseqCaptureStderr(const std::string& args) {
        std::string cmd = reseq_bin_.string() + " " + args + " 2>&1 1>/dev/null";
        std::string output;
        std::array<char, 4096> buffer;

        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            ADD_FAILURE() << "popen() failed for command: " << cmd;
            return "";
        }

        while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
            output += buffer.data();
        }

        pclose(pipe);
        return output;
    }

    /// Run reseq without suppressing stderr or forcing verbosity.
    /// Returns only the exit code.
    int RunReseqExitOnly(const std::string& args) {
        std::string cmd = reseq_bin_.string() + " " + args + " >/dev/null 2>/dev/null";
        int status = std::system(cmd.c_str());
        return WEXITSTATUS(status);
    }

    /// Generate an E. coli profile in tmp_dir_ using the 4-pair BAM.
    /// Returns the path to the generated .reseq file.
    std::filesystem::path GenerateEcoliProfile() {
        auto ref = test_dir_ / "ecoli-GCF_000005845.2_ASM584v2_genomic.fa";
        auto bam = test_dir_ / "ecoli-SRR490124-4pairs.bam";
        auto adapter_fa = adapter_dir_ / "TruSeq_single.fa";
        auto adapter_mat = adapter_dir_ / "TruSeq_single.mat";
        auto output = tmp_dir_ / "ecoli-4pairs.reseq";

        std::string args = "illuminaPE"
                           " -r " +
                           ref.string() + " -b " + bam.string() + " --adapterFile " + adapter_fa.string() +
                           " --adapterMatrix " + adapter_mat.string() +
                           " --statsOnly --noBias"
                           " -S " +
                           output.string() + " -j 1";

        int rc = RunReseq(args);
        EXPECT_EQ(0, rc) << "illuminaPE profile generation failed with exit code " << rc;

        return output;
    }

    /// Return path to Zenodo profile if it exists, otherwise empty path.
    std::filesystem::path ZenodoProfile() {
        auto path = data_dir_ / "Hs-Nova-TruSeq.reseq";
        if (std::filesystem::exists(path)) {
            return path;
        }
        return {};
    }

    /// Compare two text files line by line with descriptive error messages.
    void ExpectTextFilesEqual(const std::filesystem::path& expected, const std::filesystem::path& actual) {
        std::ifstream exp_stream(expected);
        std::ifstream act_stream(actual);
        ASSERT_TRUE(exp_stream.is_open()) << "Cannot open expected file: " << expected;
        ASSERT_TRUE(act_stream.is_open()) << "Cannot open actual file: " << actual;

        std::string exp_line, act_line;
        int line_num = 0;
        while (std::getline(exp_stream, exp_line)) {
            ++line_num;
            ASSERT_TRUE(std::getline(act_stream, act_line))
                << "Actual file has fewer lines than expected at line " << line_num << "\n"
                << "  Expected file: " << expected << "\n"
                << "  Actual file:   " << actual;
            EXPECT_EQ(exp_line, act_line) << "Line " << line_num << " differs:\n"
                                          << "  Expected file: " << expected << "\n"
                                          << "  Actual file:   " << actual;
        }

        if (std::getline(act_stream, act_line)) {
            ++line_num;
            ADD_FAILURE() << "Actual file has more lines than expected (extra content at line " << line_num << ")\n"
                          << "  Expected file: " << expected << "\n"
                          << "  Actual file:   " << actual;
        }
    }

    /// Binary comparison of two files (size check first, then byte-by-byte).
    void ExpectFilesEqual(const std::filesystem::path& expected, const std::filesystem::path& actual) {
        ASSERT_TRUE(std::filesystem::exists(expected)) << "Expected file does not exist: " << expected;
        ASSERT_TRUE(std::filesystem::exists(actual)) << "Actual file does not exist: " << actual;

        auto exp_size = std::filesystem::file_size(expected);
        auto act_size = std::filesystem::file_size(actual);
        ASSERT_EQ(exp_size, act_size) << "File sizes differ:\n"
                                      << "  Expected: " << expected << " (" << exp_size << " bytes)\n"
                                      << "  Actual:   " << actual << " (" << act_size << " bytes)";

        std::ifstream exp_stream(expected, std::ios::binary);
        std::ifstream act_stream(actual, std::ios::binary);
        ASSERT_TRUE(exp_stream.is_open()) << "Cannot open expected file: " << expected;
        ASSERT_TRUE(act_stream.is_open()) << "Cannot open actual file: " << actual;

        std::array<char, 8192> exp_buf, act_buf;
        while (exp_stream) {
            exp_stream.read(exp_buf.data(), exp_buf.size());
            act_stream.read(act_buf.data(), act_buf.size());
            auto bytes_read = exp_stream.gcount();
            ASSERT_EQ(bytes_read, act_stream.gcount());
            if (std::memcmp(exp_buf.data(), act_buf.data(), static_cast<size_t>(bytes_read)) != 0) {
                FAIL() << "Binary content differs between:\n"
                       << "  Expected: " << expected << "\n"
                       << "  Actual:   " << actual;
            }
        }
    }

  public:
    static void Register();
};

} // namespace reseq

#endif // REGRESSIONTEST_H
