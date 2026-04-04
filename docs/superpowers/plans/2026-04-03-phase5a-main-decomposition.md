# Phase 5a: main.cpp Decomposition — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decompose `reseq/main.cpp` (1,167 lines) into a ~100-line dispatcher + focused CLI command files under `reseq/cli/`, preserving exact observable behavior.

**Architecture:** Behavior-preserving extraction with adapter wrappers. Helper functions moved with exact signatures. Each command wrapped in a `Run*()` free function. File-local helpers use unnamed namespaces. Self-contained headers. New regression tests lock down CLI behavior before any code moves.

**Tech Stack:** C++20, Boost.program_options, GoogleTest, CMake `target_sources`.

**Prerequisite:** Phase 4 complete. `make build && make test` passes.

**Design spec:** `docs/superpowers/specs/2026-04-03-phase5a-main-decomposition-design.md` (rev 4)

---

## Task 1: Add CLI Behavior Regression Tests (Prerequisite)

**Files:**
- Modify: `reseq/RegressionTest.h` (add `RunReseqExitOnly` helper)
- Modify: `reseq/RegressionTest.cpp` (add new tests)

These tests lock down exact CLI behavior BEFORE any extraction begins.

- [ ] **Step 1: Add RunReseqExitOnly helper to RegressionTest.h**

Add after the existing `RunReseqCaptureStderr` method (around line 125):

```cpp
/// Run reseq without suppressing stderr or forcing verbosity.
/// Returns only the exit code.
int RunReseqExitOnly(const std::string& args) {
    std::string cmd = reseq_bin_.string() + " " + args + " >/dev/null 2>/dev/null";
    int status = std::system(cmd.c_str());
    return WEXITSTATUS(status);
}
```

- [ ] **Step 2: Add CLI behavior tests to RegressionTest.cpp**

Add after the existing `VersionOutput` test (before the closing `} // namespace reseq`):

```cpp
// --- CLI behavior tests (Phase 5a prerequisite) ---

TEST_F(RegressionTest, BareReseqExitCode) {
    int rc = RunReseqExitOnly("");
    EXPECT_EQ(0, rc) << "Bare reseq (no args) should exit 0";
}

TEST_F(RegressionTest, BareReseqOutput) {
    std::string stderr_out = RunReseqCaptureStderr("");
    EXPECT_NE(std::string::npos, stderr_out.find("reseq <command>"))
        << "Bare reseq should print usage containing 'reseq <command>', got:\n" << stderr_out;
}

TEST_F(RegressionTest, HelpSameAsBare) {
    std::string bare = RunReseqCaptureStderr("");
    std::string help = RunReseqCaptureStderr("--help");
    EXPECT_EQ(bare, help) << "reseq --help should produce identical output to bare reseq";
}

TEST_F(RegressionTest, UnknownCommandExitCode) {
    int rc = RunReseqExitOnly("nonsenseCommand123");
    EXPECT_EQ(1, rc) << "Unknown command should exit 1";
}

TEST_F(RegressionTest, UnknownCommandStderr) {
    std::string stderr_out = RunReseqCaptureStderr("nonsenseCommand123");
    EXPECT_NE(std::string::npos, stderr_out.find("Unrecognized command: 'nonsenseCommand123'"))
        << "Should contain exact error message, got:\n" << stderr_out;
}

TEST_F(RegressionTest, CommandHelpExitCode) {
    int rc = RunReseqExitOnly("replaceN --help");
    EXPECT_EQ(0, rc) << "replaceN --help should exit 0";
}

TEST_F(RegressionTest, VersionExitCode) {
    int rc = RunReseqExitOnly("--version");
    EXPECT_EQ(0, rc) << "reseq --version should exit 0";
}

TEST_F(RegressionTest, SeqToIlluminaHelpExitCode) {
    int rc = RunReseqExitOnly("seqToIllumina --help");
    EXPECT_EQ(0, rc) << "seqToIllumina --help should exit 0";
}

TEST_F(RegressionTest, IlluminaPEHelpExitCode) {
    int rc = RunReseqExitOnly("illuminaPE --help");
    EXPECT_EQ(0, rc) << "illuminaPE --help should exit 0";
}
```

- [ ] **Step 3: Build and run tests**

```bash
make build && make test
```

Expected: all existing + new tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/RegressionTest.h reseq/RegressionTest.cpp
git commit -m "test(5a): add CLI behavior regression tests before extraction

Lock down exit codes and stderr output for bare reseq, --help, --version,
unknown command, and per-command --help. Prerequisite for main.cpp
decomposition — these catch behavioral regressions during extraction."
```

---

## Task 2: Extract Shared Helpers into cli_common

**Files:**
- Create: `reseq/cli/cli_common.h`
- Create: `reseq/cli/cli_common.cpp`
- Modify: `reseq/main.cpp` (remove helper bodies, add include)
- Modify: `reseq/CMakeLists.txt` (add target_sources)

Move `AutoDetectThreads`, `GetSeed`, `GetProbsOut`, `PrepareProbabilityEstimation`
from main.cpp into cli_common with their **exact current signatures**.

- [ ] **Step 1: Create reseq/cli/ directory**

```bash
mkdir -p reseq/cli
```

- [ ] **Step 2: Create reseq/cli/cli_common.h**

```cpp
#ifndef CLI_COMMON_H
#define CLI_COMMON_H

#include <string>

#include <boost/program_options.hpp>

#include "utilities.hpp"

namespace reseq::cli {

bool AutoDetectThreads(uintNumThreads& num_threads,
                       const boost::program_options::options_description& opt_desc,
                       const std::string& usage_str);

uintSeed GetSeed(const boost::program_options::variables_map& opts_map);

void GetProbsOut(std::string& probs_out, const std::string& fallback_out,
                 const boost::program_options::variables_map& opts_map);

void PrepareProbabilityEstimation(std::string& probs_in, std::string& probs_out,
                                  const std::string& standard_probs_out,
                                  bool loaded_stats,
                                  const boost::program_options::variables_map& opts_map);

} // namespace reseq::cli

#endif // CLI_COMMON_H
```

- [ ] **Step 3: Create reseq/cli/cli_common.cpp**

Copy the function bodies verbatim from `reseq/main.cpp` lines 66-82, 293-300, 302-321, 323-335. Add self-contained includes:

```cpp
#include "cli/cli_common.h"

#include <cstdio>
#include <string>
#include <thread>

#include <boost/program_options.hpp>

#include "logging.hpp"
#include "utilities.hpp"

using std::string;
using reseq::uintNumThreads;
using reseq::uintSeed;
using reseq::utilities::TrueRandom;

namespace reseq::cli {

bool AutoDetectThreads(uintNumThreads& num_threads,
                       const boost::program_options::options_description& opt_desc,
                       const string& usage_str) {
    // EXACT COPY from main.cpp lines 67-81
    if (0 == num_threads) {
        num_threads = std::thread::hardware_concurrency();
        if (0 == num_threads) {
            printErr << "Automatic detection of available cores failed." << std::endl;
            if (0 < kVerbosityLevel) {
                std::cerr << usage_str;
                std::cerr << opt_desc << std::endl;
            }
            return false;
        } else {
            printInfo << "Detected " << num_threads << " cores to be used." << std::endl;
        }
    }
    return true;
}

uintSeed GetSeed(const boost::program_options::variables_map& opts_map) {
    // EXACT COPY from main.cpp lines 323-335
    uintSeed seed;
    auto it = opts_map.find("seed");
    if (opts_map.end() == it) {
        seed = TrueRandom();
        printInfo << "Randomly generated seed is " << seed << std::endl;
    } else {
        seed = it->second.as<uintSeed>();
        printInfo << "Using seed " << seed << std::endl;
    }
    return seed;
}

void GetProbsOut(string& probs_out, const string& fallback_out,
                 const boost::program_options::variables_map& opts_map) {
    // EXACT COPY from main.cpp lines 293-300
    auto it_probs_out = opts_map.find("probabilitiesOut");
    if (opts_map.end() == it_probs_out) {
        probs_out = fallback_out;
    } else {
        probs_out = it_probs_out->second.as<string>();
    }
}

void PrepareProbabilityEstimation(string& probs_in, string& probs_out,
                                  const string& standard_probs_out,
                                  bool loaded_stats,
                                  const boost::program_options::variables_map& opts_map) {
    // EXACT COPY from main.cpp lines 302-321
    auto it_probs_in = opts_map.find("probabilitiesIn");
    if (opts_map.end() == it_probs_in) {
        GetProbsOut(probs_out, standard_probs_out, opts_map);

        probs_in = "";
        if (loaded_stats) {
            if (FILE* file = fopen(standard_probs_out.c_str(), "r")) {
                probs_in = standard_probs_out;
                fclose(file);
            }
        }
    } else {
        probs_in = it_probs_in->second.as<string>();
        GetProbsOut(probs_out, probs_in, opts_map);
    }
}

} // namespace reseq::cli
```

- [ ] **Step 4: Update main.cpp to include and call cli_common**

Add `#include "cli/cli_common.h"` after the existing includes (around line 48). Add `using reseq::cli::AutoDetectThreads;` etc. to maintain unqualified calls.

Remove the function bodies from main.cpp (lines 66-82, 293-335) but keep the call sites unchanged.

- [ ] **Step 5: Add to CMakeLists.txt**

After the `add_executable(reseq main.cpp)` line, add:

```cmake
target_sources(reseq PRIVATE cli/cli_common.cpp)
```

- [ ] **Step 6: Build and run tests**

```bash
make build && make test
```

- [ ] **Step 7: Commit**

```bash
git add reseq/cli/cli_common.h reseq/cli/cli_common.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "refactor(5a): extract shared CLI helpers into cli_common

Move AutoDetectThreads, GetSeed, GetProbsOut, PrepareProbabilityEstimation
from main.cpp to reseq/cli/cli_common with exact original signatures.
Self-contained header with own includes."
```

---

## Task 3: Extract replaceN Command

**Files:**
- Create: `reseq/cli/replace_n.h`
- Create: `reseq/cli/replace_n.cpp`
- Modify: `reseq/main.cpp` (replace inline code with function call)
- Modify: `reseq/CMakeLists.txt`

- [ ] **Step 1: Create reseq/cli/replace_n.h**

```cpp
#ifndef CLI_REPLACE_N_H
#define CLI_REPLACE_N_H

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "utilities.hpp"

namespace reseq::cli {

int RunReplaceN(const std::vector<std::string>& args,
                uintNumThreads num_threads,
                const boost::program_options::variables_map& general_opts,
                boost::program_options::options_description& opt_desc_full);

} // namespace reseq::cli

#endif // CLI_REPLACE_N_H
```

- [ ] **Step 2: Create reseq/cli/replace_n.cpp**

Move the replaceN command block (main.cpp lines 606-681) into `RunReplaceN`. The mode banner suffix and command-name erase are now handled by main.cpp before calling, so the Run function starts with option parsing. Copy the exact code:

```cpp
#include "cli/replace_n.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "cli/cli_common.h"
#include "logging.hpp"
#include "Reference.h"
#include "utilities.hpp"

using std::cerr;
using std::exception;
using std::string;
using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using reseq::Reference;
using reseq::uintNumThreads;
using reseq::uintSeed;

namespace reseq::cli {

int RunReplaceN(const std::vector<std::string>& args,
                uintNumThreads num_threads,
                const variables_map& general_opts,
                options_description& opt_desc_full) {
    // EXACT COPY from main.cpp lines 612-681
    options_description opt_desc("ReplaceN");
    opt_desc.add_options()
        ("refIn,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")
        ("refSim,R", value<string>(),
         "File to where reference sequences in fasta format with N's randomly replace should be written to")
        ("seed", value<uintSeed>(), "Seed used for replacing N, if none is given random seed will be used");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq replaceN -r <refIn.fa> -R <refSim.fa> [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse replaceN command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
    } else if (!AutoDetectThreads(num_threads, opt_desc_full, usage_str)) {
        return 1;
    } else {
        auto it_ref_in = opts_map.find("refIn");
        auto it_ref_out = opts_map.find("refSim");
        string ref_input, ref_output;
        if (opts_map.end() == it_ref_in) {
            printErr << "refIn option is mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc_full << std::endl;
            }
            return 1;
        } else {
            ref_input = it_ref_in->second.as<string>();
            printInfo << "Reading reference from " << ref_input << std::endl;

            if (opts_map.end() == it_ref_out) {
                printErr << "refSim option is mandatory." << std::endl;
                if (0 < kVerbosityLevel) {
                    cerr << usage_str;
                    cerr << opt_desc_full << std::endl;
                }
                return 1;
            } else {
                ref_output = it_ref_out->second.as<string>();
                printInfo << "Writing reference without N to " << ref_output << std::endl;

                Reference species_reference;
                if (species_reference.ReadFasta(ref_input.c_str())) {
                    auto seed = GetSeed(opts_map);
                    species_reference.ReplaceN(seed);

                    if (species_reference.WriteFasta(ref_output.c_str())) {
                        printInfo << "Finished replacing N's." << std::endl;
                    } else {
                        return 1;
                    }
                } else {
                    return 1;
                }
            }
        }
    }

    return 0;
}

} // namespace reseq::cli
```

- [ ] **Step 3: Update main.cpp dispatch for replaceN**

Replace the replaceN block (lines 606-681) in main.cpp with:

```cpp
        } else if ("replaceN" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in replaceN mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = cli::RunReplaceN(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
```

Add `#include "cli/replace_n.h"` at the top.

- [ ] **Step 4: Add to CMakeLists.txt**

Add `cli/replace_n.cpp` to `target_sources`.

- [ ] **Step 5: Build and run tests**

```bash
make build && make test
```

- [ ] **Step 6: Commit**

```bash
git add reseq/cli/replace_n.h reseq/cli/replace_n.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "refactor(5a): extract replaceN command to cli/replace_n

Move replaceN option parsing and logic into RunReplaceN().
main.cpp handles mode banner and command-name erase before calling."
```

---

## Task 4: Extract queryProfile Command

**Files:**
- Create: `reseq/cli/query_profile.h`
- Create: `reseq/cli/query_profile.cpp`
- Modify: `reseq/main.cpp`
- Modify: `reseq/CMakeLists.txt`

Same pattern as Task 3. Move main.cpp lines 455-605 into `RunQueryProfile`.

- [ ] **Step 1: Create reseq/cli/query_profile.h**

Same interface pattern as replace_n.h but with `RunQueryProfile`.

- [ ] **Step 2: Create reseq/cli/query_profile.cpp**

Move the queryProfile block verbatim. Self-contained includes — this command uses `DataStats`, `FragmentDistributionStats`, `ProbabilityEstimates`, `Reference`, and `Simulator` headers. Add all needed includes explicitly.

The queryProfile command calls no shared helpers from cli_common — it's self-contained.

- [ ] **Step 3: Update main.cpp dispatch**

Replace the queryProfile block with:

```cpp
        if ("queryProfile" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in queryProfile mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = cli::RunQueryProfile(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
```

- [ ] **Step 4: Add to CMakeLists.txt, build and test**

```bash
make build && make test
```

- [ ] **Step 5: Commit**

```bash
git add reseq/cli/query_profile.h reseq/cli/query_profile.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "refactor(5a): extract queryProfile command to cli/query_profile"
```

---

## Task 5: Extract seqToIllumina Command

**Files:**
- Create: `reseq/cli/seq_to_illumina.h`
- Create: `reseq/cli/seq_to_illumina.cpp`
- Modify: `reseq/main.cpp`
- Modify: `reseq/CMakeLists.txt`

Move main.cpp lines 1019-1153. This command uses `PrepareProbabilityEstimation` and `GetSeed` from cli_common.

- [ ] **Step 1: Create header and implementation**

Same pattern. Self-contained includes. The implementation calls `cli::PrepareProbabilityEstimation` and `cli::GetSeed`.

Note: seqToIllumina's missing-statsIn error path (line 1080-1085) prints `usage_str` + `opt_desc` (command-only), NOT `opt_desc_full`. Preserve this exactly.

- [ ] **Step 2: Update main.cpp dispatch**

```cpp
        } else if ("seqToIllumina" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in seqToIllumina mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = cli::RunSeqToIllumina(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
```

- [ ] **Step 3: Add to CMakeLists.txt, build and test**

```bash
make build && make test
```

- [ ] **Step 4: Commit**

```bash
git add reseq/cli/seq_to_illumina.h reseq/cli/seq_to_illumina.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "refactor(5a): extract seqToIllumina command to cli/seq_to_illumina"
```

---

## Task 6: Extract illuminaPE Command

**Files:**
- Create: `reseq/cli/illumina_pe.h`
- Create: `reseq/cli/illumina_pe.cpp`
- Modify: `reseq/main.cpp`
- Modify: `reseq/CMakeLists.txt`

This is the largest extraction (~337 lines). It carries 4 file-local helpers:
`DefaultExtensionFile`, `GetDataStats`, `WriteSysError`, `PrepareSimulation`.

- [ ] **Step 1: Create reseq/cli/illumina_pe.h**

Same interface pattern.

- [ ] **Step 2: Create reseq/cli/illumina_pe.cpp**

Self-contained includes. The 4 private helpers go in an unnamed namespace at the top:

```cpp
#include "cli/illumina_pe.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
// ... all needed includes ...

#include "cli/cli_common.h"
#include "DataStats.h"
#include "FragmentDistributionStats.h"
#include "logging.hpp"
#include "ProbabilityEstimates.h"
#include "Reference.h"
#include "Simulator.h"
#include "utilities.hpp"

// ... using declarations ...

namespace {

// EXACT COPY from main.cpp lines 84-117
bool DefaultExtensionFile(std::string& file_name, const std::string extension) {
    // ... verbatim ...
}

// EXACT COPY from main.cpp lines 119-291
void GetDataStats(reseq::DataStats& real_data_stats, std::string& stats_file,
                  bool& loaded_stats, bool stats_only,
                  reseq::uintSeqLen max_ref_seq_bin_size,
                  const boost::program_options::variables_map& opts_map,
                  const boost::program_options::options_description& opt_desc,
                  const std::string& usage_str,
                  reseq::uintNumThreads num_threads) {
    // ... verbatim ...
}

// EXACT COPY from main.cpp lines 337-388
void WriteSysError(...) {
    // ... verbatim ...
}

// EXACT COPY from main.cpp lines 390-405
void PrepareSimulation(...) {
    // ... verbatim ...
}

} // anonymous namespace

namespace reseq::cli {

int RunIlluminaPE(const std::vector<std::string>& args,
                  uintNumThreads num_threads,
                  const variables_map& general_opts,
                  options_description& opt_desc_full) {
    // EXACT COPY from main.cpp lines 688-1018
    // IMPORTANT: Initialize loaded_stats defensively at line 832 equivalent:
    bool loaded_stats = false;  // defensive init (see design spec)
    // ... rest verbatim ...
}

} // namespace reseq::cli
```

**Critical:** When copying `GetDataStats`, ensure `loaded_stats` is initialized to `false`
at the call site (line 832 equivalent) before passing to `GetDataStats`. This is a
defensive initialization that does not change behavior.

- [ ] **Step 3: Update main.cpp dispatch**

```cpp
        } else if ("illuminaPE" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in illuminaPE mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = cli::RunIlluminaPE(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
```

- [ ] **Step 4: Add to CMakeLists.txt, build and test**

```bash
make build && make test
```

- [ ] **Step 5: Commit**

```bash
git add reseq/cli/illumina_pe.h reseq/cli/illumina_pe.cpp reseq/main.cpp reseq/CMakeLists.txt
git commit -m "refactor(5a): extract illuminaPE command to cli/illumina_pe

Carries DefaultExtensionFile, GetDataStats, WriteSysError, PrepareSimulation
as file-local helpers in unnamed namespace. Initialize loaded_stats=false
defensively at call site."
```

---

## Task 7: Clean Up main.cpp and Final Verification

**Files:**
- Modify: `reseq/main.cpp` (remove dead includes, dead helper forward declarations)

- [ ] **Step 1: Remove dead includes and using declarations from main.cpp**

After all commands are extracted, main.cpp no longer needs most of its includes.
Remove includes that are only used by extracted command code (DataStats.h, Simulator.h,
ProbabilityEstimates.h, FragmentDistributionStats.h, etc.). Keep:
- `<atomic>`, `<iostream>`, `<string>`, `<vector>`
- `<boost/program_options.hpp>`
- `"logging.hpp"`, `"CMakeConfig.h"`, `"utilities.hpp"`
- The 4 `cli/*.h` headers
- SeqAn includes for the complement statics

Remove any remaining helper function bodies or forward declarations that were already
moved to cli_common or the command files.

- [ ] **Step 2: Verify main.cpp line count**

```bash
wc -l reseq/main.cpp
```

Expected: ~100 lines (globals + general options + dispatch).

- [ ] **Step 3: Run format check**

```bash
make format-check
```

If formatting drift, run `make format`.

- [ ] **Step 4: Build and run full test suite**

```bash
make build && make test
```

All tests must pass.

- [ ] **Step 5: Run ASan+UBSan**

```bash
rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build -j$(nproc) && cd build && ctest --output-on-failure -j1
```

- [ ] **Step 6: Manual smoke tests**

```bash
# Verify help output for commands without strong regression coverage
build/bin/reseq illuminaPE --help 2>&1 | head -5
build/bin/reseq seqToIllumina --help 2>&1 | head -5
build/bin/reseq --verbosity 3 replaceN --help 2>&1 | head -3
build/bin/reseq 2>&1 | head -3
build/bin/reseq nonsenseCmd 2>&1 | head -3
```

Verify output matches expectations from the design spec behavior table.

- [ ] **Step 7: Commit**

```bash
git add reseq/main.cpp
git commit -m "refactor(5a): clean up main.cpp after extraction

Remove dead includes and declarations. main.cpp is now a slim
~100-line dispatcher: globals, general options, version/help, and
Run* calls for each command."
```

---

## Verification Protocol

**After every task:**

```bash
make build && make test
```

All unit tests + regression tests must pass.

**After Task 7 (final):**

- Normal build + test
- ASan+UBSan clean
- Format check clean
- Manual smoke tests for all 4 commands + help + version + unknown command
- `wc -l reseq/main.cpp` ≈ 100 lines
