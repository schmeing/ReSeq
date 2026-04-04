# Phase 5a: main.cpp Decomposition — Design Specification

**Date:** 2026-04-03
**Status:** Approved (rev 4)
**Phase:** 5a (God Object Decomposition — main.cpp)
**Prerequisite:** Phase 4 complete (PR #7, `2e65afb`)

---

## 1. Goal

Decompose `reseq/main.cpp` (1,167 lines) into a slim dispatcher (~100 lines) with
command logic extracted into focused files under `reseq/cli/`.

This is a **behavior-preserving extraction with adapter wrappers** — helper functions
are moved with their exact current signatures, command logic is wrapped in Run*
functions, and each new `.cpp` file gets self-contained includes. It is NOT a verbatim
copy-paste across translation units (the monolithic include preamble must be decomposed
per file).

## 2. Current Architecture

### File Structure

```
reseq/main.cpp (1,167 lines)
  Lines 1-64:     Includes, globals, SeqAn statics
  Lines 66-82:    AutoDetectThreads(num_threads, opt_desc, usage_str) → bool
  Lines 84-117:   DefaultExtensionFile(file_name, extension) → bool (mutates in-place)
  Lines 119-291:  GetDataStats(...) — illuminaPE only
  Lines 293-300:  GetProbsOut(probs_out, fallback_out, opts_map) → void
  Lines 302-321:  PrepareProbabilityEstimation(probs_in, probs_out, standard_probs_out,
                                               loaded_stats, opts_map) → void
  Lines 323-335:  GetSeed(opts_map) → uintSeed
  Lines 337-388:  WriteSysError(...) — illuminaPE only
  Lines 390-405:  PrepareSimulation(...) — illuminaPE only
  Lines 408-453:  main() — general options, version, help, dispatch
  Lines 455-605:  queryProfile command (~150 lines)
  Lines 606-681:  replaceN command (~75 lines)
  Lines 682-1018: illuminaPE command (~337 lines)
  Lines 1019-1153: seqToIllumina command (~135 lines)
  Lines 1154-1167: unknown command error, return
```

### Exact Helper Signatures (moved verbatim)

```cpp
// Line 66 — returns bool, prints help on failure internally
bool AutoDetectThreads(uintNumThreads& num_threads,
                       const options_description& opt_desc,
                       const string& usage_str);

// Line 84 — returns bool, mutates file_name in-place
bool DefaultExtensionFile(string& file_name, const string extension);

// Line 293 — void, output via probs_out reference
void GetProbsOut(string& probs_out, const string& fallback_out,
                 const variables_map& opts_map);

// Line 302 — void, needs standard_probs_out AND loaded_stats
void PrepareProbabilityEstimation(string& probs_in, string& probs_out,
                                  const string& standard_probs_out,
                                  bool loaded_stats,
                                  const variables_map& opts_map);

// Line 323 — takes only opts_map, returns uintSeed
uintSeed GetSeed(const variables_map& opts_map);
```

### Global Symbol Ownership

| Symbol | Declared In | Defined In |
|--------|-------------|------------|
| `kVerbosityLevel` | `logging.hpp:24` | `main.cpp:20`, `test_main.cpp:8` |
| `kNoDebugOutput` | `logging.hpp:25` | `main.cpp:22`, `test_main.cpp:9` |
| SeqAn complement statics | N/A | `main.cpp:62-63`, `test_main.cpp:34-35` |

Not changed in this phase.

### Exact Current Dispatch Behavior (verbosity-aware)

The version banner and mode suffix use DIFFERENT verbosity thresholds:
- `printInfo` (line 452) emits at verbosity >= 4 (`logging.hpp:243-244`)
- Mode suffix (line 456-458) and unknown-command newline (line 1155-1156) check
  `2 < kVerbosityLevel` (verbosity >= 3) and write to `cerr` directly

This means at `--verbosity 3`: no banner prefix (printInfo → kNullStream),
but `" in <mode> mode\n"` DOES go to cerr as an orphaned suffix. This is the
current behavior and must be preserved exactly.

| Scenario | Verbosity 4+ | Verbosity 3 | Verbosity 0-2 |
|----------|-------------|-------------|---------------|
| Valid command | banner + " in mode\n" | orphaned " in mode\n" | nothing |
| Unknown command | banner + newline + error + usage | newline (no banner) + error + usage | error only |
| No command | general_usage, exit 0 | general_usage, exit 0 | nothing, exit 0 |

Verified at runtime on 2026-04-03:
- `reseq` (no args) and `reseq --help` produce **identical output**, both exit 0
- No separate `--help` branch at top level (line 449-450 just prints general_usage)

### Known Latent Issue: loaded_stats

`GetDataStats` (`main.cpp:119`) does not initialize `loaded_stats` on all error paths.
The illuminaPE caller (`main.cpp:832-844`) is currently safe because error paths exit
before consuming it. During extraction, **initialize `loaded_stats = false` at the
call site** defensively. This does not change behavior but prevents future UB if
control flow is ever reordered.

### Helper Placement

| Helper | Callers | Placement |
|--------|---------|-----------|
| `AutoDetectThreads` | replaceN, illuminaPE, seqToIllumina | `cli_common.cpp` |
| `GetSeed` | replaceN, WriteSysError, seqToIllumina | `cli_common.cpp` |
| `GetProbsOut` | PrepareProbabilityEstimation | `cli_common.cpp` (co-located) |
| `PrepareProbabilityEstimation` | illuminaPE, seqToIllumina | `cli_common.cpp` |
| `DefaultExtensionFile` | GetDataStats only | `illumina_pe.cpp` (unnamed namespace) |
| `GetDataStats` | illuminaPE only | `illumina_pe.cpp` (unnamed namespace) |
| `WriteSysError` | illuminaPE only | `illumina_pe.cpp` (unnamed namespace) |
| `PrepareSimulation` | illuminaPE only | `illumina_pe.cpp` (unnamed namespace) |

File-local helpers use unnamed namespaces in `.cpp` files, NOT header declarations.

### Per-Command Error/Help Patterns (not normalized)

Each command has its own help/error output pattern. These are NOT identical and must
NOT be normalized during extraction:
- queryProfile: prints `usage_str` + `opt_desc` (command-only options)
- replaceN: prints `usage_str` + `opt_desc` (command-only options)
- illuminaPE: prints `usage_str` + `opt_desc_full` (all option groups combined)
- seqToIllumina: prints `usage_str` + `opt_desc` on missing-statsIn (line 1080-1085),
  `opt_desc_full` for `--help` (line 1066)

Each Run function preserves its own pattern verbatim.

## 3. Target Architecture

### File Structure

```
reseq/main.cpp                (~100 lines) — globals, SeqAn statics, general opts, dispatch
reseq/cli/cli_common.h        (~40 lines)  — self-contained: own includes for exposed types
reseq/cli/cli_common.cpp      (~100 lines) — implementations
reseq/cli/query_profile.h     (~15 lines)  — self-contained header
reseq/cli/query_profile.cpp   (~160 lines) — own includes, queryProfile implementation
reseq/cli/replace_n.h         (~15 lines)  — self-contained header
reseq/cli/replace_n.cpp       (~85 lines)  — own includes, replaceN implementation
reseq/cli/illumina_pe.h       (~15 lines)  — self-contained header
reseq/cli/illumina_pe.cpp     (~420 lines) — own includes, illuminaPE + 4 unnamed-ns helpers
reseq/cli/seq_to_illumina.h   (~15 lines)  — self-contained header
reseq/cli/seq_to_illumina.cpp (~140 lines) — own includes, seqToIllumina implementation
```

**Self-contained headers:** Each `.h` file includes exactly the types it exposes
(e.g., `<vector>`, `<string>`, `<boost/program_options.hpp>`, `"utilities.hpp"`).
Each `.cpp` file includes exactly what it needs (not the full main.cpp preamble).

### Command Interface

```cpp
namespace reseq::cli {

// args:          unrecognized_opts with command name already erased by main
// num_threads:   from general options (may be 0 = auto-detect)
// general_opts:  for checking help flag
// opt_desc_full: mutable ref — command adds its options for combined help display
// Returns:       0 success, 1 error

int RunQueryProfile(const std::vector<std::string>& args,
                    uintNumThreads num_threads,
                    const boost::program_options::variables_map& general_opts,
                    boost::program_options::options_description& opt_desc_full);

// Same signature for RunReplaceN, RunIlluminaPE, RunSeqToIllumina
}
```

### Dispatch Ownership

| Action | Who Does It |
|--------|-------------|
| Start version banner (`printInfo << "Running ReSeq version..."`) | main.cpp |
| Erase command name from `unrecognized_opts` | main.cpp (before calling Run) |
| Append mode suffix (`cerr << " in <mode> mode" << endl`) | Run function |
| Per-command option parsing | Run function |
| Per-command help/error output (each pattern preserved verbatim) | Run function |

### CMake Changes

Use `target_sources` for incremental addition rather than rewriting `add_executable`:

```cmake
# Existing:
add_executable(reseq main.cpp)
target_link_libraries(reseq PRIVATE reseq_lib)

# Add after each extraction step:
target_sources(reseq PRIVATE
  cli/cli_common.cpp
  cli/query_profile.cpp
  cli/replace_n.cpp
  cli/illumina_pe.cpp
  cli/seq_to_illumina.cpp
)
```

The `reseq` target already inherits include paths from `reseq_lib` via
`target_link_libraries(reseq PRIVATE reseq_lib)` — `reseq_lib` has PUBLIC
`target_include_directories` (`reseq/CMakeLists.txt:25-29`) which propagate to
`reseq`. This means `cli/*.cpp` can use `#include "Simulator.h"` etc. without
additional CMake configuration.

## 4. Migration Strategy

### Step 0: Strengthen regression tests (prerequisite)

Before any code moves, add:
1. CLI behavior tests (exit codes, stderr checks)
2. At least one direct `seqToIllumina` smoke test
3. One `illuminaPE` test that exercises beyond `--statsOnly --noBias`

See Testing section for details.

### Steps 1-5: Extract one piece at a time

Each step: create files, move code, add self-contained includes, build, test, commit.

1. `cli_common` (shared helpers)
2. `replace_n` (simplest command)
3. `query_profile` (self-contained, no shared helpers beyond basics)
4. `seq_to_illumina` (uses PrepareProbabilityEstimation from cli_common)
5. `illumina_pe` (largest, carries 4 unnamed-ns helpers, initialize `loaded_stats = false`)

### Step 6: Clean up main.cpp

Remove dead includes, verify line count, format.

## 5. Testing

### Actual Regression Coverage (honest assessment)

| Command | Test Coverage | Level |
|---------|--------------|-------|
| `replaceN` | `ReplaceN` (full stdout diff) | **Strong** |
| `queryProfile` | 6 tests (3 generated, 3 Zenodo) | **Strong** |
| `illuminaPE` | `GenerateEcoliProfile()` with `--statsOnly --noBias` | **Partial** — stats path only |
| `seqToIllumina` | **None** | **No coverage** |
| CLI behavior | `ErrorBadCommand` (exit != 0), `VersionOutput` (stderr substring) | **Weak** |

### New Tests (prerequisite — BEFORE extraction)

**CLI behavior tests** (exit codes + stderr checks):

| Test | Assertion |
|------|-----------|
| `BareReseqExitCode` | `reseq` → exit 0 |
| `BareReseqOutput` | `reseq` stderr contains "reseq <command>" (usage presence) |
| `HelpSameAsBare` | `reseq --help` stderr == bare `reseq` stderr |
| `UnknownCommandExitCode` | `reseq badcmd` → exit 1 |
| `UnknownCommandStderr` | stderr contains "Unrecognized command: 'badcmd'" |
| `CommandHelpExitCode` | `reseq replaceN --help` → exit 0 |
| `VersionExitCode` | `reseq --version` → exit 0 (strengthens existing VersionOutput) |

**seqToIllumina smoke test:**

| Test | What It Covers |
|------|---------------|
| `SeqToIlluminaSmoke` | Run `seqToIllumina` with test data, verify exit 0 and output files created |

**illuminaPE extended coverage** (if feasible with test data):

| Test | What It Covers |
|------|---------------|
| `IlluminaPESimulationSmoke` | Run `illuminaPE` with `--stopAfterEstimation` or minimal simulation, verify exit 0 |

If the extended illuminaPE test requires too much test data or runtime, document it as
a manual smoke test step in the verification protocol instead.

### Test Helpers Needed

Add to `RegressionTest.h`:
- `RunReseqExitOnly(args)` — runs WITHOUT `--verbosity 0` or `2>/dev/null`, returns exit code
- Leverage existing `RunReseqCaptureStderr(args)` for stderr checks

## 6. Verification Criteria

At every commit:
- `make build && make test` — all tests pass
- `make format-check` — clean

After final commit:
- `wc -l reseq/main.cpp` — target ~100 lines
- All new + existing regression tests pass
- Manual: `reseq illuminaPE --help` and `reseq seqToIllumina --help` produce correct output
- Manual: `reseq --verbosity 3 replaceN --help` preserves orphaned suffix behavior
