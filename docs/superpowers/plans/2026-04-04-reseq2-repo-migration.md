# ReSeq2 Repository Migration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the berntpopp/ReSeq fork into a standalone berntpopp/reseq2 repository with updated branding, versioning, licensing, and packaging as a new independent project.

**Architecture:** Create a new GitHub repo from the existing fork (preserving full git history), bump to version 2.0.0, update all branding/attribution, and prepare a Bioconda recipe. The existing berntpopp/ReSeq fork stays as-is for reference.

**Tech Stack:** Git, GitHub CLI (`gh`), CMake, Bioconda, JOSS (future)

---

## Context

The upstream `schmeing/ReSeq` has been inactive since April 2021 with unanswered issues from 2023-2024. The Bioconda package is stale. This fork has undergone significant refactoring (C++20, CMake modernization, SeqAn externalized, CI/CD, test coverage). The MIT license permits unrestricted forking and renaming.

## Pre-Migration Checklist

Before starting, ensure all in-progress refactoring work on `berntpopp/ReSeq` is complete and merged to `master`.

---

### Task 1: Create the New GitHub Repository

**Files:** None (GitHub operations only)

- [ ] **Step 1: Create a bare clone of the current fork**

```bash
cd /tmp
git clone --bare https://github.com/berntpopp/ReSeq.git reseq2-bare
```

- [ ] **Step 2: Create the new repo on GitHub**

```bash
gh repo create berntpopp/reseq2 --public --description "ReSeq2 - Realistic Illumina sequencing simulator (successor to ReSeq)"
```

- [ ] **Step 3: Push the bare clone to the new repo**

```bash
cd /tmp/reseq2-bare
git push --mirror https://github.com/berntpopp/reseq2.git
```

- [ ] **Step 4: Clone the new repo for further work**

```bash
cd ~/development
git clone https://github.com/berntpopp/reseq2.git
cd reseq2
```

- [ ] **Step 5: Verify git history is intact**

```bash
git log --oneline | head -20
git log --oneline | wc -l
```

Expected: Full commit history matching the original fork.

- [ ] **Step 6: Commit** (no commit needed -- repo setup only)

---

### Task 2: Bump Version to 2.0.0

**Files:**
- Modify: `VERSION`

- [ ] **Step 1: Update VERSION file**

Change the content of `VERSION` from:
```
1.1.0
```
to:
```
2.0.0
```

- [ ] **Step 2: Verify CMake reads it correctly**

```bash
mkdir -p build && cd build
cmake .. 2>&1 | head -5
```

Expected: No errors. The version macros in `CMakeConfig.h.in` will pick up `RESEQ_VERSION_MAJOR=2`, `RESEQ_VERSION_MINOR=0`, `RESEQ_VERSION_PATCH=0` at configure time.

- [ ] **Step 3: Commit**

```bash
git add VERSION
git commit -m "build: bump version to 2.0.0 for ReSeq2 launch"
```

---

### Task 3: Update LICENSE with Dual Copyright

**Files:**
- Modify: `LICENSE`

- [ ] **Step 1: Update LICENSE**

Replace the copyright line:
```
Copyright (c) 2019 schmeing
```
with:
```
Copyright (c) 2019 Stephan Schmeing (original ReSeq)
Copyright (c) 2024-2026 Bernt Popp (ReSeq2 continuation)
```

Leave the rest of the MIT license text unchanged.

- [ ] **Step 2: Verify the full LICENSE reads correctly**

```bash
cat LICENSE
```

Expected: Dual copyright header followed by standard MIT license text.

- [ ] **Step 3: Commit**

```bash
git add LICENSE
git commit -m "docs: add dual copyright for ReSeq2 continuation"
```

---

### Task 4: Update README.md

**Files:**
- Modify: `README.md`

This task rewrites the README to reflect ReSeq2 branding while preserving all usage documentation. The key changes are: title, badges, attribution section, installation instructions, and repository URLs.

- [ ] **Step 1: Update the title and description**

Change the first lines:
```markdown
# ReSeq

[![codecov](https://codecov.io/gh/berntpopp/ReSeq/graph/badge.svg)](https://codecov.io/gh/berntpopp/ReSeq)

More realistic simulator for genomic DNA sequences from Illumina machines that achieves a similar k-mer spectrum as the original sequences.
```
to:
```markdown
# ReSeq2

[![codecov](https://codecov.io/gh/berntpopp/reseq2/graph/badge.svg)](https://codecov.io/gh/berntpopp/reseq2)

Realistic simulator for genomic DNA sequences from Illumina paired-end sequencers. ReSeq2 learns error, quality, and coverage profiles from real data and uses them to generate synthetic reads with matching k-mer spectra.

> **ReSeq2** is a maintained continuation of [ReSeq](https://github.com/schmeing/ReSeq) by Schmeing & Robinson, originally published in [Genome Biology (2021)](https://doi.org/10.1186/s13059-021-02265-7). The original repository is no longer maintained. See [Attribution](#attribution) for details.
```

- [ ] **Step 2: Add Attribution section before Publication**

Insert a new section before the existing `## Publication` section:

```markdown
## <a name="attribution"></a>Attribution

ReSeq2 is a fork and continuation of the original [ReSeq](https://github.com/schmeing/ReSeq) project by Stephan Schmeing and Mark D. Robinson at the University of Zurich. The original project is no longer actively maintained.

**What changed in ReSeq2:**
- Migrated to C++20
- Modernized CMake build system with external dependency management (SeqAn 2.5.2 via FetchContent)
- Added CI/CD pipeline (GitHub Actions) with multi-compiler builds, sanitizers, and code coverage
- Expanded test suite with golden-file, edge-case, and component-level tests
- Added binary/text profile format conversion (`convertProfile` command)
- Code quality improvements: clang-format, clang-tidy, pre-commit hooks

The original work is gratefully acknowledged. Please cite both the original paper and this project when using ReSeq2.
```

- [ ] **Step 3: Update Publication section**

Change the existing Publication section to:

```markdown
## <a name="publication"></a>Publication

**Original ReSeq paper (please cite):**
Schmeing, S., Robinson, M.D. ReSeq simulates realistic Illumina high-throughput sequencing data. Genome Biol 22, 67 (2021). https://doi.org/10.1186/s13059-021-02265-7
```

- [ ] **Step 4: Update installation instructions**

In the Installation section, change the clone URL:
```
git clone https://github.com/schmeing/ReSeq.git
```
to:
```
git clone https://github.com/berntpopp/reseq2.git
```

And update the Bioconda section to note the new package name:
```markdown
## <a name="conda"></a>Bioconda
ReSeq2 will be available via Bioconda (package name: `reseq2`). The original `reseq` package on Bioconda installs the unmaintained upstream version.
```

- [ ] **Step 5: Add Attribution to Table of Contents**

Add `- [Attribution](#attribution)` to the table of contents, before the Publication entry.

- [ ] **Step 6: Update Requirements table**

Change the compiler requirement from:
```
Compiler supporting C++14
```
to:
```
Compiler supporting C++20
```

Remove the CentOS 7 / old GCC workaround text. Update the tested version column where applicable (e.g., CMake 3.16+, GCC 10+/Clang 12+).

- [ ] **Step 7: Update all CLI command references**

Replace all `reseq ` command invocations with `reseq2 ` throughout the README (quick start examples, parameter tables, FAQ, etc.). This includes:
- `reseq illuminaPE` -> `reseq2 illuminaPE`
- `reseq seqToIllumina` -> `reseq2 seqToIllumina`
- `reseq queryProfile` -> `reseq2 queryProfile`
- `reseq replaceN` -> `reseq2 replaceN`
- `reseq convertProfile` -> `reseq2 convertProfile`
- `reseq test` -> `reseq2 test`
- Python script references (`reseq-prepare-names.py` -> `reseq2-prepare-names.py` if renamed, or keep if not)

- [ ] **Step 8: Verify README renders correctly**

```bash
# Quick check for broken links or markdown issues
grep -n 'schmeing/ReSeq' README.md
```

Expected: Only the attribution/publication sections should reference the original repo. No stale URLs in installation or badge sections.

- [ ] **Step 9: Commit**

```bash
git add README.md
git commit -m "docs: rebrand README for ReSeq2 with attribution"
```

---

### Task 5: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update repository description**

Change:
```markdown
This is a fork (`berntpopp/ReSeq`) of upstream (`schmeing/ReSeq`). The upstream remote has been removed.

**All work stays in the fork.** Never add the upstream remote, create PRs on `schmeing/ReSeq`, or push to it. Use `--repo berntpopp/ReSeq` with all `gh` commands.
```
to:
```markdown
This is **ReSeq2** (`berntpopp/reseq2`), a maintained continuation of the original ReSeq (`schmeing/ReSeq`). The original project is no longer maintained.

Use `--repo berntpopp/reseq2` with all `gh` commands.
```

- [ ] **Step 2: Update the Architecture section header**

Change:
```markdown
ReSeq is a bioinformatics tool that learns error/quality profiles
```
to:
```markdown
ReSeq2 is a bioinformatics tool that learns error/quality profiles
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for ReSeq2 repo"
```

---

### Task 6: Update CMake Project Name and CI

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Update CMake project name**

In `CMakeLists.txt` line 12, change:
```cmake
project(reseq LANGUAGES CXX)
```
to:
```cmake
project(reseq2 LANGUAGES CXX)
```

- [ ] **Step 2: Update binary/target names to reseq2**

```bash
grep -rn 'add_executable\|add_library\|OUTPUT_NAME\|reseq_lib\|reseq_test\|target_link' CMakeLists.txt reseq/CMakeLists.txt
```

Rename targets and output binaries: `reseq` -> `reseq2`, `reseq_lib` -> `reseq2_lib`, `reseq_test` -> `reseq2_test`. Update all `target_link_libraries`, `add_executable`, and `add_library` references accordingly. The installed binary should be `reseq2` so it doesn't conflict with the old `reseq` Bioconda package.

- [ ] **Step 3: Update CI badge references if present**

Check `.github/workflows/ci.yml` for any hardcoded references to `berntpopp/ReSeq` and update to `berntpopp/reseq2`.

- [ ] **Step 4: Build and test**

```bash
cd build
cmake ..
make -j$(nproc)
ctest --output-on-failure
```

Expected: Clean build, all tests pass, binary is named `reseq2`.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt .github/workflows/ci.yml
git commit -m "build: rename CMake project to reseq2"
```

---

### Task 7: Update Codecov and CI Integration

**Files:**
- Modify: `.github/workflows/ci.yml` (if Codecov slug is hardcoded)

- [ ] **Step 1: Search for hardcoded repo references**

```bash
grep -rn 'berntpopp/ReSeq' .github/ README.md .codecov.yml 2>/dev/null
```

- [ ] **Step 2: Update all references to the new repo name**

Replace all instances of `berntpopp/ReSeq` with `berntpopp/reseq2` in CI and configuration files. Be careful not to change references to `schmeing/ReSeq` in attribution sections.

- [ ] **Step 3: Commit**

```bash
git add -A .github/ .codecov.yml
git commit -m "ci: update repo references to berntpopp/reseq2"
```

---

### Task 8: Create Git Tag and GitHub Release

**Files:** None (git/GitHub operations only)

- [ ] **Step 1: Create annotated tag**

```bash
git tag -a v2.0.0 -m "ReSeq2 v2.0.0 - first independent release

Maintained continuation of ReSeq by Schmeing & Robinson.
Major changes from upstream v1.1:
- C++20 migration
- Modernized CMake build with external dependency management
- CI/CD with multi-compiler builds, sanitizers, coverage
- Expanded test suite
- Binary/text profile format conversion
- Code quality tooling (clang-format, clang-tidy, pre-commit)"
```

- [ ] **Step 2: Push tag and branch**

```bash
git push origin master --tags
```

- [ ] **Step 3: Create GitHub release**

```bash
gh release create v2.0.0 \
  --repo berntpopp/reseq2 \
  --title "ReSeq2 v2.0.0" \
  --notes "$(cat <<'EOF'
## ReSeq2 v2.0.0

First independent release of ReSeq2, a maintained continuation of [ReSeq](https://github.com/schmeing/ReSeq) by Schmeing & Robinson ([Genome Biology, 2021](https://doi.org/10.1186/s13059-021-02265-7)).

### What's new compared to upstream ReSeq v1.1

- **C++20** — modernized language standard
- **CMake overhaul** — external dependency management via FetchContent (SeqAn 2.5.2, GoogleTest, NLopt)
- **CI/CD** — GitHub Actions with GCC 13 / Clang 17, ASan+UBSan, code coverage via Codecov
- **Test suite** — golden-file tests, edge-case tests, component-level tests
- **`convertProfile` command** — convert between binary and portable text profile formats
- **Code quality** — clang-format, clang-tidy, ruff, pre-commit hooks

### Installation

```bash
git clone https://github.com/berntpopp/reseq2.git
cd reseq2
mkdir build && cd build
cmake ..
make -j$(nproc)
```

Bioconda package (`reseq2`) coming soon.

### Attribution

ReSeq2 builds on the original work of Stephan Schmeing and Mark D. Robinson. Please cite the original paper when using ReSeq2.
EOF
)"
```

- [ ] **Step 4: Verify release page**

```bash
gh release view v2.0.0 --repo berntpopp/reseq2
```

Expected: Release page shows correctly with all notes.

---

### Task 9: Prepare Bioconda Recipe (Draft)

**Files:**
- Create: `conda/meta.yaml` (draft, to be submitted to bioconda-recipes later)

- [ ] **Step 1: Create draft Bioconda recipe**

Create `conda/meta.yaml`:

```yaml
{% set name = "reseq2" %}
{% set version = "2.0.0" %}

package:
  name: {{ name }}
  version: {{ version }}

source:
  url: https://github.com/berntpopp/reseq2/archive/v{{ version }}.tar.gz
  sha256: PLACEHOLDER_UPDATE_AFTER_RELEASE

build:
  number: 0
  skip: true  # [not linux]
  run_exports:
    - {{ pin_subpackage(name, max_pin="x.x") }}

requirements:
  build:
    - {{ compiler('cxx') }}
    - cmake >=3.16
    - make
  host:
    - boost-cpp
    - zlib
    - bzip2
  run:
    - boost-cpp
    - zlib
    - bzip2

test:
  commands:
    - reseq2 --version

about:
  home: https://github.com/berntpopp/reseq2
  license: MIT
  license_family: MIT
  license_file: LICENSE
  summary: "ReSeq2 - Realistic Illumina sequencing simulator"
  description: |
    ReSeq2 is a maintained continuation of ReSeq (Schmeing & Robinson, Genome Biology 2021).
    It simulates realistic Illumina paired-end sequencing data by learning error, quality,
    and coverage profiles from real data.
  dev_url: https://github.com/berntpopp/reseq2
  doc_url: https://github.com/berntpopp/reseq2#readme

extra:
  recipe-maintainers:
    - berntpopp
```

- [ ] **Step 2: Commit**

```bash
git add conda/meta.yaml
git commit -m "build: add draft Bioconda recipe for reseq2"
```

Note: The actual Bioconda submission happens as a PR to the `bioconda/bioconda-recipes` repository after the GitHub release is published and the source tarball sha256 is known.

---

### Task 10: Optional -- Send Courtesy FYI to Original Author

This is non-blocking and has no code changes. Do it whenever convenient.

- [ ] **Step 1: Open an issue on schmeing/ReSeq (informational, not asking permission)**

Title: `FYI: Maintained continuation available as ReSeq2`

Body:
```
Hi @schmeing,

I've been maintaining a fork of ReSeq with significant modernization work (C++20, CMake overhaul, CI/CD, expanded tests). Since the upstream appears inactive, I've published it as an independent project:

https://github.com/berntpopp/reseq2

The LICENSE retains your original copyright and the README credits your work. The original Genome Biology paper is cited prominently.

If you have any concerns, please let me know. Thanks for creating ReSeq -- it's great software.

Best,
Bernt Popp
```

---

## Future Work (Not Part of This Plan)

- **Bioconda submission:** After v2.0.0 release, update `conda/meta.yaml` with real sha256 and submit PR to `bioconda/bioconda-recipes`
- **JOSS paper:** Write and submit a JOSS paper describing ReSeq2 improvements (750-1750 words)
- **bio.tools registration:** Register ReSeq2 at https://bio.tools/
- **Archive old fork:** Add a notice to berntpopp/ReSeq README pointing to reseq2, then archive the repo
