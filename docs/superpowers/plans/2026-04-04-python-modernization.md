# Python Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Modernize ReSeq's optional Python scripts and SWIG build so they work on current Python/CMake and are covered by automated tests.

**Architecture:** Keep Python optional behind `RESEQ_BUILD_PYTHON`, modernize the helper scripts in place, and convert the Python build subtree to target-based `Python3` and current `UseSWIG` APIs. Verification runs through repository tests by combining focused Python behavior tests with CTest smoke coverage.

**Tech Stack:** Python 3, argparse, gzip, CMake, UseSWIG, CTest, SWIG, GoogleTest-era repo build conventions

---

### Task 1: Add failing Python tests for script behavior

**Files:**
- Create: `python/tests/test_reseq_prepare_names.py`
- Create: `python/tests/test_plot_data_stats.py`
- Modify: `python/CMakeLists.txt`

- [ ] **Step 1: Write the failing tests**

```python
from pathlib import Path
import gzip
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python" / "reseq-prepare-names.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_fastq(path: Path, header: str) -> None:
    path.write_text(f"{header}\nACGT\n+\n!!!!\n", encoding="utf-8")


def write_fastq_gz(path: Path, header: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"{header}\nTGCA\n+\n####\n")


def test_reseq_prepare_names_rewrites_first_file_headers_for_plain_inputs(tmp_path):
    file1 = tmp_path / "reads_1.fastq"
    file2 = tmp_path / "reads_2.fastq"
    write_fastq(file1, "@INST:1:FCID:2:2104:15343:197393 1:N:0:ATCACG")
    write_fastq(file2, "@INST:1:FCID:2:2104:15343:197393 2:N:0:ATCACG")

    result = run_script(str(file1), str(file2))

    assert result.returncode == 0
    assert result.stdout.splitlines()[0] == "@INST:1:FCID:2:2104:15343:197393"


def test_reseq_prepare_names_supports_gzip_inputs(tmp_path):
    file1 = tmp_path / "reads_1.fastq.gz"
    file2 = tmp_path / "reads_2.fastq.gz"
    write_fastq_gz(file1, "@INST:1:FCID:2:2104:15343:197393 1:N:0:ATCACG")
    write_fastq_gz(file2, "@INST:1:FCID:2:2104:15343:197393 2:N:0:ATCACG")

    result = run_script(str(file1), str(file2))

    assert result.returncode == 0
    assert result.stdout.splitlines()[0] == "@INST:1:FCID:2:2104:15343:197393"
```

```python
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python" / "plotDataStats.py"


def test_plot_data_stats_help_exits_successfully():
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "Usage" in result.stdout or "usage" in result.stdout
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m pytest python/tests/test_reseq_prepare_names.py python/tests/test_plot_data_stats.py -v`
Expected: FAIL because the tests are not yet wired into the repo and at least one script behavior is still on legacy implementation paths.

- [ ] **Step 3: Add test registration hooks**

```cmake
find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_test(
  NAME python_prepare_names_pytest
  COMMAND ${Python3_EXECUTABLE} -m pytest
          ${PROJECT_SOURCE_DIR}/python/tests/test_reseq_prepare_names.py
          ${PROJECT_SOURCE_DIR}/python/tests/test_plot_data_stats.py
)
set_tests_properties(python_prepare_names_pytest PROPERTIES
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
```

- [ ] **Step 4: Run the targeted test command again**

Run: `python3 -m pytest python/tests/test_reseq_prepare_names.py python/tests/test_plot_data_stats.py -v`
Expected: FAIL with concrete assertions against the current script implementations, not import or collection errors.

- [ ] **Step 5: Commit**

```bash
git add python/tests/test_reseq_prepare_names.py python/tests/test_plot_data_stats.py python/CMakeLists.txt
git commit -m "test: add python script smoke coverage"
```

### Task 2: Modernize `reseq-prepare-names.py`

**Files:**
- Modify: `python/reseq-prepare-names.py`
- Test: `python/tests/test_reseq_prepare_names.py`

- [ ] **Step 1: Write one more failing CLI test for argument validation**

```python
def test_reseq_prepare_names_requires_exactly_two_inputs():
    result = run_script()

    assert result.returncode != 0
    assert "usage" in result.stderr.lower()
```

- [ ] **Step 2: Run the targeted test to verify it fails**

Run: `python3 -m pytest python/tests/test_reseq_prepare_names.py::test_reseq_prepare_names_requires_exactly_two_inputs -v`
Expected: FAIL because the current script prints custom usage to stdout and exits via legacy `getopt` handling.

- [ ] **Step 3: Write the minimal implementation**

```python
import argparse
import gzip
from pathlib import Path
import sys


def open_text(path):
    file_path = Path(path)
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt", encoding="utf-8")
    return file_path.open("r", encoding="utf-8")


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite FASTQ names in File1 so paired reads share a common name "
            "without losing tile information."
        )
    )
    parser.add_argument("file1")
    parser.add_argument("file2")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    prepareNames(args.file1, args.file2)
```

- [ ] **Step 4: Run the script tests to verify they pass**

Run: `python3 -m pytest python/tests/test_reseq_prepare_names.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python/reseq-prepare-names.py python/tests/test_reseq_prepare_names.py
git commit -m "feat: modernize reseq prepare names script"
```

### Task 3: Modernize `plotDataStats.py`

**Files:**
- Modify: `python/plotDataStats.py`
- Test: `python/tests/test_plot_data_stats.py`

- [ ] **Step 1: Add a failing import-level compatibility test**

```python
def test_plot_data_stats_compiles_on_modern_python():
    result = subprocess.run(
        [sys.executable, "-m", "py_compile", str(SCRIPT)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
```

- [ ] **Step 2: Run the targeted test to verify it fails**

Run: `python3 -m pytest python/tests/test_plot_data_stats.py::test_plot_data_stats_compiles_on_modern_python -v`
Expected: FAIL on modern Python because `time.clock` has been removed.

- [ ] **Step 3: Write the minimal implementation**

```python
from time import perf_counter
```

```python
start_time = perf_counter()
```

- [ ] **Step 4: Run the plot script tests to verify they pass**

Run: `python3 -m pytest python/tests/test_plot_data_stats.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python/plotDataStats.py python/tests/test_plot_data_stats.py
git commit -m "fix: modernize plot data stats timing"
```

### Task 4: Modernize Python/SWIG CMake integration

**Files:**
- Modify: `python/CMakeLists.txt`
- Modify: `reseq/CMakeLists.txt`

- [ ] **Step 1: Write the failing build/import smoke test registration**

```cmake
add_test(
  NAME python_datastats_import
  COMMAND ${CMAKE_COMMAND}
          -DRESEQ_PYMODS=${PROJECT_BINARY_DIR}/pyMods
          -DPython3_EXECUTABLE=${Python3_EXECUTABLE}
          -P ${PROJECT_SOURCE_DIR}/python/tests/cmake/import_datastats.cmake
)
```

- [ ] **Step 2: Run the Python-enabled configure/build to verify the current path fails or remains on legacy APIs**

Run: `cmake -S . -B build-python -DRESEQ_BUILD_PYTHON=ON`
Expected: current legacy Python/SWIG configuration is still in use and the new smoke test is not yet satisfiable.

- [ ] **Step 3: Write the minimal modern CMake implementation**

```cmake
find_package(SWIG 3 REQUIRED)
include(UseSWIG)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

set(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/pyMods)
set_source_files_properties(DataStats.i PROPERTIES
  CPLUSPLUS ON
  SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE
)

swig_add_library(DataStats
  TYPE MODULE
  LANGUAGE python
  SOURCES DataStats.i
)

target_include_directories(${SWIG_MODULE_DataStats_REAL_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/reseq
  ${PROJECT_SOURCE_DIR}/python
)

target_link_libraries(${SWIG_MODULE_DataStats_REAL_NAME} PRIVATE
  Python3::Module
  reseq_lib
)
```

- [ ] **Step 4: Run the Python-enabled configure/build and import smoke test**

Run: `cmake -S . -B build-python -DRESEQ_BUILD_PYTHON=ON && cmake --build build-python --target _DataStats && ctest --test-dir build-python --output-on-failure -R 'python_'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add python/CMakeLists.txt reseq/CMakeLists.txt python/tests/cmake/import_datastats.cmake
git commit -m "build: modernize python swig integration"
```
