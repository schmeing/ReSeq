# Phase 3: Safety Stabilization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate memory safety hazards, modernize null handling, add exception-safe locking, absorb the vendored ROOTPWA dependency, and convert plain enums to scoped enums — creating the safety foundation required before Phase 4 concurrency work.

**Architecture:** The vendored `2016-05-15_ROOTPWA/reportingUtils.hpp` (404 lines) is split into `reseq/logging.hpp` (logging macros, verbosity, terminal color) and `reseq/format_utils.hpp` (precision printing, stream operators, indenting). All 85 `NULL` literals become `nullptr`. ~23 manual `mutex.lock()`/`unlock()` pairs become `std::scoped_lock` or `std::unique_lock`. Raw `new`/`delete` across Simulator, CoverageStats, DataStats, and test fixtures become `std::unique_ptr`/`std::make_unique`. Three plain enums become `enum class`. One throwing copy constructor becomes `= delete`.

**Tech Stack:** C++20, CMake 3.16+, GoogleTest, std::unique_ptr, std::scoped_lock, std::unique_lock

---

## File Structure

### New Files

| File | Responsibility |
|------|---------------|
| `reseq/logging.hpp` | `NullBuffer`, `NullStream`, `kNullStream`, verbosity globals (`kVerbosityLevel`, `kNoDebugOutput`), terminal color detection, VT100 escape codes, logging macros (`printErr`/`printWarn`/`printInfo`/`printSucc`/`printDebug`), `getClassMethod__`, `log_time` |
| `reseq/format_utils.hpp` | `omanip`, `maxPrecisionValue__`, `maxPrecision`, `maxPrecisionAlign`, `maxPrecisionDouble`, `indent`, stream operators for `std::pair`/`std::vector`, `nmbOfDigits` |

### Modified Files (by sub-phase)

**3a (ROOTPWA absorption):** 18 files that `#include "reportingUtils.hpp"` → include `logging.hpp` and/or `format_utils.hpp`. Delete `2016-05-15_ROOTPWA/` directory. Update `reseq/CMakeLists.txt` to remove `vendored_rootpwa`.

**3b (NULL → nullptr):** 32 files across `reseq/` containing 85 `NULL` literals.

**3c (scoped_lock):** `ProbabilityEstimates.h`, `ProbabilityEstimates.cpp`, `FragmentDistributionStats.cpp`, `CoverageStats.cpp`, `Simulator.cpp`.

**3d (unique_ptr):** `Simulator.h`, `Simulator.cpp`, `CoverageStats.h`, `CoverageStats.cpp`, `DataStats.cpp`, plus 10 test files.

**3e (minor modernizations):** `ErrorStats.h`, `FragmentDistributionStats.h`, `ProbabilityEstimates.h`, `CoverageStats.h`, `Simulator.cpp`.

---

## Task 1: Create `reseq/logging.hpp` (from ROOTPWA)

**Files:**
- Create: `reseq/logging.hpp`

- [ ] **Step 1: Create `reseq/logging.hpp` with logging infrastructure**

```cpp
#ifndef RESEQ_LOGGING_HPP
#define RESEQ_LOGGING_HPP

#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace reseq {
extern uint16_t kVerbosityLevel;
extern bool kNoDebugOutput;

class NullBuffer : public std::streambuf {
  public:
    int overflow(int c) { return c; }
};

class NullStream : public std::ostream {
  public:
    NullStream() : std::ostream(&m_sb_) {}

  private:
    NullBuffer m_sb_;
};

static NullStream kNullStream;
} // namespace reseq

namespace reseq::logging {

//////////////////////////////////////////////////////////////////////////////
// functions for colored output on interactive terminals

// based on code taken from cmake (www.cmake.org)
// Source/kwsys/Terminal.c
inline bool terminalSupportsColor() {
    // terminals running inside emacs do not support escape sequences
    const char* envEmacs = getenv("EMACS");
    if (envEmacs and (*envEmacs == 't'))
        return false;

    // check terminal name
    const std::string colorTermNames[] = {
        "Eterm",
        "ansi",
        "color-xterm",
        "con132x25",
        "con132x30",
        "con132x43",
        "con132x60",
        "con80x25",
        "con80x28",
        "con80x30",
        "con80x43",
        "con80x50",
        "con80x60",
        "cons25",
        "console",
        "cygwin",
        "dtterm",
        "eterm-color",
        "gnome",
        "gnome-256color",
        "konsole",
        "konsole-256color",
        "kterm",
        "linux",
        "msys",
        "linux-c",
        "mach-color",
        "mlterm",
        "putty",
        "rxvt",
        "rxvt-256color",
        "rxvt-cygwin",
        "rxvt-cygwin-native",
        "rxvt-unicode",
        "rxvt-unicode-256color",
        "screen",
        "screen-256color",
        "screen-256color-bce",
        "screen-bce",
        "screen-w",
        "screen.linux",
        "vt100",
        "xterm",
        "xterm-16color",
        "xterm-256color",
        "xterm-88color",
        "xterm-color",
        "xterm-debian",
    };
    const char* envTerm = getenv("TERM");
    if (envTerm) {
        const std::string t = envTerm;
        for (unsigned i = 0; i < std::size(colorTermNames); ++i)
            if (colorTermNames[i] == t)
                return true;
    }
    return false;
}
#ifndef __CINT__
const bool isColorTerminal = terminalSupportsColor();
#endif

inline bool streamIsNotInteractive(const int fileDescriptor) {
    struct stat streamStatus;
    if (fstat(fileDescriptor, &streamStatus) == 0)
        // check whether stream is a regular file
        if (streamStatus.st_mode & S_IFREG)
            return true;
    return false;
}

inline bool stdoutIsColorTerminal() {
    return isColorTerminal and (isatty(STDOUT_FILENO) == 1) and not streamIsNotInteractive(STDOUT_FILENO);
}

inline bool stderrIsColorTerminal() {
    return isColorTerminal and (isatty(STDERR_FILENO) == 1) and not streamIsNotInteractive(STDERR_FILENO);
}

// VT100 escape sequences
enum class Vt100EscapeCode {
    NORMAL = 0,
    BOLD = 1,
    UNDERLINE = 4,
    BLINK = 5,
    INVERSE = 7,
    FG_BLACK = 30,
    FG_RED = 31,
    FG_GREEN = 32,
    FG_YELLOW = 33,
    FG_BLUE = 34,
    FG_MAGENTA = 35,
    FG_CYAN = 36,
    FG_WHITE = 37,
    BG_BLACK = 40,
    BG_RED = 41,
    BG_GREEN = 42,
    BG_YELLOW = 43,
    BG_BLUE = 44,
    BG_MAGENTA = 45,
    BG_CYAN = 46,
    BG_WHITE = 47
};

// ostream manipulator function that inserts VT100 escape sequence into stream
inline std::ostream& vt100SequenceFor(std::ostream& out, const Vt100EscapeCode& vt100Code) {
    bool isVt100 = false;
    if (out.rdbuf() == std::cout.rdbuf())
        isVt100 = stdoutIsColorTerminal();
    else if (out.rdbuf() == std::cerr.rdbuf())
        isVt100 = stderrIsColorTerminal();
    if (isVt100)
        out << "\33[" << static_cast<int>(vt100Code) << "m";
    return out;
}

// general output stream manipulator with one argument
template <typename T> class omanip {
  private:
    typedef std::ostream& (*funcPointer)(std::ostream&, const T&);

  public:
    omanip(funcPointer func, const T& val) : _func(func), _val(val) {}

    friend std::ostream& operator<<(std::ostream& out, const omanip& manip) { return manip._func(out, manip._val); }

  private:
    funcPointer _func;
    T _val;
};

// ostream manipulator for VT100 codes
inline omanip<Vt100EscapeCode> setStreamTo(const Vt100EscapeCode vt100Code) {
    return omanip<Vt100EscapeCode>(&vt100SequenceFor, vt100Code);
}

//////////////////////////////////////////////////////////////////////////////
// macros for printing errors, warnings, and infos

// cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
inline std::string getClassMethod__(std::string prettyFunction) {
    size_t pos = prettyFunction.find("(");
    if (pos == std::string::npos)
        return prettyFunction; // something is not right
    prettyFunction.erase(pos); // cut away signature
    pos = prettyFunction.rfind(" ");
    if (pos == std::string::npos)
        return prettyFunction; // something is not right
    prettyFunction.erase(0, pos + 1); // cut away return type
    return prettyFunction;
}

inline std::tm* log_time() {
    auto t = std::time(nullptr);
    return std::localtime(&t);
}

} // namespace reseq::logging

#ifndef __CINT__
// clang-format off
#define printErr   ((1>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_RED    ) << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "   << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printWarn  ((2>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_YELLOW ) << "??? " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: " << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printSucc  ((3>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_GREEN  ) << "*** " << std::put_time(reseq::logging::log_time(), "%d-%m-%y %H:%M:%S") << ": success: " << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printInfo  ((4>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::BOLD      ) << ">>> " << std::put_time(reseq::logging::log_time(), "%d-%m-%y %H:%M:%S") << ": info: "    << ' ' << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printDebug ((reseq::kNoDebugOutput)? reseq::kNullStream : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_MAGENTA) << "+++ " << reseq::logging::getClassMethod__(__PRETTY_FUNCTION__) << "(): debug: "   << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
// clang-format on
#else
// rootcint crashes in dictionary generation (sigh)
#define printErr   ((1>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << std::flush)
#define printWarn  ((2>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << std::flush)
#define printSucc  ((3>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << std::flush)
#define printInfo  ((4>reseq::kVerbosityLevel)? reseq::kNullStream : std::cerr << std::flush)
#define printDebug ((reseq::kNoDebugOutput)? reseq::kNullStream : std::cerr << std::flush)
#endif

#endif // RESEQ_LOGGING_HPP
```

Key changes from original:
- `rpwa::` namespace → `reseq::logging::`
- `<stdlib.h>` → `<cstdlib>`, `<stdint.h>` → `<cstdint>`
- `enum vt100EscapeCodesEnum` → `enum class Vt100EscapeCode` (requires `static_cast<int>` in escape sequence output)
- `sizeof(array)/sizeof(element)` → `std::size()`
- `std::time(NULL)` → `std::time(nullptr)`
- All macros now use `reseq::logging::` qualification instead of `rpwa::`

- [ ] **Step 2: Build to verify the new header compiles**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build (header not yet included anywhere)

- [ ] **Step 3: Commit**

```bash
git add reseq/logging.hpp
git commit -m "refactor(3a): create reseq/logging.hpp from ROOTPWA reportingUtils

Extract logging infrastructure (NullStream, verbosity globals, terminal
color detection, VT100 escape codes, logging macros) from vendored
2016-05-15_ROOTPWA/utilities/reportingUtils.hpp into reseq/logging.hpp.

Changes during absorption:
- rpwa:: namespace -> reseq::logging::
- <stdlib.h> -> <cstdlib>, <stdint.h> -> <cstdint>
- enum vt100EscapeCodesEnum -> enum class Vt100EscapeCode
- sizeof(array)/sizeof(element) -> std::size()
- std::time(NULL) -> std::time(nullptr)"
```

---

## Task 2: Create `reseq/format_utils.hpp` (from ROOTPWA)

**Files:**
- Create: `reseq/format_utils.hpp`

- [ ] **Step 1: Create `reseq/format_utils.hpp` with formatting utilities**

```cpp
#ifndef RESEQ_FORMAT_UTILS_HPP
#define RESEQ_FORMAT_UTILS_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace reseq::format_utils {

//////////////////////////////////////////////////////////////////////////////
// output stream manipulator that prints a value with its maximum precision

template <typename T> class maxPrecisionValue__;

template <typename T> inline maxPrecisionValue__<T> maxPrecision(const T& value) {
    return maxPrecisionValue__<T>(value);
}

// output stream manipulator that prints a value with its maximum precision
// in addition manipulator reserves space so that values will align
template <typename T> inline maxPrecisionValue__<T> maxPrecisionAlign(const T& value) {
    return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::ALIGN);
}

// output stream manipulator that prints a value with maximum precision for double
template <typename T> inline maxPrecisionValue__<T> maxPrecisionDouble(const T& value) {
    return maxPrecisionValue__<T>(value, maxPrecisionValue__<T>::DOUBLE);
}

// general helper class that encapsulates a value of type T
template <typename T> class maxPrecisionValue__ {
  public:
    enum modeEnum { PLAIN, ALIGN, DOUBLE }; // forces precision for double
    maxPrecisionValue__(const T& value, const modeEnum mode = PLAIN) : _value(value), _mode(mode) {}
    std::ostream& print(std::ostream& out) const {
        const int nmbDigits =
            (_mode != DOUBLE) ? std::numeric_limits<T>::digits10 + 1 : std::numeric_limits<double>::digits10 + 1;
        std::ostringstream s;
        s.precision(nmbDigits);
        s.setf(std::ios_base::scientific, std::ios_base::floatfield);
        s << _value;
        switch (_mode) {
        case ALIGN:
            return out << std::setw(nmbDigits + 7) << s.str(); // make space for sign, dot, and exponent
        case PLAIN:
        case DOUBLE:
        default:
            return out << s.str();
        }
    }

  private:
    const T& _value;
    modeEnum _mode;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& out, const maxPrecisionValue__<T>& value) {
    return value.print(out);
}

//////////////////////////////////////////////////////////////////////////////
// indenting

inline void indent(std::ostream& out, const unsigned int offset) {
    for (unsigned int i = 0; i < offset; ++i)
        out << " ";
}

//////////////////////////////////////////////////////////////////////////////
// simple stream operators for some STL classes

template <typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& pair) {
    return out << "(" << pair.first << ", " << pair.second << ")";
}

template <typename T> inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    if (vec.size() == 0) {
        return out << "{}";
    }
    out << "{";
    for (unsigned int i = 0; i < (vec.size() - 1); ++i)
        out << "[" << i << "] = " << vec[i] << ", ";
    return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
}

//////////////////////////////////////////////////////////////////////////////
// various stuff

template <typename T> inline unsigned int nmbOfDigits(const T& val) {
    double logVal = 0;
    if (val > 0)
        logVal = log(val);
    if (val < 0)
        logVal = log(-val);
    return (unsigned int)(logVal / log(10)) + 1;
}

} // namespace reseq::format_utils

#endif // RESEQ_FORMAT_UTILS_HPP
```

Key change: `rpwa::` namespace → `reseq::format_utils::`

- [ ] **Step 2: Build to verify the new header compiles**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build

- [ ] **Step 3: Commit**

```bash
git add reseq/format_utils.hpp
git commit -m "refactor(3a): create reseq/format_utils.hpp from ROOTPWA reportingUtils

Extract formatting utilities (maxPrecision, indent, stream operators for
pair/vector, nmbOfDigits) from vendored reportingUtils.hpp into
reseq/format_utils.hpp. Namespace rpwa:: -> reseq::format_utils::."
```

---

## Task 3: Update all include sites to use new headers

**Files:**
- Modify: `reseq/utilities.hpp:15`
- Modify: `reseq/Vect.hpp:14`
- Modify: `reseq/SeqQualityStats.hpp:8`
- Modify: `reseq/BasicTestClass.hpp:10`
- Modify: `reseq/CoverageStats.h:13`
- Modify: `reseq/test_main.cpp:10`
- Modify: `reseq/main.cpp:23`
- Modify: `reseq/Simulator.cpp:45`
- Modify: `reseq/Reference.cpp:37`
- Modify: `reseq/FragmentDistributionStats.cpp:58`
- Modify: `reseq/DataStats.cpp:37`
- Modify: `reseq/TileStats.cpp:18`
- Modify: `reseq/Surrounding.cpp:12`
- Modify: `reseq/QualityStats.cpp:9`
- Modify: `reseq/FragmentDuplicationStats.cpp:9`
- Modify: `reseq/ProbabilityEstimates.cpp:40`
- Modify: `reseq/ErrorStats.cpp:8`
- Modify: `reseq/AdapterStats.cpp:27`

- [ ] **Step 1: Determine which header each file needs**

Most files only use logging macros (`printErr`, `printWarn`, `printInfo`, `printSucc`, `printDebug`) → they need `logging.hpp` only.

Files that also use formatting utilities (`maxPrecision`, `maxPrecisionAlign`, `maxPrecisionDouble`, `indent`, `nmbOfDigits`, or stream operators for `pair`/`vector`):

Search for `rpwa::maxPrecision`, `rpwa::indent`, `rpwa::nmbOfDigits` usage to determine which files need `format_utils.hpp`. The format utilities are used via ADL through the `rpwa::` `operator<<` overloads for `std::pair` and `std::vector`, plus explicit calls.

Run: `grep -rn 'maxPrecision\|rpwa::indent\|nmbOfDigits' reseq/ --include='*.cpp' --include='*.h' --include='*.hpp' | grep -v reportingUtils`

Based on codebase analysis, the mapping is:

**Only `logging.hpp`:** `Simulator.cpp`, `Reference.cpp`, `FragmentDistributionStats.cpp`, `DataStats.cpp`, `TileStats.cpp`, `Surrounding.cpp`, `QualityStats.cpp`, `FragmentDuplicationStats.cpp`, `ErrorStats.cpp`, `AdapterStats.cpp`, `test_main.cpp`, `main.cpp`, `BasicTestClass.hpp`, `SeqQualityStats.hpp`

**Both `logging.hpp` and `format_utils.hpp`:** `utilities.hpp` (uses stream operators), `Vect.hpp` (uses stream operators via printErr), `CoverageStats.h` (uses logging), `ProbabilityEstimates.cpp` (uses maxPrecision in output)

- [ ] **Step 2: Replace all `#include "reportingUtils.hpp"` with appropriate new headers**

In each of the 18 files, replace:
```cpp
#include "reportingUtils.hpp"
```
with one of:
```cpp
#include "logging.hpp"
```
or:
```cpp
#include "format_utils.hpp"
#include "logging.hpp"
```

For files that use `rpwa::` qualified names (format utility functions), also replace the namespace qualification. Search for any explicit `rpwa::` usage in the codebase outside of `reportingUtils.hpp` — the exploration found **zero** direct `rpwa::` references in `reseq/` files, because the macros in `reportingUtils.hpp` already expand with full qualification. The ADL-found `operator<<` overloads for `std::pair`/`std::vector` also don't require explicit `rpwa::` qualification.

**Important:** The logging macros (`printErr`, etc.) reference `reseq::logging::setStreamTo` and `reseq::logging::Vt100EscapeCode` — these are defined inside `logging.hpp`, so any file using logging macros only needs `#include "logging.hpp"`.

Replace in all 18 files. For example, in `reseq/utilities.hpp:15`:
```cpp
// Before:
#include "reportingUtils.hpp"
// After:
#include "format_utils.hpp"
#include "logging.hpp"
```

In `reseq/Simulator.cpp:45`:
```cpp
// Before:
#include "reportingUtils.hpp"
// After:
#include "logging.hpp"
```

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build with no errors.

If there are compile errors about missing symbols from `format_utils.hpp`, add that include to the affected file.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/utilities.hpp reseq/Vect.hpp reseq/SeqQualityStats.hpp reseq/BasicTestClass.hpp \
        reseq/CoverageStats.h reseq/test_main.cpp reseq/main.cpp reseq/Simulator.cpp \
        reseq/Reference.cpp reseq/FragmentDistributionStats.cpp reseq/DataStats.cpp \
        reseq/TileStats.cpp reseq/Surrounding.cpp reseq/QualityStats.cpp \
        reseq/FragmentDuplicationStats.cpp reseq/ProbabilityEstimates.cpp \
        reseq/ErrorStats.cpp reseq/AdapterStats.cpp
git commit -m "refactor(3a): replace reportingUtils.hpp includes with logging/format_utils

Update all 18 include sites to use the new reseq/logging.hpp and
reseq/format_utils.hpp headers instead of the vendored ROOTPWA
reportingUtils.hpp."
```

---

## Task 4: Delete vendored ROOTPWA and remove CMake target

**Files:**
- Delete: `2016-05-15_ROOTPWA/` (entire directory)
- Modify: `CMakeLists.txt` (top-level, remove `vendored_rootpwa` target)
- Modify: `reseq/CMakeLists.txt:31` (remove `vendored_rootpwa` from link)

- [ ] **Step 1: Identify the vendored_rootpwa CMake target definition**

Run: `grep -rn 'vendored_rootpwa\|rootpwa' CMakeLists.txt */CMakeLists.txt`

This will show where the INTERFACE target is defined and where it's linked.

- [ ] **Step 2: Remove `vendored_rootpwa` from `reseq/CMakeLists.txt` link line**

In `reseq/CMakeLists.txt:31`, change:
```cmake
# Before:
target_link_libraries(reseq_lib PUBLIC
  skewer_matrix nlopt vendored_seqan vendored_rootpwa GTest::gtest
  ${Boost_LIBRARIES} ZLIB::ZLIB BZip2::BZip2 rt
)
# After:
target_link_libraries(reseq_lib PUBLIC
  skewer_matrix nlopt vendored_seqan GTest::gtest
  ${Boost_LIBRARIES} ZLIB::ZLIB BZip2::BZip2 rt
)
```

- [ ] **Step 3: Remove `vendored_rootpwa` INTERFACE target from top-level CMakeLists.txt**

Find and delete the `add_library(vendored_rootpwa INTERFACE)` block and its associated `target_include_directories` in the top-level `CMakeLists.txt`.

- [ ] **Step 4: Delete the vendored directory**

```bash
rm -rf 2016-05-15_ROOTPWA/
```

- [ ] **Step 5: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make clean && make build 2>&1 | tail -30`
Expected: Clean build succeeds with no references to ROOTPWA.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor(3a): delete vendored ROOTPWA, remove CMake target

Remove 2016-05-15_ROOTPWA/ directory and vendored_rootpwa INTERFACE
target. All functionality has been absorbed into reseq/logging.hpp and
reseq/format_utils.hpp."
```

---

## Task 5: NULL → nullptr (production headers)

**Files:**
- Modify: `reseq/Simulator.h` (lines 120, 131)
- Modify: `reseq/CoverageStats.h` (lines 93, 385)
- Modify: `reseq/ProbabilityEstimates.h` (lines 265, 1392, 1409, 1429)
- Modify: `reseq/FragmentDistributionStats.h` (lines 94, 95)
- Modify: `reseq/TileStats.h` (lines 45, 84, 87)
- Modify: `reseq/DataStats.h` (1 occurrence)
- Modify: `reseq/AdapterStats.h` (1 occurrence)

- [ ] **Step 1: Replace all NULL with nullptr in production headers**

This is a mechanical find-and-replace within each file. In each file, replace every occurrence of `NULL` with `nullptr`. Examples:

`reseq/Simulator.h:120`:
```cpp
// Before:
            : id_(id), start_pos_(start_pos), finished_(false), next_block_(NULL), partner_block_(partner_block),
// After:
            : id_(id), start_pos_(start_pos), finished_(false), next_block_(nullptr), partner_block_(partner_block),
```

`reseq/Simulator.h:131`:
```cpp
// Before:
            : ref_seq_id_(ref_seq_id), first_block_(NULL), last_block_(NULL), next_unit_(NULL) {}
// After:
            : ref_seq_id_(ref_seq_id), first_block_(nullptr), last_block_(nullptr), next_unit_(nullptr) {}
```

`reseq/CoverageStats.h:93`:
```cpp
// Before:
            : sequence_id_(seq_id), start_pos_(start_pos), previous_block_(prev_block), next_block_(NULL),
// After:
            : sequence_id_(seq_id), start_pos_(start_pos), previous_block_(prev_block), next_block_(nullptr),
```

`reseq/ProbabilityEstimates.h:265`:
```cpp
// Before:
        const std::array<std::vector<uintMatrixIndex>, N>* dim_indices_new_count = NULL)
// After:
        const std::array<std::vector<uintMatrixIndex>, N>* dim_indices_new_count = nullptr)
```

`reseq/ProbabilityEstimates.h` (lines 1392, 1409, 1429):
```cpp
// Before:
            margins.at(N) = {NULL, true};
// After:
            margins.at(N) = {nullptr, true};
```

`reseq/FragmentDistributionStats.h` (lines 94-95):
```cpp
// Before:
    static constexpr const char* kParameterInfoFile = NULL;  //"maxlike_fit.csv";
    static constexpr const char* kDispersionInfoFile = NULL; //"dispersion_fit.csv";
// After:
    static constexpr const char* kParameterInfoFile = nullptr;  //"maxlike_fit.csv";
    static constexpr const char* kDispersionInfoFile = nullptr; //"dispersion_fit.csv";
```

`reseq/TileStats.h` (lines 45, 84, 87) — replace `= NULL)` with `= nullptr)`.

`reseq/DataStats.h` and `reseq/AdapterStats.h` — replace their `NULL` occurrences similarly.

- [ ] **Step 2: Build to verify**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.h reseq/CoverageStats.h reseq/ProbabilityEstimates.h \
        reseq/FragmentDistributionStats.h reseq/TileStats.h reseq/DataStats.h reseq/AdapterStats.h
git commit -m "fix(3b): replace NULL with nullptr in production headers

Mechanical replacement of NULL -> nullptr across 7 header files.
Covers default parameters, member initializer lists, and constexpr
pointer constants."
```

---

## Task 6: NULL → nullptr (production .cpp files)

**Files:**
- Modify: `reseq/Simulator.cpp` (13 occurrences: lines 549, 683, 932, 936, 1246, 1317, 1328, 2511, 2546, 2737, 3005, 3013, 3014)
- Modify: `reseq/CoverageStats.cpp` (6 occurrences: lines 458, 895, 995, 996, 1071, 1072)
- Modify: `reseq/ProbabilityEstimates.cpp` (3 occurrences: lines 936, 953, 974)
- Modify: `reseq/DataStats.cpp` (2 occurrences: lines 931, 1035)
- Modify: `reseq/main.cpp` (3 occurrences: lines 538, 825, 1074)

- [ ] **Step 1: Replace all NULL with nullptr in each .cpp file**

Mechanical replacement. In each file, replace every `NULL` with `nullptr`. Examples:

`reseq/Simulator.cpp:2737`:
```cpp
// Before:
Simulator::Simulator() : written_records_(0), last_unit_(NULL), deletion_buffer_(0), rdist_zero_to_one_(0, 1) {
// After:
Simulator::Simulator() : written_records_(0), last_unit_(nullptr), deletion_buffer_(0), rdist_zero_to_one_(0, 1) {
```

`reseq/Simulator.cpp:932`:
```cpp
// Before:
    auto block = new SimBlock(first_block_id, 0, NULL, block_seed_gen_());
// After:
    auto block = new SimBlock(first_block_id, 0, nullptr, block_seed_gen_());
```

`reseq/CoverageStats.cpp:1071-1072`:
```cpp
// Before:
            first_block_ = NULL;
            last_block_ = NULL;
// After:
            first_block_ = nullptr;
            last_block_ = nullptr;
```

`reseq/main.cpp:538`:
```cpp
// Before:
        DataStats real_data_stats(NULL);
// After:
        DataStats real_data_stats(nullptr);
```

Continue for all remaining occurrences in each file.

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.cpp reseq/CoverageStats.cpp reseq/ProbabilityEstimates.cpp \
        reseq/DataStats.cpp reseq/main.cpp
git commit -m "fix(3b): replace NULL with nullptr in production .cpp files

Mechanical replacement of NULL -> nullptr across 5 implementation files
(27 occurrences). Covers constructor args, member assignments, pointer
comparisons, and function call arguments."
```

---

## Task 7: NULL → nullptr (test files)

**Files:**
- Modify: `reseq/CoverageStatsTest.cpp` (10 occurrences)
- Modify: `reseq/ProbabilityEstimatesTest.cpp` (7 occurrences)
- Modify: `reseq/SimulatorTest.cpp` (6 occurrences)
- Modify: `reseq/DataStatsTest.cpp` (2 occurrences)
- Modify: `reseq/DataStatsTest.h` (1 occurrence)
- Modify: `reseq/AdapterStatsTest.h` (1 occurrence)
- Modify: `reseq/AdapterStatsTest.cpp` (1 occurrence)
- Modify: `reseq/ErrorStatsTest.h` (1 occurrence)
- Modify: `reseq/ErrorStatsTest.cpp` (1 occurrence)
- Modify: `reseq/FragmentDistributionStatsTest.cpp` (1 occurrence)
- Modify: `reseq/FragmentDuplicationStatsTest.h` (1 occurrence)
- Modify: `reseq/FragmentDuplicationStatsTest.cpp` (1 occurrence)
- Modify: `reseq/QualityStatsTest.h` (1 occurrence)
- Modify: `reseq/QualityStatsTest.cpp` (1 occurrence)
- Modify: `reseq/TileStatsTest.h` (1 occurrence)
- Modify: `reseq/TileStatsTest.cpp` (1 occurrence)

- [ ] **Step 1: Replace all NULL with nullptr in test files**

Mechanical replacement in each file. Notable patterns:

`reseq/CoverageStatsTest.cpp` — test assertions:
```cpp
// Before:
    EXPECT_TRUE(NULL == test.first_block_);
// After:
    EXPECT_TRUE(nullptr == test.first_block_);
```

`reseq/ProbabilityEstimatesTest.cpp` — constructor calls:
```cpp
// Before:
    DataStats stats(NULL);
// After:
    DataStats stats(nullptr);
```

Test fixture destructors:
```cpp
// Before:
    test_ = NULL;
// After:
    test_ = nullptr;
```

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 3: Verify no NULL remains in reseq/**

Run: `grep -rn '\bNULL\b' reseq/ --include='*.cpp' --include='*.h' --include='*.hpp' | grep -v '//.*NULL'`
Expected: Zero matches (no remaining NULL literals in reseq/ source).

- [ ] **Step 4: Commit**

```bash
git add reseq/*Test*.cpp reseq/*Test*.h
git commit -m "fix(3b): replace NULL with nullptr in test files

Mechanical replacement of NULL -> nullptr across 16 test files
(~35 occurrences). Covers test fixture cleanup, constructor arguments,
and pointer assertions."
```

---

## Task 8: scoped_lock — ProbabilityEstimates print_mutex

**Files:**
- Modify: `reseq/ProbabilityEstimates.h` (7 lock/unlock pairs at lines 319-322, 949-951, 1096-1100, 1126-1129, 1139-1143, 1173-1176, 1197-1206)
- Modify: `reseq/ProbabilityEstimates.cpp` (2 lock/unlock pairs at lines 425-438, 446-458)

- [ ] **Step 1: Replace lock/unlock pairs with scoped_lock in ProbabilityEstimates.h**

Each pair follows the pattern: `print_mutex.lock(); ... print_mutex.unlock();`. Replace with a block-scoped `std::scoped_lock`.

Example at line 319-322:
```cpp
// Before:
                    print_mutex.lock();
                    printErr << "dim2_[" << n << "][" << i / dim_size_.at(dim_b) << "][" << i % dim_size_.at(dim_b)
                             << "] is NaN for " << descriptor << std::endl;
                    print_mutex.unlock();
// After:
                    {
                        std::scoped_lock lock(print_mutex);
                        printErr << "dim2_[" << n << "][" << i / dim_size_.at(dim_b) << "][" << i % dim_size_.at(dim_b)
                                 << "] is NaN for " << descriptor << std::endl;
                    }
```

Example at line 949-951:
```cpp
// Before:
            print_mutex.lock();
            printErr << "precision_ is NaN for " << descriptor << std::endl;
            print_mutex.unlock();
// After:
            {
                std::scoped_lock lock(print_mutex);
                printErr << "precision_ is NaN for " << descriptor << std::endl;
            }
```

Example at line 1197-1206 (multi-statement block):
```cpp
// Before:
                    print_mutex.lock();
                    if (precision_ > precision_aim) {
                        printWarn << descriptor << " did not reach precision aim: " << precision_ * 100 << "%"
                                  << std::endl;
                    } else {
                        printInfo << "Finished iterative proportional fitting for " << descriptor << " after step "
                                  << steps_ << " with precision " << precision_ * 100 << "%" << std::endl;
                    }
                    print_mutex.unlock();
// After:
                    {
                        std::scoped_lock lock(print_mutex);
                        if (precision_ > precision_aim) {
                            printWarn << descriptor << " did not reach precision aim: " << precision_ * 100 << "%"
                                      << std::endl;
                        } else {
                            printInfo << "Finished iterative proportional fitting for " << descriptor << " after step "
                                      << steps_ << " with precision " << precision_ * 100 << "%" << std::endl;
                        }
                    }
```

Apply the same pattern to all 7 pairs in the header.

- [ ] **Step 2: Replace lock/unlock pairs in ProbabilityEstimates.cpp**

At lines 425-438 and 446-458, the lock/unlock guards multi-line error reporting blocks. Wrap each in a scoped block:

```cpp
// Before (line 425):
                print_mutex.lock();
                printErr << "The dimension " << dim_a << " is inconsistent in margin " << n << ":" << std::endl;
                if (0 < kVerbosityLevel) {
                    for (auto element : dim_control_a) {
                        std::cerr << element << ' ';
                    }
                    std::cerr << std::endl;
                    for (auto element : dim_control.at(dim_a)) {
                        std::cerr << element << ' ';
                    }
                    std::cerr << std::endl;
                }
                print_mutex.unlock();
// After:
                {
                    std::scoped_lock lock(print_mutex);
                    printErr << "The dimension " << dim_a << " is inconsistent in margin " << n << ":" << std::endl;
                    if (0 < kVerbosityLevel) {
                        for (auto element : dim_control_a) {
                            std::cerr << element << ' ';
                        }
                        std::cerr << std::endl;
                        for (auto element : dim_control.at(dim_a)) {
                            std::cerr << element << ' ';
                        }
                        std::cerr << std::endl;
                    }
                }
```

Same pattern for the dim_b check at line 446.

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/ProbabilityEstimates.h reseq/ProbabilityEstimates.cpp
git commit -m "fix(3c): replace manual lock/unlock with scoped_lock in ProbabilityEstimates

Convert 9 manual print_mutex.lock()/unlock() pairs to std::scoped_lock
for exception safety. All are logging-only critical sections."
```

---

## Task 9: scoped_lock — FragmentDistributionStats result_mutex

**Files:**
- Modify: `reseq/FragmentDistributionStats.cpp` (1 lock/unlock pair at lines 3114-3123)

- [ ] **Step 1: Replace lock/unlock with scoped_lock**

```cpp
// Before (line 3114-3123):
    result_mutex.lock();
    for (auto cov_group = max_bias.size(); cov_group--;) {
        for (auto frag_len = max_bias.at(cov_group).size(); frag_len--;) {
            SetToMax(max_bias.at(cov_group).at(frag_len).at(0), tmp_max_bias.at(cov_group).at(frag_len));
        }
    }
    for (auto frag_len = norm.size(); frag_len--;) {
        norm.at(frag_len) += tmp_norm.at(frag_len);
    }
    result_mutex.unlock();
// After:
    {
        std::scoped_lock lock(result_mutex);
        for (auto cov_group = max_bias.size(); cov_group--;) {
            for (auto frag_len = max_bias.at(cov_group).size(); frag_len--;) {
                SetToMax(max_bias.at(cov_group).at(frag_len).at(0), tmp_max_bias.at(cov_group).at(frag_len));
            }
        }
        for (auto frag_len = norm.size(); frag_len--;) {
            norm.at(frag_len) += tmp_norm.at(frag_len);
        }
    }
```

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add reseq/FragmentDistributionStats.cpp
git commit -m "fix(3c): replace manual lock/unlock with scoped_lock in FragmentDistributionStats

Convert result_mutex lock/unlock pair in bias calculation lambda to
std::scoped_lock for exception safety."
```

---

## Task 10: scoped_lock — Simulator print_mutex_ (simple logging locks)

**Files:**
- Modify: `reseq/Simulator.cpp` (lines 142-145, 174-177, 209-211)

- [ ] **Step 1: Replace print_mutex_ lock/unlock pairs with scoped_lock**

These are all simple logging-only critical sections.

Line 142-145:
```cpp
// Before:
        print_mutex_.lock();
        printErr << "Could not write records " << written_records_ + 1 << " to "
                 << (written_records_ + length(*old_output_ids)) << ": " << e.what() << std::endl;
        print_mutex_.unlock();
// After:
        {
            std::scoped_lock lock(print_mutex_);
            printErr << "Could not write records " << written_records_ + 1 << " to "
                     << (written_records_ + length(*old_output_ids)) << ": " << e.what() << std::endl;
        }
```

Line 174-177:
```cpp
// Before:
        print_mutex_.lock();
        printInfo << "Generated " << written_records_ << " read pairs ("
                  << static_cast<uintPercentPrint>(Percent(written_records_, total_pairs_)) << "%)." << std::endl;
        print_mutex_.unlock();
// After:
        {
            std::scoped_lock lock(print_mutex_);
            printInfo << "Generated " << written_records_ << " read pairs ("
                      << static_cast<uintPercentPrint>(Percent(written_records_, total_pairs_)) << "%)." << std::endl;
        }
```

Line 209-211:
```cpp
// Before:
        print_mutex_.lock();
        printInfo << "Generated " << written_records_ << " reads." << std::endl;
        print_mutex_.unlock();
// After:
        {
            std::scoped_lock lock(print_mutex_);
            printInfo << "Generated " << written_records_ << " reads." << std::endl;
        }
```

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "fix(3c): replace print_mutex_ lock/unlock with scoped_lock in Simulator

Convert 3 manual print_mutex_ lock/unlock pairs in Simulator to
std::scoped_lock for exception safety."
```

---

## Task 11: scoped_lock — Simulator block_creation_mutex_

**Files:**
- Modify: `reseq/Simulator.cpp` (lines 2678-2693)

- [ ] **Step 1: Replace block_creation_mutex_ lock/unlock with scoped_lock**

```cpp
// Before (line 2678-2693):
        self.block_creation_mutex_.lock();

        if (atEnd(org_seq_reader)) {
            keep_running = false;
        } else {
            rgen.seed(self.block_seed_gen_());
            rdist.Reset();
            cur_block = self.read_blocks_++;

            readRecords(input_ids, input_seqs, org_seq_reader, self.kBatchSizeErrorModelOnly);
            if (0 == length(input_ids)) {
                keep_running = false;
            }
        }

        self.block_creation_mutex_.unlock();
// After:
        {
            std::scoped_lock lock(self.block_creation_mutex_);

            if (atEnd(org_seq_reader)) {
                keep_running = false;
            } else {
                rgen.seed(self.block_seed_gen_());
                rdist.Reset();
                cur_block = self.read_blocks_++;

                readRecords(input_ids, input_seqs, org_seq_reader, self.kBatchSizeErrorModelOnly);
                if (0 == length(input_ids)) {
                    keep_running = false;
                }
            }
        }
```

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "fix(3c): replace block_creation_mutex_ lock/unlock with scoped_lock in Simulator

Convert manual block_creation_mutex_ lock/unlock to std::scoped_lock
for exception safety in seqToIllumina read processing loop."
```

---

## Task 12: scoped_lock — CoverageStats try_lock patterns

**Files:**
- Modify: `reseq/CoverageStats.cpp` (lines 449-472, 1036-1043)

**Note:** The `try_lock` patterns in CoverageStats are **not** simple lock/unlock pairs — they use `mutex.try_lock()` which returns immediately if the lock is unavailable. These cannot be directly converted to `std::scoped_lock`. Instead, use `std::unique_lock` with `std::try_to_lock`.

- [ ] **Step 1: Convert reuse_mutex_ try_lock pattern (lines 449-472)**

```cpp
// Before:
CoverageStats::CoverageBlock* CoverageStats::CreateBlock(uintRefSeqId seq_id, uintSeqLen start_pos) {
    CoverageBlock* new_block;
    if (reuse_mutex_.try_lock()) {
        if (reusable_blocks_.size()) {
            new_block = reusable_blocks_.back();
            reusable_blocks_.pop_back();
            reuse_mutex_.unlock();

            new_block->sequence_id_ = seq_id;
            new_block->start_pos_ = start_pos;
            new_block->previous_block_ = last_block_;
            new_block->next_block_ = nullptr;
            new_block->coverage_.clear();
            new_block->previous_coverage_.clear();
            new_block->reads_.clear();
            new_block->scheduled_for_processing_.clear();
            new_block->processed_ = false;
        } else {
            reuse_mutex_.unlock();
            new_block = new CoverageBlock(seq_id, start_pos, last_block_);
            new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
        }
    } else {
        new_block = new CoverageBlock(seq_id, start_pos, last_block_);
        new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
    }
// After:
CoverageStats::CoverageBlock* CoverageStats::CreateBlock(uintRefSeqId seq_id, uintSeqLen start_pos) {
    CoverageBlock* new_block;
    if (std::unique_lock lock(reuse_mutex_, std::try_to_lock); lock.owns_lock()) {
        if (reusable_blocks_.size()) {
            new_block = reusable_blocks_.back();
            reusable_blocks_.pop_back();
            lock.unlock();

            new_block->sequence_id_ = seq_id;
            new_block->start_pos_ = start_pos;
            new_block->previous_block_ = last_block_;
            new_block->next_block_ = nullptr;
            new_block->coverage_.clear();
            new_block->previous_coverage_.clear();
            new_block->reads_.clear();
            new_block->scheduled_for_processing_.clear();
            new_block->processed_ = false;
        } else {
            lock.unlock();
            new_block = new CoverageBlock(seq_id, start_pos, last_block_);
            new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
        }
    } else {
        new_block = new CoverageBlock(seq_id, start_pos, last_block_);
        new_block->previous_coverage_.reserve(maximum_read_length_on_reference_);
    }
```

Note: The early `lock.unlock()` calls are kept because the code deliberately releases the lock before doing further work. The `unique_lock` destructor will be a no-op for the already-unlocked lock.

- [ ] **Step 2: Convert variant_loading_mutex_ try_lock pattern (lines 1036-1043)**

```cpp
// Before:
        if (variant_loading_mutex_.try_lock()) {
            if (!reference.ReadVariantPositions((*last_block_).sequence_id_ + 2)) {
                variant_loading_mutex_.unlock();
                return false;
            }

            variant_loading_mutex_.unlock();
        }
// After:
        if (std::unique_lock lock(variant_loading_mutex_, std::try_to_lock); lock.owns_lock()) {
            if (!reference.ReadVariantPositions((*last_block_).sequence_id_ + 2)) {
                return false;
            }
        }
```

Here the `unique_lock` destructor handles unlock on both the success and error paths, simplifying the code.

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/CoverageStats.cpp
git commit -m "fix(3c): replace try_lock patterns with unique_lock in CoverageStats

Convert reuse_mutex_ and variant_loading_mutex_ try_lock/unlock patterns
to std::unique_lock with std::try_to_lock for exception safety.
Simplifies error-path unlocking in PreLoadVariants."
```

---

## Task 13: scoped_lock — Simulator try_lock patterns

**Files:**
- Modify: `reseq/Simulator.cpp` (lines 1338-1347, 1355-1361)

- [ ] **Step 1: Convert var_read_mutex_ try_lock pattern (lines 1338-1347)**

```cpp
// Before:
        if (var_read_mutex_.try_lock()) {
            uintSeqLen max_del_shift = 0;
            if (!ref.ReadVariants(
                    max_del_shift, current_unit_->ref_seq_id_ + 2,
                    2 * stats.MaxReadLenOnReference())) {
                var_read_mutex_.unlock();
                return false;
            }
            RequestBufferSize(max_del_shift);
            var_read_mutex_.unlock();
// After:
        if (std::unique_lock lock(var_read_mutex_, std::try_to_lock); lock.owns_lock()) {
            uintSeqLen max_del_shift = 0;
            if (!ref.ReadVariants(
                    max_del_shift, current_unit_->ref_seq_id_ + 2,
                    2 * stats.MaxReadLenOnReference())) {
                return false;
            }
            RequestBufferSize(max_del_shift);
```

- [ ] **Step 2: Convert methylation_read_mutex_ try_lock pattern (lines 1355-1361)**

```cpp
// Before:
        if (methylation_read_mutex_.try_lock()) {
            if (!ref.ReadMethylation(current_unit_->ref_seq_id_ + 2)) {
                methylation_read_mutex_.unlock();
                return false;
            }
            methylation_read_mutex_.unlock();
// After:
        if (std::unique_lock lock(methylation_read_mutex_, std::try_to_lock); lock.owns_lock()) {
            if (!ref.ReadMethylation(current_unit_->ref_seq_id_ + 2)) {
                return false;
            }
```

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "fix(3c): replace try_lock patterns with unique_lock in Simulator

Convert var_read_mutex_ and methylation_read_mutex_ try_lock/unlock
patterns to std::unique_lock with std::try_to_lock. Eliminates manual
unlock on error paths."
```

---

## Task 14: scoped_lock — Simulator complex locking (output_mutex_, flush_mutex_)

**Files:**
- Modify: `reseq/Simulator.cpp` (lines 152-170, 221-232, 2989)

**Note:** These locking patterns are more complex and cannot all be trivially converted. The `Flush()` method has a carefully ordered multi-mutex protocol (flush_mutex_[0] → flush_mutex_[1] with overlap). The `Output()` method acquires `output_mutex_` and conditionally calls `Flush()` which releases it. These patterns are **deferred to Phase 4** where the entire output pipeline is redesigned with condition variables.

For now, only convert the simple `print_mutex_` cases (already done in Task 10) and the final cleanup lock:

- [ ] **Step 1: Convert output_mutex_ cleanup lock (line 2989)**

This lock is held during the final flush and cleanup. It cannot be scoped_lock because `Flush()` releases it internally. Leave this as-is — it will be redesigned in Phase 4.

**Decision: Skip this task.** The remaining Simulator locking patterns (`output_mutex_`, `flush_mutex_`) involve complex multi-mutex ordering and conditional release that would require redesigning the output pipeline. This is exactly what Phase 4 addresses. Mark as intentionally deferred.

- [ ] **Step 2: Document the deferral**

Add a comment at `Simulator.cpp:152` (Flush method):

```cpp
// NOTE: output_mutex_ and flush_mutex_ use complex multi-mutex ordering
// that cannot be safely converted to scoped_lock without redesigning the
// output pipeline. Deferred to Phase 4 (concurrency modernization).
```

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "docs(3c): document deferred Simulator output mutex conversion

The output_mutex_/flush_mutex_ multi-mutex ordering protocol in
Simulator::Flush() requires Phase 4's output pipeline redesign before
it can use RAII locking. Added explanatory comment."
```

---

## Task 15: unique_ptr — Test fixture new/delete

**Files:**
- Modify: `reseq/ErrorStatsTest.cpp` (new line 5, delete line 10)
- Modify: `reseq/ErrorStatsTest.h` (member type)
- Modify: `reseq/AdapterStatsTest.cpp` (new line 18, delete line 23)
- Modify: `reseq/AdapterStatsTest.h` (member type)
- Modify: `reseq/FragmentDistributionStatsTest.cpp` (new lines 71/1275, delete line 91)
- Modify: `reseq/FragmentDistributionStatsTest.h` (member type)
- Modify: `reseq/CoverageStatsTest.cpp` (new line 8, delete line 13)
- Modify: `reseq/CoverageStatsTest.h` (member type if exists)
- Modify: `reseq/DataStatsTest.cpp` (new line 36, delete line 41)
- Modify: `reseq/DataStatsTest.h` (member type)
- Modify: `reseq/QualityStatsTest.cpp` (new line 10, delete line 15)
- Modify: `reseq/QualityStatsTest.h` (member type)
- Modify: `reseq/SimulatorTest.cpp` (new line 39, delete line 44)
- Modify: `reseq/SimulatorTest.h` (member type)
- Modify: `reseq/FragmentDuplicationStatsTest.cpp` (new line 14, delete line 20)
- Modify: `reseq/FragmentDuplicationStatsTest.h` (member type)
- Modify: `reseq/TileStatsTest.cpp` (new line 19, delete line 24)
- Modify: `reseq/TileStatsTest.h` (member type)

**Strategy:** Each test fixture has a `T* test_` member, a `CreateTestObject()` that does `test_ = new T(...)`, and a destructor that does `delete test_; test_ = NULL;`. Convert to `std::unique_ptr<T> test_` with `test_ = std::make_unique<T>(...)` and remove the explicit delete.

- [ ] **Step 1: Read one test header/cpp pair to understand the pattern**

Run: Read `reseq/ErrorStatsTest.h` and `reseq/ErrorStatsTest.cpp` to see the exact pattern.

- [ ] **Step 2: Convert each test fixture**

For each test fixture, make two changes:

**In the header (.h):** Change the member from raw pointer to unique_ptr:
```cpp
// Before:
    ErrorStats* test_;
// After:
    std::unique_ptr<ErrorStats> test_;
```

Add `#include <memory>` if not already present.

**In the .cpp:** Change `new` to `make_unique`, remove `delete`, remove null assignment in destructor:
```cpp
// Before (constructor/create):
    ASSERT_TRUE(test_ = new ErrorStats()) << "Could not allocate memory for ErrorStats object\n";
// After:
    test_ = std::make_unique<ErrorStats>();

// Before (destructor):
    delete test_;
    test_ = nullptr;
// After:
    test_.reset();
```

Note: `ASSERT_TRUE(test_ = new ...)` used the assignment-in-condition pattern. With `make_unique`, the allocation either succeeds or throws, so the ASSERT is unnecessary. Replace with a plain assignment.

Apply this pattern to all 10 test fixtures.

- [ ] **Step 3: Update test code that accesses test_ via pointer**

Any test code that uses `test_->Method()` continues to work because `unique_ptr` supports `->`. Any code that passes `test_` as a raw pointer needs `.get()`:

Search for patterns like `SomeFunction(test_)` (passing raw pointer) — these need `SomeFunction(test_.get())`.

- [ ] **Step 4: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add reseq/*Test*.h reseq/*Test*.cpp
git commit -m "fix(3d): convert test fixture raw new/delete to unique_ptr

Replace raw new/delete pattern in 10 test fixture classes with
std::unique_ptr and std::make_unique. Eliminates leak risk in test
teardown."
```

---

## Task 16: unique_ptr — DataStats FullRecord allocation

**Files:**
- Modify: `reseq/DataStats.cpp` (lines 936, 940, 1058, 1059 for new/delete of FullRecord)
- Modify: `reseq/CoverageStats.h` (line 87: `std::vector<FullRecord*> reads_` → `std::vector<FullRecord*>` stays raw — see note)
- Modify: `reseq/CoverageStats.cpp` (line 208: `delete record`)

**Analysis:** `FullRecord` objects are allocated in `DataStats::CollectStats()` and stored in `CoverageBlock::reads_` vectors. They are deleted either when the record is unpaired (immediate delete in `CoverageStats::EvalRecord`) or when the block is processed. The ownership is: DataStats allocates → CoverageBlock holds → CoverageStats deletes.

**Strategy:** Use `std::unique_ptr<FullRecord>` for the initial allocation. Transfer ownership into the `reads_` vector via `std::move`. Change `reads_` to `std::vector<std::unique_ptr<FullRecord>>`.

- [ ] **Step 1: Change CoverageBlock::reads_ to hold unique_ptr**

In `reseq/CoverageStats.h:87`:
```cpp
// Before:
        std::vector<FullRecord*> reads_;
// After:
        std::vector<std::unique_ptr<FullRecord>> reads_;
```

Add `#include <memory>` to `CoverageStats.h` if not present.

- [ ] **Step 2: Update all reads_ access sites**

Search for all uses of `reads_` in CoverageStats.cpp and DataStats.cpp. Key changes:

In `DataStats.cpp` where records are created and pushed:
```cpp
// Before:
    record = new CoverageStats::FullRecord;
    ...
    block->reads_.push_back(record);
// After:
    auto record = std::make_unique<CoverageStats::FullRecord>();
    ...
    block->reads_.push_back(std::move(record));
```

In `CoverageStats.cpp` where records are iterated:
```cpp
// Before:
    for (auto rec : first_block->reads_) {
        EvalRead(rec, first_block, reference, qualities, errors, phred_quality_offset);
    }
// After:
    for (auto& rec : first_block->reads_) {
        EvalRead(rec.get(), first_block, reference, qualities, errors, phred_quality_offset);
    }
```

In `CoverageStats::EvalRecord` where paired record is deleted:
```cpp
// Before:
    delete record;
// After:
    // unique_ptr handles deallocation — just let it go out of scope or reset
```

**Important:** Check all `EvalRead` call sites to understand whether they take `FullRecord*` — if so, pass `.get()`. The function signatures stay as raw pointer parameters (they don't own the record).

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/CoverageStats.h reseq/CoverageStats.cpp reseq/DataStats.cpp
git commit -m "fix(3d): convert FullRecord allocation to unique_ptr

Change CoverageBlock::reads_ from vector<FullRecord*> to
vector<unique_ptr<FullRecord>>. DataStats creates records with
make_unique and transfers ownership via std::move. Eliminates manual
delete tracking across DataStats/CoverageStats boundary."
```

---

## Task 17: unique_ptr — Simulator StringSet buffers

**Files:**
- Modify: `reseq/Simulator.h` (lines 301-303: `StringSet<>*` members)
- Modify: `reseq/Simulator.cpp` (lines 124/128/132 FlushCopyValues, 180-183 delete, 2826-2830 Initialize, 3075-3079 Initialize)

**Strategy:** Convert the three `std::array<StringSet<>*, 2>` members to `std::array<std::unique_ptr<StringSet<>>, 2>`. Update `FlushCopyValues` to use `std::move` + `std::make_unique`. Remove explicit `delete` calls.

- [ ] **Step 1: Change Simulator.h member types**

```cpp
// Before:
    std::array<seqan::StringSet<seqan::CharString>*, 2> output_ids_;
    std::array<seqan::StringSet<seqan::Dna5String>*, 2> output_seqs_;
    std::array<seqan::StringSet<seqan::CharString>*, 2> output_quals_;
// After:
    std::array<std::unique_ptr<seqan::StringSet<seqan::CharString>>, 2> output_ids_;
    std::array<std::unique_ptr<seqan::StringSet<seqan::Dna5String>>, 2> output_seqs_;
    std::array<std::unique_ptr<seqan::StringSet<seqan::CharString>>, 2> output_quals_;
```

Add `#include <memory>` to `Simulator.h`.

- [ ] **Step 2: Update FlushCopyValues signature and body**

The function signature takes `StringSet<>*&` (reference to raw pointer) for the old values. Since the caller (`Flush()`) needs raw pointers to pass to `FlushWriteValues` (which doesn't own them), change FlushCopyValues to return the old pointers:

```cpp
// Before:
void Simulator::FlushCopyValues(uintTempSeq template_segment, StringSet<CharString>*& old_output_ids,
                                StringSet<Dna5String>*& old_output_seqs, StringSet<CharString>*& old_output_quals) {
    old_output_ids = output_ids_.at(template_segment);
    output_ids_.at(template_segment) = new StringSet<CharString>;
    reserve(*output_ids_.at(template_segment), kBatchSize, Exact());
    // ... same for seqs and quals
}
// After:
void Simulator::FlushCopyValues(uintTempSeq template_segment,
                                std::unique_ptr<StringSet<CharString>>& old_output_ids,
                                std::unique_ptr<StringSet<Dna5String>>& old_output_seqs,
                                std::unique_ptr<StringSet<CharString>>& old_output_quals) {
    old_output_ids = std::move(output_ids_.at(template_segment));
    output_ids_.at(template_segment) = std::make_unique<StringSet<CharString>>();
    reserve(*output_ids_.at(template_segment), kBatchSize, Exact());

    old_output_seqs = std::move(output_seqs_.at(template_segment));
    output_seqs_.at(template_segment) = std::make_unique<StringSet<Dna5String>>();
    reserve(*output_seqs_.at(template_segment), kBatchSize, Exact());

    old_output_quals = std::move(output_quals_.at(template_segment));
    output_quals_.at(template_segment) = std::make_unique<StringSet<CharString>>();
    reserve(*output_quals_.at(template_segment), kBatchSize, Exact());
}
```

- [ ] **Step 3: Update Flush() to use unique_ptr**

```cpp
// Before:
bool Simulator::Flush() {
    array<StringSet<CharString>*, 2> old_output_ids;
    array<StringSet<Dna5String>*, 2> old_output_seqs;
    array<StringSet<CharString>*, 2> old_output_quals;
    ...
    for (uintTempSeq template_segment = 2; template_segment--;) {
        delete old_output_ids.at(template_segment);
        delete old_output_seqs.at(template_segment);
        delete old_output_quals.at(template_segment);
    }
// After:
bool Simulator::Flush() {
    std::array<std::unique_ptr<StringSet<CharString>>, 2> old_output_ids;
    std::array<std::unique_ptr<StringSet<Dna5String>>, 2> old_output_seqs;
    std::array<std::unique_ptr<StringSet<CharString>>, 2> old_output_quals;
    ...
    // Remove explicit delete loop — unique_ptr handles cleanup
```

Update FlushWriteValues calls to pass `.get()`:
```cpp
    bool success = FlushWriteValues(0, old_output_ids.at(0).get(), old_output_seqs.at(0).get(), old_output_quals.at(0).get());
```

- [ ] **Step 4: Update Initialize() allocation sites**

In `Simulator.cpp` around lines 2826 and 3075, replace `new` with `make_unique`:
```cpp
// Before:
    output_ids_.at(i) = new StringSet<CharString>;
// After:
    output_ids_.at(i) = std::make_unique<StringSet<CharString>>();
```

- [ ] **Step 5: Update all dereference sites**

All code that accesses `*output_ids_.at(template_segment)` continues to work because `unique_ptr` supports `operator*`. Code that uses `output_ids_.at(template_segment)` as a raw pointer (e.g., `appendValue(*output_ids_.at(template_segment), ...)`) also works unchanged via `operator*`.

- [ ] **Step 6: Update header declarations to match new signatures**

In `reseq/Simulator.h`, update the `FlushCopyValues` declaration:
```cpp
// Before:
    void FlushCopyValues(uintTempSeq template_segment, seqan::StringSet<seqan::CharString>*& old_output_ids,
                         seqan::StringSet<seqan::Dna5String>*& old_output_seqs,
                         seqan::StringSet<seqan::CharString>*& old_output_quals);
// After:
    void FlushCopyValues(uintTempSeq template_segment,
                         std::unique_ptr<seqan::StringSet<seqan::CharString>>& old_output_ids,
                         std::unique_ptr<seqan::StringSet<seqan::Dna5String>>& old_output_seqs,
                         std::unique_ptr<seqan::StringSet<seqan::CharString>>& old_output_quals);
```

- [ ] **Step 7: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add reseq/Simulator.h reseq/Simulator.cpp
git commit -m "fix(3d): convert Simulator StringSet buffers to unique_ptr

Replace raw StringSet<>* members with unique_ptr. FlushCopyValues uses
std::move to transfer ownership. Eliminates manual delete in Flush()
and potential leaks in Initialize()."
```

---

## Task 18: unique_ptr — CoverageStats CoverageBlock (partial, prepare for Phase 4)

**Files:**
- Modify: `reseq/CoverageStats.h` (reusable_blocks_ type, CoverageBlock creation)
- Modify: `reseq/CoverageStats.cpp` (CreateBlock, PrepareForReading, Finalize)

**Note:** The spec says Phase 3d must convert CoverageBlock to `std::vector<std::unique_ptr<>>` as a **prerequisite for Phase 4**. However, CoverageBlock uses `std::atomic<CoverageBlock*> next_block_` and raw `previous_block_` in a doubly-linked list traversed by concurrent threads. Full container normalization (replacing the linked list with a vector) is Phase 4a work. Phase 3d converts the **ownership** model: the reusable_blocks_ pool and the initial allocation/deletion patterns.

- [ ] **Step 1: Convert reusable_blocks_ to unique_ptr**

In `reseq/CoverageStats.h`, find the declaration:
```cpp
// Before:
    std::vector<CoverageBlock*> reusable_blocks_;
// After:
    std::vector<std::unique_ptr<CoverageBlock>> reusable_blocks_;
```

- [ ] **Step 2: Update Finalize() to remove explicit block deletion**

In `reseq/CoverageStats.cpp`, the Finalize method deletes blocks in two places:

Line 1052-1053 (reusable blocks cleanup):
```cpp
// Before:
    for (auto block : reusable_blocks_) {
        delete block;
    }
// After:
    // unique_ptr handles cleanup automatically
```
Just call `reusable_blocks_.clear()` (which destroys the unique_ptrs).

Lines 1082-1083 and 1092 (linked list traversal cleanup): These delete blocks during traversal. Since the linked list nodes are not yet owned by a container (they're raw-allocated and linked), we need to wrap them in unique_ptr at deletion time or keep raw delete for now.

**Decision:** The linked list traversal in Finalize (lines 1074-1092) deletes nodes as it walks `next_block_`. Since the nodes are connected via `std::atomic<CoverageBlock*>`, converting them to unique_ptr ownership requires the Phase 4a container redesign. For Phase 3d, convert only the `reusable_blocks_` pool. Leave the linked list node allocation/deletion as-is with a comment noting Phase 4 dependency.

- [ ] **Step 3: Update CreateBlock() for reusable_blocks_ unique_ptr**

In `reseq/CoverageStats.cpp:449-472`, when popping from reusable_blocks_:
```cpp
// Before:
            new_block = reusable_blocks_.back();
            reusable_blocks_.pop_back();
// After:
            new_block = reusable_blocks_.back().release();
            reusable_blocks_.pop_back();
```

When blocks are returned to the pool (search for `reusable_blocks_.push_back`):
```cpp
// Before:
    reusable_blocks_.push_back(block);
// After:
    reusable_blocks_.push_back(std::unique_ptr<CoverageBlock>(block));
```

- [ ] **Step 4: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add reseq/CoverageStats.h reseq/CoverageStats.cpp
git commit -m "fix(3d): convert CoverageStats reusable_blocks_ to unique_ptr

Convert reusable block pool from vector<CoverageBlock*> to
vector<unique_ptr<CoverageBlock>>. Linked list node ownership
(next_block_/previous_block_) deferred to Phase 4a container
redesign."
```

---

## Task 19: unique_ptr — Simulator SimBlock/SimUnit (partial, prepare for Phase 4)

**Files:**
- Modify: `reseq/Simulator.cpp` (CreateUnit, CreateBlock/GetNextBlock, Finalize cleanup)

**Note:** Like CoverageBlock, SimBlock/SimUnit use intrusive linked lists with atomic pointers. Full conversion to `vector<unique_ptr<>>` is Phase 4a. Phase 3d focuses on making ownership explicit where possible without changing the data structure.

The spec says: "Phase 3d must convert these to `std::vector<std::unique_ptr<>>` or similar owning containers, eliminating the intrusive linked-list patterns entirely."

However, this is a large structural change that interleaves with concurrency (atomic pointers, lock-free traversal). The safe approach for Phase 3d is to document the current ownership model and prepare the groundwork.

**Decision:** The SimBlock/SimUnit linked list is deeply intertwined with atomic operations and concurrent traversal. Converting to vector<unique_ptr<>> requires the Phase 4a container redesign. Phase 3d adds ownership documentation comments. The actual conversion is Phase 4a.

- [ ] **Step 1: Add ownership documentation to Simulator.h**

Add comments to SimBlock and SimUnit documenting the ownership model:

```cpp
    // SimBlock ownership: Allocated in CreateUnit() and CreateBlock(). Blocks form a
    // singly-linked list via atomic next_block_. partner_block_ links forward↔reverse
    // strand blocks. Deleted in GetNextBlock() (incremental cleanup) and Finalize()
    // (final cleanup). Phase 4 will convert to vector<unique_ptr<SimBlock>>.
    struct SimBlock {
```

```cpp
    // SimUnit ownership: Allocated in CreateUnit(). Units form a singly-linked list
    // via next_unit_. Each unit owns a range of SimBlocks (first_block_ to last_block_).
    // Deleted in GetNextBlock() and Finalize(). Phase 4 will convert to
    // vector<unique_ptr<SimUnit>>.
    struct SimUnit {
```

- [ ] **Step 2: Commit**

```bash
git add reseq/Simulator.h
git commit -m "docs(3d): document SimBlock/SimUnit ownership model

Add ownership documentation comments to SimBlock and SimUnit structs.
Actual conversion to unique_ptr containers deferred to Phase 4a
(data structure redesign required due to atomic pointer traversal)."
```

---

## Task 20: enum class — ErrorStats::InDelDef

**Files:**
- Modify: `reseq/ErrorStats.h:16`
- Modify: All files that reference `kNoInDel`, `kDeletion`, `kInsertionA`, etc.

- [ ] **Step 1: Convert enum to enum class**

In `reseq/ErrorStats.h:16`:
```cpp
// Before:
    enum InDelDef { kNoInDel, kDeletion, kInsertionA, kInsertionC, kInsertionG, kInsertionT, kInsertionN };
// After:
    enum class InDelDef { kNoInDel, kDeletion, kInsertionA, kInsertionC, kInsertionG, kInsertionT, kInsertionN };
```

- [ ] **Step 2: Update all usage sites to qualify enum values**

Search for all uses of `kNoInDel`, `kDeletion`, `kInsertionA`, etc. across the codebase:

Run: `grep -rn 'kNoInDel\|kDeletion\|kInsertionA\|kInsertionC\|kInsertionG\|kInsertionT\|kInsertionN' reseq/ --include='*.cpp' --include='*.h' --include='*.hpp'`

Each unqualified use must be prefixed with `ErrorStats::InDelDef::` or just `InDelDef::` if within the ErrorStats class. For example:
```cpp
// Before:
    if (indel == kNoInDel) {
// After:
    if (indel == InDelDef::kNoInDel) {
```

Or if outside the class:
```cpp
// Before:
    ErrorStats::kNoInDel
// After:
    ErrorStats::InDelDef::kNoInDel
```

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/ErrorStats.h reseq/ErrorStats.cpp reseq/Simulator.cpp reseq/QualityStats.cpp
git commit -m "refactor(3e): convert ErrorStats::InDelDef to enum class

Scoped enum prevents implicit integer conversion and namespace
pollution. All usage sites updated with qualified names."
```

---

## Task 21: enum class — RefSeqBiasSimulation

**Files:**
- Modify: `reseq/FragmentDistributionStats.h:25`
- Modify: All files that reference `kKeep`, `kNo`, `kDraw`, `kFile`, `kError`

- [ ] **Step 1: Convert enum to enum class**

In `reseq/FragmentDistributionStats.h:25`:
```cpp
// Before:
enum RefSeqBiasSimulation { kKeep, kNo, kDraw, kFile, kError };
// After:
enum class RefSeqBiasSimulation { kKeep, kNo, kDraw, kFile, kError };
```

- [ ] **Step 2: Update all usage sites**

Search for all uses of the enum values. These are likely in `main.cpp` (command-line parsing), `Simulator.cpp`, and `FragmentDistributionStats.cpp`:

Run: `grep -rn '\bkKeep\b\|\bkNo\b\|\bkDraw\b\|\bkFile\b\|\bkError\b' reseq/ --include='*.cpp' --include='*.h'`

**Note:** `kNo` and `kError` are very common words — filter results to only those in the context of `RefSeqBiasSimulation`. Each must become `RefSeqBiasSimulation::kKeep`, etc.

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/FragmentDistributionStats.h reseq/FragmentDistributionStats.cpp \
        reseq/Simulator.cpp reseq/main.cpp
git commit -m "refactor(3e): convert RefSeqBiasSimulation to enum class

Scoped enum prevents namespace pollution for generic names like kKeep,
kNo, kError. All usage sites updated with qualified names."
```

---

## Task 22: enum class — ProbabilityEstimates::IPFDataSelector

**Files:**
- Modify: `reseq/ProbabilityEstimates.h:1303-1310`
- Modify: `reseq/ProbabilityEstimates.cpp` (usage sites)

- [ ] **Step 1: Convert enum to enum class**

In `reseq/ProbabilityEstimates.h:1303`:
```cpp
// Before:
    enum IPFDataSelector {
        kIPFQuality,
        kIPFSequenceQuality,
        kIPFBaseCall,
        kIPFDominantError,
        kIPFErrorRate,
        kIPFInDels
    };
// After:
    enum class IPFDataSelector {
        kIPFQuality,
        kIPFSequenceQuality,
        kIPFBaseCall,
        kIPFDominantError,
        kIPFErrorRate,
        kIPFInDels
    };
```

- [ ] **Step 2: Update usage sites**

Search for `kIPFQuality`, `kIPFSequenceQuality`, etc.:

Run: `grep -rn 'kIPF' reseq/ --include='*.cpp' --include='*.h'`

All are within ProbabilityEstimates, so qualify with `IPFDataSelector::`:
```cpp
// Before:
    params.selected_data = kIPFQuality;
// After:
    params.selected_data = IPFDataSelector::kIPFQuality;
```

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -30`
Expected: Successful build.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add reseq/ProbabilityEstimates.h reseq/ProbabilityEstimates.cpp
git commit -m "refactor(3e): convert ProbabilityEstimates::IPFDataSelector to enum class

Scoped enum for IPF data selection. All usage sites within
ProbabilityEstimates updated with qualified names."
```

---

## Task 23: Delete CoveragePosition throwing copy constructor

**Files:**
- Modify: `reseq/CoverageStats.h:58-60`

- [ ] **Step 1: Replace throwing copy constructor with = delete**

In `reseq/CoverageStats.h:58-60`:
```cpp
// Before:
        CoveragePosition(const CoveragePosition& UNUSED(right)) {
            throw std::runtime_error("This function should never be called.");
        }
// After:
        CoveragePosition(const CoveragePosition&) = delete;
```

This is safer because it catches misuse at compile time rather than runtime.

- [ ] **Step 2: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make build 2>&1 | tail -20`
Expected: Successful build. If there's a compile error, it means something was actually calling this copy constructor — investigate.

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add reseq/CoverageStats.h
git commit -m "refactor(3e): delete CoveragePosition copy constructor instead of throwing

Replace runtime-throwing copy constructor with = delete for compile-time
detection of accidental copies."
```

---

## Task 24: Document Simulator failure state (line ~1295)

**Files:**
- Modify: `reseq/Simulator.cpp` (around line 1285-1330)

The spec mentions "make failure state explicit (don't continue through failed block creation)" at line 1295. Looking at the code, after `CreateBlock` sets up a block, if `simulation_error_` is true, it falls through to the cleanup section (line 1302+) which deletes blocks. The `else` at line 1327 sets `current_unit_ = nullptr` and returns false.

- [ ] **Step 1: Review the actual issue**

Read `Simulator.cpp:1285-1330` carefully. The concern is that when `simulation_error_` is true at line 1285, the code skips systematic error setup but still enters the block cleanup loop. This is actually the intended behavior — cleanup happens even on error. The issue is subtle: if CreateBlock returns false due to completion (line 1257) vs error (line 1328), the caller can't distinguish them.

**Decision:** This is a design improvement that touches the simulation control flow and is better addressed as part of Phase 5c (Simulator decomposition). Add a comment documenting the current behavior.

- [ ] **Step 2: Add clarifying comment**

```cpp
    if (!simulation_error_) {
        // Set systematic errors for new block
        ...
    } else {
        // simulation_error_ was set by another thread — skip systematic error setup
        // and signal the caller to stop. Block cleanup happens in Finalize().
        current_unit_ = nullptr;
        return false;
    }
```

- [ ] **Step 3: Commit**

```bash
git add reseq/Simulator.cpp
git commit -m "docs(3e): clarify Simulator::CreateBlock failure state handling

Document the behavior when simulation_error_ is set during block
creation. Structural improvement deferred to Phase 5c decomposition."
```

---

## Task 25: Run format check and fix any violations

**Files:**
- Potentially any file modified in this phase

- [ ] **Step 1: Run format check**

Run: `cd /home/bernt-popp/development/ReSeq && make format-check 2>&1`
Expected: May show formatting violations in new/modified files.

- [ ] **Step 2: Fix formatting**

Run: `cd /home/bernt-popp/development/ReSeq && make format`

- [ ] **Step 3: Build and run tests**

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "style: fix formatting violations in Phase 3 changes"
```

---

## Task 26: Final verification gate

- [ ] **Step 1: Clean rebuild**

Run: `cd /home/bernt-popp/development/ReSeq && make clean && make build 2>&1 | tail -30`
Expected: Clean build with no warnings related to Phase 3 changes.

- [ ] **Step 2: Run all tests**

Run: `cd /home/bernt-popp/development/ReSeq && make test 2>&1`
Expected: All tests pass (57 unit tests + regression tests).

- [ ] **Step 3: Verify no NULL remains**

Run: `grep -rn '\bNULL\b' reseq/ --include='*.cpp' --include='*.h' --include='*.hpp' | grep -v '//.*NULL' | wc -l`
Expected: 0

- [ ] **Step 4: Verify no reportingUtils.hpp includes remain**

Run: `grep -rn 'reportingUtils' reseq/ --include='*.cpp' --include='*.h' --include='*.hpp'`
Expected: 0 matches

- [ ] **Step 5: Verify ROOTPWA directory is gone**

Run: `ls -d 2016-05-15_ROOTPWA/ 2>/dev/null && echo "STILL EXISTS" || echo "REMOVED"`
Expected: REMOVED

- [ ] **Step 6: Run format check**

Run: `cd /home/bernt-popp/development/ReSeq && make format-check`
Expected: No violations.

- [ ] **Step 7: Run pre-commit hooks**

Run: `cd /home/bernt-popp/development/ReSeq && pre-commit run --all-files`
Expected: All hooks pass.
