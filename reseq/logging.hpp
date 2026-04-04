// Logging and terminal-color infrastructure for ReSeq.
//
// Extracted from the vendored 2016-05-15_ROOTPWA/utilities/reportingUtils.hpp
// (originally by Boris Grube, TUM; adapted for ReSeq by Stephan Schmeing, UZH).
// Licensed under GPLv3 — see the original file header for full text.
//
// Namespace mapping: rpwa:: → reseq::logging::

#ifndef RESEQ_LOGGING_HPP
#define RESEQ_LOGGING_HPP

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace reseq {

extern std::atomic<uint16_t> kVerbosityLevel;
extern bool kNoDebugOutput;

class NullBuffer : public std::streambuf {
  public:
    int overflow(int c) override { return c; }
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
//
// Based on code taken from CMake (www.cmake.org), Source/kwsys/Terminal.c

inline bool terminalSupportsColor() {
    // Terminals running inside emacs do not support escape sequences.
    const char* envEmacs = getenv("EMACS");
    if (envEmacs and (*envEmacs == 't')) {
        return false;
    }

    // Check terminal name.
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
        for (const auto& colorTermName : colorTermNames) {
            if (colorTermName == t) {
                return true;
            }
        }
    }
    return false;
}

#ifndef __CINT__
const bool isColorTerminal = terminalSupportsColor();
#endif

inline bool streamIsNotInteractive(const int fileDescriptor) {
    struct stat streamStatus;
    if (fstat(fileDescriptor, &streamStatus) == 0) {
        if (streamStatus.st_mode & S_IFREG) {
            return true;
        }
    }
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
    BG_WHITE = 47,
};

//////////////////////////////////////////////////////////////////////////////
// general output stream manipulator with one argument

template <typename T> class omanip {
  private:
    using funcPointer = std::ostream& (*)(std::ostream&, const T&);

  public:
    omanip(funcPointer func, const T& val) : _func(func), _val(val) {}

    friend std::ostream& operator<<(std::ostream& out, const omanip& manip) { return manip._func(out, manip._val); }

  private:
    funcPointer _func;
    T _val;
};

// ostream manipulator function that inserts VT100 escape sequence into stream
inline std::ostream& vt100SequenceFor(std::ostream& out, const Vt100EscapeCode& vt100Code) {
    bool isVt100 = false;
    if (out.rdbuf() == std::cout.rdbuf()) {
        isVt100 = stdoutIsColorTerminal();
    } else if (out.rdbuf() == std::cerr.rdbuf()) {
        isVt100 = stderrIsColorTerminal();
    }
    if (isVt100) {
        out << "\33[" << static_cast<int>(vt100Code) << "m";
    }
    return out;
}

// ostream manipulator for VT100 codes
inline omanip<Vt100EscapeCode> setStreamTo(const Vt100EscapeCode vt100Code) {
    return omanip<Vt100EscapeCode>(&vt100SequenceFor, vt100Code);
}

// Cuts out block "className::methodName" from __PRETTY_FUNCTION__ output.
inline std::string getClassMethod_(std::string prettyFunction) {
    size_t pos = prettyFunction.find('(');
    if (pos == std::string::npos) {
        return prettyFunction;
    }
    prettyFunction.erase(pos);
    pos = prettyFunction.rfind(' ');
    if (pos == std::string::npos) {
        return prettyFunction;
    }
    prettyFunction.erase(0, pos + 1);
    return prettyFunction;
}

inline std::tm* log_time() {
    auto t = std::time(nullptr);
    return std::localtime(&t);
}

} // namespace reseq::logging

//////////////////////////////////////////////////////////////////////////////
// Logging macros
//
// These intentionally live outside any namespace so that they remain globally
// accessible exactly as before.  They use reseq::logging:: qualifications
// instead of the old rpwa:: qualifications.

#ifndef __CINT__
#define printErr                                                                                                       \
    ((1 > reseq::kVerbosityLevel)                                                                                      \
         ? reseq::kNullStream                                                                                          \
         : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_RED) << "!!! "                 \
                     << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: error: "                     \
                     << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printWarn                                                                                                      \
    ((2 > reseq::kVerbosityLevel)                                                                                      \
         ? reseq::kNullStream                                                                                          \
         : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_YELLOW) << "??? "              \
                     << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: warning: "                   \
                     << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printSucc                                                                                                      \
    ((3 > reseq::kVerbosityLevel)                                                                                      \
         ? reseq::kNullStream                                                                                          \
         : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_GREEN) << "*** "               \
                     << std::put_time(reseq::logging::log_time(), "%d-%m-%y %H:%M:%S") << ": success: "                \
                     << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printInfo                                                                                                      \
    ((4 > reseq::kVerbosityLevel)                                                                                      \
         ? reseq::kNullStream                                                                                          \
         : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::BOLD) << ">>> "                   \
                     << std::put_time(reseq::logging::log_time(), "%d-%m-%y %H:%M:%S") << ": info: " << ' '            \
                     << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#define printDebug                                                                                                     \
    ((reseq::kNoDebugOutput)                                                                                           \
         ? reseq::kNullStream                                                                                          \
         : std::cerr << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::FG_MAGENTA) << "+++ "             \
                     << reseq::logging::getClassMethod__(__PRETTY_FUNCTION__) << "(): debug: "                         \
                     << reseq::logging::setStreamTo(reseq::logging::Vt100EscapeCode::NORMAL) << std::flush)
#else
// rootcint crashes in dictionary generation
#define printErr ((1 > reseq::kVerbosityLevel) ? reseq::kNullStream : std::cerr << std::flush)
#define printWarn ((2 > reseq::kVerbosityLevel) ? reseq::kNullStream : std::cerr << std::flush)
#define printSucc ((3 > reseq::kVerbosityLevel) ? reseq::kNullStream : std::cerr << std::flush)
#define printInfo ((4 > reseq::kVerbosityLevel) ? reseq::kNullStream : std::cerr << std::flush)
#define printDebug ((reseq::kNoDebugOutput) ? reseq::kNullStream : std::cerr << std::flush)
#endif

#endif // RESEQ_LOGGING_HPP
