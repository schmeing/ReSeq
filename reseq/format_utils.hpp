// Formatting utilities for ReSeq.
//
// Extracted from the vendored 2016-05-15_ROOTPWA/utilities/reportingUtils.hpp
// (originally by Boris Grube, TUM; adapted for ReSeq by Stephan Schmeing, UZH).
// Licensed under GPLv3 — see the original file header for full text.
//
// Namespace mapping: rpwa:: → reseq::format_utils::

#ifndef RESEQ_FORMAT_UTILS_HPP
#define RESEQ_FORMAT_UTILS_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace reseq::format_utils {

//////////////////////////////////////////////////////////////////////////////
// output stream manipulator that prints a value with its maximum precision

template <typename T> class maxPrecisionValue_;

template <typename T> inline maxPrecisionValue_<T> maxPrecision(const T& value) {
    return maxPrecisionValue_<T>(value);
}

// Reserves space so that values will align.
template <typename T> inline maxPrecisionValue_<T> maxPrecisionAlign(const T& value) {
    return maxPrecisionValue_<T>(value, maxPrecisionValue_<T>::ALIGN);
}

// Forces maximum precision for double.
template <typename T> inline maxPrecisionValue_<T> maxPrecisionDouble(const T& value) {
    return maxPrecisionValue_<T>(value, maxPrecisionValue_<T>::DOUBLE);
}

template <typename T> class maxPrecisionValue_ {
  public:
    enum modeEnum { PLAIN, ALIGN, DOUBLE };
    maxPrecisionValue_(const T& value, const modeEnum mode = PLAIN) : _value(value), _mode(mode) {}
    std::ostream& print(std::ostream& out) const {
        const int nmbDigits =
            (_mode != DOUBLE) ? std::numeric_limits<T>::digits10 + 1 : std::numeric_limits<double>::digits10 + 1;
        std::ostringstream s;
        s.precision(nmbDigits);
        s.setf(std::ios_base::scientific, std::ios_base::floatfield);
        s << _value;
        switch (_mode) {
        case ALIGN:
            return out << std::setw(nmbDigits + 7) << s.str(); // space for sign, dot, exponent
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

template <typename T> inline std::ostream& operator<<(std::ostream& out, const maxPrecisionValue_<T>& value) {
    return value.print(out);
}

//////////////////////////////////////////////////////////////////////////////
// indenting

inline void indent(std::ostream& out, const unsigned int offset) {
    for (unsigned int i = 0; i < offset; ++i) {
        out << " ";
    }
}

//////////////////////////////////////////////////////////////////////////////
// simple stream operators for common STL classes

template <typename T1, typename T2> inline std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& pair) {
    return out << "(" << pair.first << ", " << pair.second << ")";
}

template <typename T> inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    if (vec.size() == 0) {
        return out << "{}";
    }
    out << "{";
    for (unsigned int i = 0; i < (vec.size() - 1); ++i) {
        out << "[" << i << "] = " << vec[i] << ", ";
    }
    return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
}

//////////////////////////////////////////////////////////////////////////////
// various utilities

template <typename T> inline unsigned int nmbOfDigits(const T& val) {
    double logVal = 0;
    if (val > 0) {
        logVal = log(val);
    }
    if (val < 0) {
        logVal = log(-val);
    }
    return (unsigned int)(logVal / log(10)) + 1;
}

} // namespace reseq::format_utils

#endif // RESEQ_FORMAT_UTILS_HPP
