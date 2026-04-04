#ifndef UTILITIES_H
#define UTILITIES_H

#include <array>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <sys/stat.h>

#include <boost/filesystem.hpp>

#include "format_utils.hpp"
#include "logging.hpp"
#include "types.hpp"

#include "CMakeConfig.h"

namespace reseq {

#ifndef SWIG // This part is not needed for the python plotting and swig can't handle the seqan stuff
namespace utilities {

// Functions needed in classes
template <typename T> struct AtDummy {
    static T dummy_;
};

template <typename T> inline const T& at(const seqan::StringSet<T>& set, size_t n) {
    if (length(set) > n) {
        return set[n];
    } else {
        printErr << "Called position " << n << ". StringSet size is only " << length(set) << "." << std::endl;
        throw std::out_of_range("Accessing StringSet after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T> inline T& at(seqan::StringSet<T>& set, size_t n) {
    if (length(set) > n) {
        return set[n];
    } else {
        printErr << "Called position " << n << ". StringSet size is only " << length(set) << "." << std::endl;
        throw std::out_of_range("Accessing StringSet after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T> inline const T& at(const seqan::String<T>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". String size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing String after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T> inline T& at(seqan::String<T>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". String size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing String after last position");
        return AtDummy<T>::dummy_;
    }
}

// Do not return references for ModifiedString types, as they generate temporary values and it does not work and we
// don't want to change the underlying string in those cases anyways
template <typename T, typename U>
inline T at(const seqan::ModifiedString<const seqan::String<T>, U>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U> inline T at(const seqan::ModifiedString<seqan::String<T>, U>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U> inline T at(seqan::ModifiedString<seqan::String<T>, U>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U, typename V>
inline T at(const seqan::ModifiedString<seqan::ModifiedString<const seqan::String<T>, U>, V>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U, typename V>
inline T at(const seqan::ModifiedString<seqan::ModifiedString<seqan::String<T>, U>, V>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U, typename V>
inline T at(seqan::ModifiedString<seqan::ModifiedString<seqan::String<T>, U>, V>& string, size_t n) {
    if (length(string) > n) {
        return string[n];
    } else {
        printErr << "Called position " << n << ". ModifiedString size is only " << length(string) << "." << std::endl;
        throw std::out_of_range("Accessing ModifiedString after last position");
        return AtDummy<T>::dummy_;
    }
}

template <typename T, typename U> inline void SetToMax(T& current_max, const U potential_max) {
    if (potential_max > current_max) {
        current_max = potential_max;
    }
}

template <typename T, typename U> inline void SetToMin(T& current_min, const U potential_min) {
    if (potential_min < current_min) {
        current_min = potential_min;
    }
}

// Mini classes
class Complement {
  public:
    static seqan::FunctorComplement<seqan::Dna5> Dna5;
    static seqan::FunctorComplement<seqan::Dna> Dna;
};

class DominantBase {
  public:
    static constexpr uintReadLen kLastX = 5;

  private:
    seqan::Dna dom_base_;
    std::array<uintReadLen, 5> seq_content_;

    template <typename S> void FindDominant(const S& seq, uintSeqLen cur_pos) {
        // Get count of most appearing base
        uintSeqLen max_content(0);
        for (uintBaseCall nucleotide = 4; nucleotide--;) { // Leave out N for this
            SetToMax(max_content, seq_content_.at(nucleotide));
        }

        if (0 == max_content) {
            // Only N's in kLastX
            if (length(seq) <= cur_pos || 4 == at(seq, cur_pos)) {
                dom_base_ = 0; // We won't use this position anyways for the statistics so just make it an A
            } else {
                dom_base_ = at(seq, cur_pos); // Use cur base as all before where N's
            }
        } else {
            // Get the base that is closest to the current base which matches max_content (in case there are multiple
            // bases with same number of appearances)
            uintSeqLen pos = cur_pos;
            while (max_content != seq_content_.at(at(seq, --pos))) {
                ;
            }
            dom_base_ = at(seq, pos);
        }
    }

  public:
    DominantBase() : dom_base_(0) { seq_content_.fill(0); }

    DominantBase& operator=(const DominantBase& rhs) {
        dom_base_ = rhs.dom_base_;
        seq_content_ = rhs.seq_content_;
        return *this;
    }

    seqan::Dna Get() { return dom_base_; }

    void Clear() { seq_content_.fill(0); }

    template <typename S> void Set(const S& seq, uintSeqLen cur_pos) {
        // Fill seq_content_
        for (uintSeqLen pos = (kLastX < cur_pos ? cur_pos - kLastX : 0); pos < cur_pos; ++pos) {
            ++seq_content_.at(at(seq, pos));
        }

        FindDominant(seq, cur_pos);
    }

    template <typename U, typename S> void Update(U base, const S& seq, uintSeqLen last_pos) {
        // Add new base
        ++seq_content_.at(base);
        // Remove old base if we reached kLastX bases, so that the sum of seq_content_ is actually kLastX
        if (kLastX <= last_pos) {
            --seq_content_.at(at(seq, last_pos - kLastX));
        }

        FindDominant(seq, last_pos + 1);
    }
};

class DominantBaseWithMemory {
  private:
    DominantBase dom_base_;
    seqan::Dna5String memory_;

  public:
    DominantBaseWithMemory() { reserve(memory_, DominantBase::kLastX + 2); }

    DominantBaseWithMemory& operator=(const DominantBaseWithMemory& rhs) {
        dom_base_ = rhs.dom_base_;
        memory_ = rhs.memory_;
        return *this;
    }

    seqan::Dna Get() { return dom_base_.Get(); }

    void Clear() {
        dom_base_.Clear();
        clear(memory_);
    }

    template <typename S> void Set(const S& seq, uintSeqLen cur_pos) {
        resize(memory_, std::min(static_cast<uintSeqLen>(DominantBase::kLastX), cur_pos) + 1);
        for (auto mem_pos = length(memory_); mem_pos--;) {
            at(memory_, mem_pos) = at(seq, cur_pos + mem_pos + 1 - length(memory_));
        }
        dom_base_.Set(memory_, length(memory_) - 1);
    }

    // Update the current base not the last base like in DominantBase as we need to already store the current base
    template <typename U> void Update(U base) {
        if (length(memory_) > DominantBase::kLastX + 1) {
            erase(memory_, 0);
        }
        memory_ += base;

        if (1 < length(memory_)) {
            // Only update when we have a last base.
            dom_base_.Update(at(memory_, length(memory_) - 2), memory_, length(memory_) - 2);
        } else {
            // When we only added the current base we don't want to add a base, but still set the dominant base to the
            // current base
            dom_base_.Set(memory_, 0);
        }
    }
};

// Functions
template <typename T> void Acquire(std::vector<T>& receiver, std::vector<VectorAtomic<T>>& donator) {
    receiver.resize(donator.size());
    for (auto i = donator.size(); i--;) {
        receiver.at(i) = donator.at(i);
    }
    donator.clear();
    donator.shrink_to_fit();
}

inline void CreateDir(const char* file) {
    auto file_path = boost::filesystem::path(file).parent_path();
    if (!file_path.empty()) {
        boost::filesystem::create_directories(file_path);
    }
}

inline bool FileExists(const std::string& file) {
    struct stat buffer;
    return (stat(file.c_str(), &buffer) == 0);
}

inline void DeleteFile(const char* file) {
    boost::filesystem::remove(boost::filesystem::path(file));
}

inline bool GetReSeqDir(std::string& full_dir, const std::string& dir, const std::string& test_file) {
    full_dir = std::string(PROJECT_SOURCE_DIR) + '/' + dir + '/';
    if (!FileExists(full_dir + test_file)) {
        const char* conda_prefix = std::getenv("CONDA_PREFIX");
        if (conda_prefix) {
            full_dir = std::string(conda_prefix) + "/etc/reseq/" + dir + '/';
        }
        if (!conda_prefix || !FileExists(full_dir + test_file)) {
            const char* reseq_folder = std::getenv("RESEQ_FOLDER");
            if (reseq_folder) {
                full_dir = std::string(reseq_folder) + '/' + dir + '/';
            }
            if (!reseq_folder || !FileExists(full_dir + test_file)) {
                printErr << "Could not find " << dir << " folder. Searching at " << std::endl
                         << (std::string(PROJECT_SOURCE_DIR) + '/' + dir + '/') << std::endl
                         << ("${CONDA_PREFIX}/etc/reseq/" + dir + '/') << std::endl
                         << ("${RESEQ_FOLDER}/" + dir + '/') << std::endl;
                return false;
            }
        }
    }

    return true;
}

template <typename T> inline T Divide(T nom, T den) { // Calculates the division: nom/den with proper rounding to int
    // NOLINTNEXTLINE(clang-analyzer-core.DivideZero) caller guarantees den != 0
    return (nom + den / static_cast<T>(2)) / den;
}
template <typename T>
inline T Divide(const std::atomic<T>& nom, T den) { // Calculates the division: nom/den with proper rounding to int
    return (nom + den / static_cast<T>(2)) / den;
}
template <typename T>
inline T Divide(T nom, const std::atomic<T>& den) { // Calculates the division: nom/den with proper rounding to int
    return (nom + den / static_cast<T>(2)) / den;
}
template <typename T>
inline T Divide(const std::atomic<T>& nom,
                const std::atomic<T>& den) { // Calculates the division: nom/den with proper rounding to int
    return (nom + den / static_cast<T>(2)) / den;
}
template <typename T> inline T DivideAndCeil(T nom, T den) { // Calculates the division: nom/den and always ceiling
    return (nom + den - static_cast<T>(1)) / den;
}

inline bool FileExists(const char* file_name) {
    std::ifstream ifile(file_name);
    return static_cast<bool>(ifile);
}

template <typename T> inline const T& getConst(T& object) {
    return object;
}

inline bool IsN(seqan::Dna5 base);
template <typename T> inline bool HasNTemplate(const T& sequence) {
    for (auto base : sequence) {
        if (IsN(base)) {
            return true;
        }
    }

    return false;
}

// Specify possible options so we do not accidentally call function with seqan::CharString etc., which results in wrong
// output
inline bool HasN(const seqan::Dna5String& sequence) {
    return HasNTemplate(sequence);
}
inline bool HasN(const seqan::Dna5StringReverseComplement& sequence) {
    return HasNTemplate(sequence);
}
inline bool HasN(const ConstDna5StringReverseComplement& sequence) {
    return HasNTemplate(sequence);
}

constexpr uint64_t IntPow(uint16_t base, uint16_t power) {
    uint64_t result = 1;
    for (auto i = power; i--;) {
        result *= base;
    }
    return result;
}

inline double InvLogit2(const double bias) {
    return 2 / (1 + exp(-bias));
}

inline bool IsGC(const seqan::Dna5String& base) {
    return base == 'G' || base == 'C';
}

// Specify possible options instead of using a template so we do not accidentally call function with seqan::CharString
// etc., which results in wrong output
inline bool IsGC(const seqan::Dna5String& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}
inline bool IsGC(const seqan::DnaString& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}
inline bool IsGC(const seqan::Dna5StringReverseComplement& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}
inline bool IsGC(const seqan::DnaStringReverseComplement& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}
inline bool IsGC(const ConstDna5StringReverseComplement& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}
inline bool IsGC(const ConstDnaStringReverseComplement& sequence, uintSeqLen pos) {
    return IsGC(at(sequence, pos));
}

inline bool IsN(seqan::Dna5 base) {
    return 3 < static_cast<uintBaseCall>(base);
}

// Specify possible options instead of using a template so we do not accidentally call function with seqan::CharString
// etc., which results in wrong output
inline bool IsN(const seqan::Dna5String& sequence, uintSeqLen pos) {
    return IsN(at(sequence, pos));
}
inline bool IsN(const seqan::Dna5StringReverseComplement& sequence, uintSeqLen pos) {
    return IsN(at(sequence, pos));
}
inline bool IsN(const ConstDna5StringReverseComplement& sequence, uintSeqLen pos) {
    return IsN(at(sequence, pos));
}

template <typename T>
inline T MeanWithRoundingToFirst(T first,
                                 T second) { // Mean between first and second and rounding in the direction of first
    return (first + second + (first > second ? static_cast<T>(1) : static_cast<T>(0))) / static_cast<T>(2);
}

template <typename T, typename U>
inline uintPercent Percent(T nom, U den) { // Calculates the percentage: nom*100/den with proper rounding to int
    return Divide(static_cast<T>(nom * 100), den);
}
template <typename T, typename U>
inline uintPercent Percent(const std::atomic<T>& nom,
                           U den) { // Calculates the percentage: nom*100/den with proper rounding to int
    return Divide(static_cast<T>(nom * 100), den);
}

inline seqan::Dna5String ReverseComplementorDna5(seqan::Dna5String org) {
    return ConstDna5StringReverseComplement(org);
}
inline seqan::DnaString ReverseComplementorDna(seqan::DnaString org) {
    return ConstDnaStringReverseComplement(org);
}

template <typename T, typename U> inline uintPercent SafePercent(T nom, U den) {
    if (den) {
        return Percent(nom, den);
    } else {
        return 50;
    }
}

template <typename T> void SetDimensions(std::vector<T>& vec, size_t dim1, size_t dim2) {
    vec.resize(dim1);
    for (auto i = dim1; i--;) {
        vec.at(i).resize(dim2);
    }
}

template <typename T, typename... Args>
void SetDimensions(std::vector<T>& vec, size_t dim1, size_t dim2, Args... args) {
    vec.resize(dim1);
    for (auto i = dim1; i--;) {
        SetDimensions(vec.at(i), dim2, args...);
    }
}

template <typename T> inline int Sign(T val) {
    return (T(0) < val) - (val < T(0));
}

inline uintSeqLen TransformDistanceToStartOfErrorRegion(uintSeqLen start_dist_error_region) {
    return (start_dist_error_region + 9) / 10;
}

inline unsigned int TrueRandom() {
    std::random_device rd;
    return rd();
}
} // namespace utilities
#endif // SWIG
} // namespace reseq

#endif // UTILITIES_H
