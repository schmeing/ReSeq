#include "cli/cli_common.h"

#include <cstdio>
#include <string>
#include <thread>

#include <boost/program_options.hpp>

#include "logging.hpp"
#include "utilities.hpp"

using reseq::uintNumThreads;
using reseq::uintSeed;
using reseq::utilities::TrueRandom;
using std::string;

namespace reseq::cli {

bool AutoDetectThreads(uintNumThreads& num_threads, const boost::program_options::options_description& opt_desc,
                       const string& usage_str) {
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

void GetProbsOut(string& probs_out, const string& fallback_out, const boost::program_options::variables_map& opts_map) {
    auto it_probs_out = opts_map.find("probabilitiesOut");
    if (opts_map.end() == it_probs_out) {
        probs_out = fallback_out;
    } else {
        probs_out = it_probs_out->second.as<string>();
    }
}

void PrepareProbabilityEstimation(string& probs_in, string& probs_out, const string& standard_probs_out,
                                  bool loaded_stats, const boost::program_options::variables_map& opts_map) {
    auto it_probs_in = opts_map.find("probabilitiesIn");
    if (opts_map.end() == it_probs_in) {
        GetProbsOut(probs_out, standard_probs_out, opts_map);

        probs_in = "";
        if (loaded_stats) {
            // Check if standard probabilities output file exists and in case it does use it as input for the
            // probabilities
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
