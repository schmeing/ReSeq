#include "cli/convert_profile.h"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "DataStats.h"
#include "logging.hpp"
#include "ProbabilityEstimates.h"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::cerr;
using std::exception;
using std::string;

namespace reseq::cli {

namespace {

/// Save to a temporary file then rename over the target for safe in-place conversion.
bool SafeRename(const string& tmp_path, const string& target_path) {
    std::error_code ec;
    std::filesystem::rename(tmp_path, target_path, ec);
    if (ec) {
        printErr << "Failed to rename '" << tmp_path << "' to '" << target_path << "': " << ec.message() << std::endl;
        std::filesystem::remove(tmp_path);
        return false;
    }
    return true;
}

} // anonymous namespace

int RunConvertProfile(const std::vector<std::string>& args, const variables_map& general_opts,
                      options_description& opt_desc_full) {
    options_description opt_desc("convertProfile");
    opt_desc.add_options()("statsIn,s", value<string>(), "Input stats file (.reseq)")(
        "statsOut,o", value<string>(), "Output stats file [overwrites input if omitted]")(
        "probsIn,p", value<string>(), "Input probabilities file (.reseq.ipf)")(
        "probsOut,P", value<string>(), "Output probabilities file [overwrites input if omitted]")(
        "textFormat",
        "Write legacy text format (default: compressed binary). "
        "Text format is portable across platforms; binary format is faster but not portable across different "
        "architectures or compilers");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq convertProfile -s <stats.reseq> [-p <probs.reseq.ipf>] [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse convertProfile command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
        return 0;
    }

    bool text_format = opts_map.count("textFormat");
    bool has_stats = opts_map.count("statsIn");
    bool has_probs = opts_map.count("probsIn");

    if (!has_stats && !has_probs) {
        printErr << "At least one of statsIn (-s) or probsIn (-p) is required." << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (has_stats) {
        string stats_in = opts_map["statsIn"].as<string>();
        string stats_out;
        bool in_place = false;
        if (opts_map.count("statsOut")) {
            stats_out = opts_map["statsOut"].as<string>();
        } else {
            stats_out = stats_in + ".tmp";
            in_place = true;
        }

        printInfo << "Loading stats from " << stats_in << std::endl;
        DataStats data_stats(nullptr);
        if (!data_stats.Load(stats_in.c_str())) {
            return 1;
        }

        printInfo << "Saving stats to " << (in_place ? stats_in : stats_out)
                  << (text_format ? " (text format)" : " (compressed binary)") << std::endl;
        if (!data_stats.Save(stats_out.c_str(), text_format)) {
            return 1;
        }

        if (in_place && !SafeRename(stats_out, stats_in)) {
            return 1;
        }
    }

    if (has_probs) {
        string probs_in = opts_map["probsIn"].as<string>();
        string probs_out;
        bool in_place = false;
        if (opts_map.count("probsOut")) {
            probs_out = opts_map["probsOut"].as<string>();
        } else {
            probs_out = probs_in + ".tmp";
            in_place = true;
        }

        printInfo << "Loading probabilities from " << probs_in << std::endl;
        ProbabilityEstimates probs;
        if (!probs.Load(probs_in.c_str())) {
            return 1;
        }

        printInfo << "Saving probabilities to " << (in_place ? probs_in : probs_out)
                  << (text_format ? " (text format)" : " (compressed binary)") << std::endl;
        if (!probs.Save(probs_out.c_str(), text_format)) {
            return 1;
        }

        if (in_place && !SafeRename(probs_out, probs_in)) {
            return 1;
        }
    }

    printInfo << "Conversion complete." << std::endl;
    return 0;
}

} // namespace reseq::cli
