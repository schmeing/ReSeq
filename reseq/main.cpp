#include <atomic>
#include <exception>
using std::exception;
#include <iostream>
using std::cerr;
#include <cstdint>
#include <string>
using std::string;
using std::to_string;
#include <vector>
using std::vector;

namespace reseq {
std::atomic<uint16_t> kVerbosityLevel{99};
bool kNoDebugOutput = false;
} // namespace reseq
#include "logging.hpp"
using reseq::kVerbosityLevel;

#include <boost/program_options.hpp>
using boost::program_options::command_line_parser;
using boost::program_options::include_positional;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::parsed_options;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;

#include "cli/convert_profile.h"
#include "cli/illumina_pe.h"
#include "cli/query_profile.h"
#include "cli/replace_n.h"
#include "cli/seq_to_illumina.h"
#include "CMakeConfig.h"
#include "utilities.hpp"
using reseq::uintNumThreads;

// Definitions so that referencing a const static is valid
seqan::FunctorComplement<seqan::Dna5> reseq::utilities::Complement::Dna5;
seqan::FunctorComplement<seqan::Dna> reseq::utilities::Complement::Dna;

// Main
int main(int argc, char* argv[]) {
    uintNumThreads num_threads;
    uint16_t verbosity_opt = 4;
    options_description opt_desc_full("General");
    opt_desc_full.add_options() // Returns a special object with defined operator ()
        ("help,h", "Prints help information and exits")(
            "threads,j", value<uintNumThreads>(&num_threads)->default_value(0), "Number of threads used (0=auto)")(
            "verbosity", value<uint16_t>(&verbosity_opt)->default_value(4),
            "Sets the level of verbosity (4=everything, 0=nothing)")("version", "Prints version info and exits");

    vector<string> unrecognized_opts;
    variables_map general_opts_map;
    try {
        parsed_options general_opts = command_line_parser(argc, argv).options(opt_desc_full).allow_unregistered().run();
        unrecognized_opts = collect_unrecognized(general_opts.options, include_positional);

        store(general_opts, general_opts_map);
        notify(general_opts_map);
        kVerbosityLevel.store(verbosity_opt);
    } catch (const exception& e) {
        printErr << "Could not parse general command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts_map.count("version")) { // Check if user only wants to know version
        cerr << "ReSeq version " << RESEQ_VERSION_MAJOR << '.' << RESEQ_VERSION_MINOR << std::endl;
        return 0;
    }

    string general_usage =
        string("\nProgram: reseq (REal SEQuence replicator)\n") + "Version: " + to_string(RESEQ_VERSION_MAJOR) + '.' +
        to_string(RESEQ_VERSION_MINOR) + '\n' + "Contact: Stephan Schmeing <stephan.schmeing@uzh.ch>\n\n" +
        "Usage:  reseq <command> [options]\n" + "Commands:\n" + "  illuminaPE\t\t" +
        "simulates illumina paired-end data\n" + "  queryProfile\t\t" +
        "queries reseq statistic files for information\n" + "  replaceN\t\t" + "replaces N's in reference\n" +
        "  seqToIllumina\t\t" + "applies illumina quality and error model to input sequences\n" + "  convertProfile\t" +
        "converts profiles between text and binary formats\n";

    int return_code = 0;
    if (0 == unrecognized_opts.size()) {
        cerr << general_usage << std::endl;
    } else {
        printInfo << "Running ReSeq version " << RESEQ_VERSION_MAJOR << '.'
                  << RESEQ_VERSION_MINOR; // Always show version

        if ("queryProfile" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in queryProfile mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunQueryProfile(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
        } else if ("replaceN" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in replaceN mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunReplaceN(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
        } else if ("illuminaPE" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in illuminaPE mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunIlluminaPE(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
        } else if ("seqToIllumina" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in seqToIllumina mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunSeqToIllumina(unrecognized_opts, num_threads, general_opts_map, opt_desc_full);
        } else if ("convertProfile" == unrecognized_opts.at(0)) {
            if (2 < kVerbosityLevel) {
                cerr << " in convertProfile mode" << std::endl;
            }
            unrecognized_opts.erase(unrecognized_opts.begin());
            return_code = reseq::cli::RunConvertProfile(unrecognized_opts, general_opts_map, opt_desc_full);
        } else {
            if (2 < kVerbosityLevel) {
                cerr << std::endl;
            }
            printErr << "Unrecognized command: '" << unrecognized_opts.at(0) << "'" << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << general_usage << std::endl;
            }
            return 1;
        }
    }

    return return_code;
}
