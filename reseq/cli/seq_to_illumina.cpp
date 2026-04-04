#include "cli/seq_to_illumina.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "cli/cli_common.h"
#include "DataStats.h"
#include "logging.hpp"
#include "ProbabilityEstimates.h"
#include "Simulator.h"
#include "utilities.hpp"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using reseq::DataStats;
using reseq::ProbabilityEstimates;
using reseq::Simulator;
using reseq::uintNumFits;
using reseq::uintNumThreads;
using reseq::uintSeed;
using std::cerr;
using std::exception;
using std::string;

namespace reseq::cli {

int RunSeqToIllumina(const std::vector<std::string>& args, uintNumThreads num_threads,
                     const variables_map& general_opts, options_description& opt_desc_full) {
    double error_multiplier;
    uintNumFits ipf_iterations;
    double ipf_precision;

    options_description opt_desc("seqToIllumina");
    opt_desc.add_options()(
        "errorMutliplier", value<double>(&error_multiplier)->default_value(1.0),
        "Divides the original probability of correct base calls(no substitution error) by this value "
        "and renormalizes")("input,i", value<string>(), "Input file (fasta format, gz and bz2 supported) [stdin]")(
        "ipfIterations", value<uintNumFits>(&ipf_iterations)->default_value(200),
        "Maximum number of iterations for iterative proportional fitting")(
        "ipfPrecision", value<double>(&ipf_precision)->default_value(5),
        "Iterative proportional fitting procedure stops after reaching this precision (%)")(
        "noInDelErrors", "Simulate reads without InDel errors")("noSubstitutionErrors",
                                                                "Simulate reads without substitution errors")(
        "output,o", value<string>(), "Output file (fastq format, gz and bz2 supported) [stdout]")(
        "probabilitiesIn,p", value<string>(),
        "Loads last estimated probabilities and continues from there if precision is not met [<statsIn>.ipf]")(
        "probabilitiesOut,P", value<string>(),
        "Stores the probabilities estimated by iterative proportional fitting [<probabilitiesIn>]")(
        "seed", value<uintSeed>(), "Seed used for simulation, if none is given random seed will be used")(
        "statsIn,s", value<string>(), "Profile file that contains the statistics used for simulation")(
        "textFormat",
        "Write profile files in legacy text format instead of compressed binary. "
        "Text format is portable across platforms; binary format is faster but not portable across different "
        "architectures or compilers")("bothFormats", "Write profile files in both binary and text formats. The "
                                                     "alternate format is saved with a .text or .bin suffix");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq seqToIllumina -i <input.fa> -o <output.fq> -s <stats.reseq> [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse seqToIllumina command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
    } else if (ipf_precision < 0.0) {
        printErr << "ipfPrecision must be positive." << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    } else if (!AutoDetectThreads(num_threads, opt_desc_full, usage_str)) {
        return 1;
    } else {
        bool text_format = opts_map.count("textFormat");
        bool both_formats = opts_map.count("bothFormats");
        DataStats real_data_stats(nullptr);
        string probs_in, probs_out;

        auto it_stats_in = opts_map.find("statsIn");
        if (opts_map.end() == it_stats_in) { // "statsIn" hasn't been found
            printErr << "statsIn option is mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc << std::endl;
            }
        } else {
            auto stats_file = it_stats_in->second.as<string>();
            printInfo << "Reading real data statistics from " << stats_file << std::endl;

            if (real_data_stats.Load(stats_file.c_str())) {
                real_data_stats.PrepareProcessing();

                PrepareProbabilityEstimation(probs_in, probs_out, stats_file + ".ipf", true, opts_map);
            }
        }

        if (0 == real_data_stats.TotalNumberReads()) {
            return 1;
        } else {
            ProbabilityEstimates probabilities;
            if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads, probs_out.c_str(),
                                        probs_in.c_str(), text_format)) {
                return 1;
            } else {
                if (both_formats && probs_out.size()) {
                    string alt_probs = probs_out + (text_format ? ".bin" : ".text");
                    probabilities.Save(alt_probs.c_str(), !text_format);
                }
                probabilities.PrepareResult();

                string org_seq_file;
                auto it = opts_map.find("input");
                if (opts_map.end() == it) {
                    printInfo << "Reading original sequences from stdin" << std::endl;
                } else {
                    org_seq_file = it->second.as<string>();
                    printInfo << "Reading original sequences from " << org_seq_file << std::endl;
                }

                string destination_file;
                it = opts_map.find("output");
                if (opts_map.end() == it) {
                    printInfo << "Writing simulated data to stdout" << std::endl;
                } else {
                    destination_file = it->second.as<string>();
                    printInfo << "Writing simulated data to " << destination_file << std::endl;
                }

                auto seed = GetSeed(opts_map);

                if (opts_map.count("noInDelErrors")) {
                    probabilities.RemoveInDelErrors();
                }
                if (opts_map.count("noSubstitutionErrors")) {
                    probabilities.RemoveSubstitutionErrors();
                }

                if (1.0 != error_multiplier) {
                    if (0.0 == error_multiplier) {
                        probabilities.RemoveSubstitutionErrors();
                    } else if (0.0 > error_multiplier) {
                        printErr << "--error_multiplier must be a positive value." << std::endl;
                        return 1;
                    } else {
                        probabilities.ChangeErrorRate(error_multiplier);
                    }
                }

                Simulator sim;
                if (!sim.SimulateErrorModelOnly(destination_file, org_seq_file, real_data_stats, probabilities,
                                                num_threads, seed)) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

} // namespace reseq::cli
