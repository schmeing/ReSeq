#include "cli/query_profile.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "DataStats.h"
#include "FragmentDistributionStats.h"
#include "logging.hpp"
#include "Reference.h"
#include "utilities.hpp"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using reseq::DataStats;
using reseq::Reference;
using reseq::uintNumThreads;
using std::cerr;
using std::cout;
using std::exception;
using std::max;
using std::string;

namespace reseq::cli {

int RunQueryProfile(const std::vector<std::string>& args, uintNumThreads num_threads, const variables_map& general_opts,
                    options_description& opt_desc_full) {
    options_description opt_desc("queryProfile");
    opt_desc.add_options()("fragLenBias", value<string>(),
                           "Output fragment length bias to file (tsv format; - for stdout)")(
        "maxLenDeletion", "Output lengths of longest detected deletion to stdout")(
        "maxReadLength", "Output lengths of longest detected read to stdout")(
        "ref,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")(
        "refSeqBias", value<string>(), "Output reference sequence bias to file (tsv format; - for stdout)")(
        "stats,s", value<string>(), "Reseq statistics file to extract reference sequence bias");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq queryProfile -r <ref.fa> -s <stats.reseq>\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse queryProfile command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
    } else {
        auto it_ref = opts_map.find("ref");
        auto it_refseq_bias = opts_map.find("refSeqBias");
        string ref_file = "";
        if (opts_map.end() == it_ref) {
            if (opts_map.end() != it_refseq_bias) {
                printErr << "ref option is required for refSeqBias output." << std::endl;
                if (0 < kVerbosityLevel) {
                    cerr << usage_str;
                    cerr << opt_desc_full << std::endl;
                }
                return 1;
            }
        } else {
            ref_file = it_ref->second.as<string>();
            printInfo << "Reading reference from " << ref_file << std::endl;
        }

        auto it_stats = opts_map.find("stats");
        if (opts_map.end() == it_stats) {
            printErr << "stats option is mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc_full << std::endl;
            }
            return 1;
        } else {
            auto stats_file = it_stats->second.as<string>();
            printInfo << "Reading reference sequence biases from " << stats_file << std::endl;

            string fraglen_bias_file = "";
            auto it_fraglen_bias = opts_map.find("fragLenBias");
            if (opts_map.end() != it_fraglen_bias) {
                fraglen_bias_file = it_fraglen_bias->second.as<string>();
                if ("-" == fraglen_bias_file) {
                    fraglen_bias_file = "";
                    printInfo << "Writing fragment length biases to stdout" << std::endl;
                } else {
                    printInfo << "Writing fragment length biases to " << fraglen_bias_file << std::endl;
                }
            }

            string refseq_bias_file = "";
            if (opts_map.end() != it_refseq_bias) {
                refseq_bias_file = it_refseq_bias->second.as<string>();
                if ("-" == refseq_bias_file) {
                    refseq_bias_file = "";
                    printInfo << "Writing reference sequence biases to stdout" << std::endl;
                } else {
                    printInfo << "Writing reference sequence biases to " << refseq_bias_file << std::endl;
                }
            }

            DataStats real_data_stats(nullptr);
            Reference species_reference;
            if (opts_map.end() != it_ref) {
                if (species_reference.ReadFasta(ref_file.c_str())) {
                    real_data_stats.SetReference(&species_reference);
                } else {
                    return 1;
                }
            }

            if (real_data_stats.Load(stats_file.c_str())) {
                bool error = false;
                bool no_output = true;

                if (opts_map.end() != it_fraglen_bias) {
                    if ("" == fraglen_bias_file) {
                        cout << "fragLenBias:" << std::endl;
                    }
                    if (!real_data_stats.FragmentDistribution().WriteFragLenBias(fraglen_bias_file)) {
                        error = true;
                    }
                    no_output = false;
                }

                if (opts_map.count("maxLenDeletion")) {
                    real_data_stats.CalculateMaxLenDeletion();
                    cout << "maxLenDeletion: " << real_data_stats.Errors().MaxLenDeletion() << std::endl;
                    no_output = false;
                }

                if (opts_map.count("maxReadLength")) {
                    cout << "maxReadLength: "
                         << max(real_data_stats.ReadLengths(0).to(), real_data_stats.ReadLengths(1).to()) - 1
                         << std::endl;
                    no_output = false;
                }

                if (opts_map.end() != it_refseq_bias) {
                    if ("" == refseq_bias_file) {
                        cout << "refSeqBias:" << std::endl;
                    }
                    if (!real_data_stats.FragmentDistribution().WriteRefSeqBias(refseq_bias_file, species_reference)) {
                        error = true;
                    }
                    no_output = false;
                }

                if (no_output) {
                    printErr << "No output option was selected." << std::endl;
                    if (0 < kVerbosityLevel) {
                        cerr << usage_str;
                        cerr << opt_desc_full << std::endl;
                    }
                    return 1;
                }

                if (error) {
                    return 1;
                }
            } else {
                return 1;
            }
        }
    }

    return 0;
}

} // namespace reseq::cli
