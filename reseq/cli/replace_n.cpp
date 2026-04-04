#include "cli/replace_n.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "cli/cli_common.h"
#include "logging.hpp"
#include "Reference.h"
#include "utilities.hpp"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using reseq::Reference;
using reseq::uintNumThreads;
using reseq::uintSeed;
using std::cerr;
using std::exception;
using std::string;

namespace reseq::cli {

int RunReplaceN(const std::vector<std::string>& args, uintNumThreads num_threads, const variables_map& general_opts,
                options_description& opt_desc_full) {
    options_description opt_desc("ReplaceN");
    opt_desc.add_options()("refIn,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")(
        "refSim,R", value<string>(),
        "File to where reference sequences in fasta format with N's randomly replace should be written to")(
        "seed", value<uintSeed>(), "Seed used for replacing N, if none is given random seed will be used");
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq replaceN -r <refIn.fa> -R <refSim.fa> [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse replaceN command line arguments: " << e.what() << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    }

    if (general_opts.count("help")) {
        cerr << usage_str;
        cerr << opt_desc_full << std::endl;
    } else if (!AutoDetectThreads(num_threads, opt_desc_full, usage_str)) {
        return 1;
    } else {
        auto it_ref_in = opts_map.find("refIn");
        auto it_ref_out = opts_map.find("refSim");
        string ref_input, ref_output;
        if (opts_map.end() == it_ref_in) {
            printErr << "refIn option is mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc_full << std::endl;
            }
            return 1;
        } else {
            ref_input = it_ref_in->second.as<string>();
            printInfo << "Reading reference from " << ref_input << std::endl;

            if (opts_map.end() == it_ref_out) {
                printErr << "refSim option is mandatory." << std::endl;
                if (0 < kVerbosityLevel) {
                    cerr << usage_str;
                    cerr << opt_desc_full << std::endl;
                }
                return 1;
            } else {
                ref_output = it_ref_out->second.as<string>();
                printInfo << "Writing reference without N to " << ref_output << std::endl;

                Reference species_reference;
                if (species_reference.ReadFasta(ref_input.c_str())) {
                    auto seed = GetSeed(opts_map);
                    species_reference.ReplaceN(seed);

                    if (species_reference.WriteFasta(ref_output.c_str())) {
                        printInfo << "Finished replacing N's." << std::endl;
                    } else {
                        return 1;
                    }
                } else {
                    return 1;
                }
            }
        }
    }

    return 0;
}

} // namespace reseq::cli
