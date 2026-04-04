#include "cli/illumina_pe.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "cli/cli_common.h"
#include "DataStats.h"
#include "FragmentDistributionStats.h"
#include "logging.hpp"
#include "ProbabilityEstimates.h"
#include "Reference.h"
#include "Simulator.h"
#include "utilities.hpp"

using boost::program_options::command_line_parser;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variables_map;
using reseq::DataStats;
using reseq::kVerbosityLevel;
using reseq::ProbabilityEstimates;
using reseq::Reference;
using reseq::RefSeqBiasSimulation;
using reseq::Simulator;
using reseq::uintFragCount;
using reseq::uintNumFits;
using reseq::uintNumThreads;
using reseq::uintQualPrint;
using reseq::uintSeed;
using reseq::uintSeqLen;
using reseq::utilities::DeleteFile;
using reseq::utilities::FileExists;
using reseq::utilities::GetReSeqDir;
using std::cerr;
using std::count;
using std::exception;
using std::string;

namespace {

bool DefaultExtensionFile(string& file_name, const string& extension) {
    if (FileExists(file_name)) {
        // File exists without modification
        return true;
    }

    if (0 == count(file_name.begin(), file_name.end(), '.')) {
        file_name += extension;
        if (FileExists(file_name)) {
            // File exists with adding extension
            return true;
        }
    }

    if (0 == count(file_name.begin(), file_name.end(), '/')) {
        string adapter_dir;
        if (GetReSeqDir(adapter_dir, "adapters", "TruSeq_single.fa")) {
            file_name = adapter_dir + file_name;
        } else {
            file_name = "";
            return false;
        }

        if (FileExists(file_name)) {
            // File exists with adding folder
            return true;
        }
    }

    printErr << "Automatic filename deduction failed for: '" << file_name << "'. Please provide valid path."
             << std::endl;
    file_name = "";
    return false;
}

void GetDataStats(DataStats& real_data_stats, string& stats_file, bool& loaded_stats, bool stats_only,
                  uintSeqLen max_ref_seq_bin_size, const variables_map& opts_map, const options_description& opt_desc,
                  const string& usage_str, uintNumThreads num_threads, bool text_format, bool both_formats) {
    auto it_bam_in = opts_map.find("bamIn");
    auto it_stats_in = opts_map.find("statsIn");
    auto it_stats_out = opts_map.find("statsOut");
    bool no_tiles = opts_map.count("noTiles");
    bool no_bias_calculation = opts_map.count("noBias");
    bool tiles = opts_map.count("tiles");
    if (no_tiles && tiles) {
        printErr << "noTiles and tiles options are exclusive." << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc << std::endl;
        }
    }

    if (opts_map.end() == it_bam_in) { // "bamIn" hasn't been found
        if (stats_only) {
            printErr << "With statsOnly option bamIn option is mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc << std::endl;
            }
        } else if (no_tiles || tiles) {
            printErr << (no_tiles ? "noTiles" : "tiles")
                     << " option can only be specified for generation of statistics. After that it depends on the "
                        "statistics loaded."
                     << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc << std::endl;
            }
        } else if (opts_map.end() == it_stats_in) { // "statsIn" hasn't been found
            printErr << "Either bamIn or statsIn option are mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc << std::endl;
            }
        } else {
            if (opts_map.end() != it_stats_out) {
                printWarn << "statsOut option ignored as statsIn is also given." << std::endl;
            }

            stats_file = it_stats_in->second.as<string>();
            printInfo << "Reading real data statistics from " << stats_file << std::endl;

            if (real_data_stats.Load(stats_file.c_str())) {
                if (no_bias_calculation) {
                    real_data_stats.SetUniformBias();
                }
                real_data_stats.PrepareProcessing();

                loaded_stats = true;
            }
        }
    } else {
        if (opts_map.end() != it_stats_in) {
            printErr << "statsIn and bamIn options are exclusive." << std::endl;
        } else if (!real_data_stats.HasReference()) {
            printErr << "Reading in statistics requires refIn option." << std::endl;
        } else {
            string bam_input(it_bam_in->second.as<string>());
            printInfo << "Reading mapping from " << bam_input << std::endl;

            if (opts_map.end() == it_stats_out) {
                stats_file = bam_input + ".reseq";
            } else {
                stats_file = it_stats_out->second.as<string>();
            }
            printInfo << "Storing real data statistics in " << stats_file << std::endl;

            // Set adapter files to input or default
            bool adapters_failed = false;
            string adapter_file, adapter_matrix;
            auto it_adapter_file = opts_map.find("adapterFile");
            auto it_adapter_matrix = opts_map.find("adapterMatrix");

            if (opts_map.end() != it_adapter_file) {
                adapter_file = it_adapter_file->second.as<string>();
                if (!DefaultExtensionFile(adapter_file, ".fa")) {
                    adapters_failed = true;
                }
            }
            if (opts_map.end() != it_adapter_matrix) {
                adapter_matrix = it_adapter_matrix->second.as<string>();
                if (!DefaultExtensionFile(adapter_matrix, ".mat")) {
                    adapters_failed = true;
                }
            }

            if (!adapters_failed) {
                if (adapter_file.size()) {
                    if (0 == adapter_matrix.size()) {
                        // Only sequences
                        if (2 < adapter_file.size()) {
                            adapter_matrix = adapter_file;
                            adapter_matrix.replace(adapter_matrix.size() - 3, 3, ".mat");
                            if (!FileExists(adapter_matrix)) {
                                adapters_failed = true;
                            }
                        } else {
                            adapters_failed = true;
                        }

                        if (adapters_failed) {
                            printErr << "Switching extension of adapter-sequence file did not result in valid matrix "
                                        "file. Please specify the matrix file."
                                     << std::endl;
                            if (0 < kVerbosityLevel) {
                                cerr << usage_str;
                                cerr << opt_desc << std::endl;
                            }
                        }
                    }
                } else if (adapter_matrix.size()) {
                    // Only matrix
                    if (3 < adapter_matrix.size()) {
                        adapter_file = adapter_matrix;
                        adapter_file.replace(adapter_file.size() - 4, 4, ".fa");
                        if (!FileExists(adapter_file)) {
                            adapters_failed = true;
                        }
                    } else {
                        adapters_failed = true;
                    }

                    if (adapters_failed) {
                        printErr << "Switching extension of adapter-matrix file did not result in valid sequence file. "
                                    "Please specify the sequence file."
                                 << std::endl;
                        if (0 < kVerbosityLevel) {
                            cerr << usage_str;
                            cerr << opt_desc << std::endl;
                        }
                    }
                }
            }

            if (!adapters_failed) {
                if (!tiles) {
                    printInfo << "Tiles will be ignored and all statistics will be generated like all reads are from "
                                 "the same tile."
                              << std::endl;
                    real_data_stats.IgnoreTiles();
                } else {
                    printInfo << "Statistics will be split into tiles if possible." << std::endl;
                }

                string variant_file = "";
                auto it_variant_file = opts_map.find("vcfIn");
                if (opts_map.end() != it_variant_file) {
                    variant_file = it_variant_file->second.as<string>();
                    printInfo << "Ignoring all variant positions listed in '" << variant_file
                              << "' for error statistics." << std::endl;
                }

                if (real_data_stats.ReadBam(bam_input.c_str(), adapter_file.c_str(), adapter_matrix.c_str(),
                                            variant_file, max_ref_seq_bin_size, num_threads, !no_bias_calculation)) {
                    real_data_stats.Save(stats_file.c_str(), text_format);
                    if (both_formats) {
                        string alt_path = stats_file + (text_format ? ".bin" : ".text");
                        real_data_stats.Save(alt_path.c_str(), !text_format);
                    }
                    real_data_stats.PrepareProcessing();
                    real_data_stats.ClearReference(); // It shouldn't be used after the read in to guarantee that we can
                                                      // remove or change the reference
                } else {
                    DeleteFile(
                        stats_file.c_str()); // Remove the stats file if one existed previously so it is clear that we
                                             // encountered an error and do not accidentally continue with the old file
                }
                loaded_stats = false;
            }
        }
    }
}

bool WriteSysError(string& sys_error_file, uintSeed& seed, bool stop_after_estimation, const variables_map& opts_map,
                   const options_description& opt_desc, const string& usage_str, const Reference& ref,
                   const DataStats& stats, const ProbabilityEstimates& estimates) {
    auto it_read = opts_map.find("readSysError");
    auto it_write = opts_map.find("writeSysError");

    if (opts_map.end() == it_write) {
        if (stop_after_estimation) {
            if (opts_map.end() != it_read) {
                printWarn << "readSysError option is only for simulation, so it's useless with the stopAfterEstimation "
                             "option."
                          << std::endl;
                if (0 < kVerbosityLevel) {
                    cerr << usage_str;
                    cerr << opt_desc << std::endl;
                }
                return false;
            }
        } else {
            if (opts_map.end() != it_read) {
                sys_error_file = it_read->second.as<string>();
            } else {
                sys_error_file = "";
            }

            seed = reseq::cli::GetSeed(opts_map);
        }
    } else {
        if (opts_map.end() != it_read) {
            printErr << "writeSysError and readSysError option are mutually exclusive. Specify the one or the other."
                     << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc << std::endl;
            }
            return false;
        } else {
            seed = reseq::cli::GetSeed(opts_map);

            sys_error_file = it_write->second.as<string>();
            Simulator sim;
            if (!sim.CreateSystematicErrorProfile(sys_error_file.c_str(), ref, stats, estimates, seed)) {
                DeleteFile(sys_error_file
                               .c_str()); // Remove the sys error file if one existed previously so it is clear that we
                                          // encountered an error and do not accidentally continue with the old file
                return false;
            }
        }
    }

    return true;
}

void PrepareSimulation(string& sim_output_first, string& sim_output_second, const variables_map& opts_map) {
    auto it = opts_map.find("firstReadsOut");
    if (opts_map.end() == it) {
        sim_output_first = "reseq-R1.fq";
    } else {
        sim_output_first = it->second.as<string>();
    }

    it = opts_map.find("secondReadsOut");
    if (opts_map.end() == it) {
        sim_output_second = "reseq-R2.fq";
    } else {
        sim_output_second = it->second.as<string>();
    }
    printInfo << "Storing simulated data in " << sim_output_first << " and " << sim_output_second << std::endl;
}

} // anonymous namespace

namespace reseq::cli {

int RunIlluminaPE(const std::vector<std::string>& args, uintNumThreads num_threads, const variables_map& general_opts,
                  options_description& opt_desc_full) {
    uintNumFits ipf_iterations;
    double ipf_precision;
    uintFragCount num_read_pairs;
    double coverage;
    double error_multiplier;
    uintSeqLen maximum_insert_length;
    uintQualPrint minimum_mapping_quality;
    uintSeqLen max_ref_seq_bin_size;
    std::string record_base_identifier, meth_file;

    options_description opt_desc("Stats");
    opt_desc.add_options()("adapterFile", value<string>(), "Fasta file with adapter sequences [(AutoDetect)]")(
        "adapterMatrix", value<string>(),
        "0/1 matrix with valid adapter pairing (first read in rows, second read in columns) "
        "[(AutoDetect)]")("bamIn,b", value<string>(), "Position sorted bam/sam file with reads mapped to refIn")(
        "binSizeBiasFit", value<uintSeqLen>(&max_ref_seq_bin_size)->default_value(100000000),
        "Reference sequences large then this are split for bias fitting to limit memory consumption")(
        "maxFragLen", value<uintSeqLen>(&maximum_insert_length)->default_value(2000),
        "Maximum fragment length to include pairs into statistics")(
        "minMapQ", value<uintQualPrint>(&minimum_mapping_quality)->default_value(2),
        "Minimum mapping quality to include pairs into statistics")(
        "noBias", "Do not perform bias fit. Results in uniform coverage if simulated from")(
        "noTiles", "Ignore tiles for the statistics [default]")(
        "refIn,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")(
        "statsOnly", "Only generate the statistics")("statsIn,s", value<string>(),
                                                     "Skips statistics generation and reads directly from stats file")(
        "statsOut,S", value<string>(), "Stores the real data statistics for reuse in given file [<bamIn>.reseq]")(
        "tiles", "Use tiles for the statistics")("vcfIn,v", value<string>(),
                                                 "Ignore all positions with a listed variant for stats generation")(
        "textFormat",
        "Write profile files in legacy text format instead of compressed binary. "
        "Text format is portable across platforms; binary format is faster but not portable across different "
        "architectures or compilers")("bothFormats", "Write profile files in both binary and text formats. The "
                                                     "alternate format is saved with a .text or .bin suffix");

    options_description opt_desc_ipf("Probabilities");
    opt_desc_ipf.add_options()("ipfIterations", value<uintNumFits>(&ipf_iterations)->default_value(200),
                               "Maximum number of iterations for iterative proportional fitting")(
        "ipfPrecision", value<double>(&ipf_precision)->default_value(5),
        "Iterative proportional fitting procedure stops after reaching this precision (%)")(
        "probabilitiesIn,p", value<string>(),
        "Loads last estimated probabilities and continues from there if precision is not met [<statsIn>.ipf]")(
        "probabilitiesOut,P", value<string>(),
        "Stores the probabilities estimated by iterative proportional fitting [<probabilitiesIn>]")(
        "stopAfterEstimation", "Stop after estimating the probabilities");
    opt_desc.add(opt_desc_ipf);

    options_description opt_desc_sim("Simulation");
    opt_desc_sim.add_options()("firstReadsOut,1", value<string>(),
                               "Writes the simulated first reads into this file [reseq-R1.fq]")(
        "secondReadsOut,2", value<string>(), "Writes the simulated second reads into this file [reseq-R2.fq]")(
        "coverage,c", value<double>(&coverage)->default_value(0.0),
        "Approximate average read depth simulated (0 = Corrected original coverage)")(
        "errorMutliplier", value<double>(&error_multiplier)->default_value(1.0),
        "Divides the original probability of correct base calls(no substitution error) by this value and renormalizes")(
        "methylation", value<string>(&meth_file)->default_value(""),
        "Extended bed graph file specifying methylation for regions. Multiple score columns for individual "
        "alleles are possible, but must match with vcfSim. C->T conversions for 1-specified value in region.")(
        "noInDelErrors", "Simulate reads without InDel errors")("noSubstitutionErrors",
                                                                "Simulate reads without substitution errors")(
        "numReads", value<uintFragCount>(&num_read_pairs)->default_value(0),
        "Approximate number of read pairs simulated (0 = Use <coverage>)")(
        "readSysError", value<string>(),
        "Read systematic errors from file in fastq format (seq=dominant error, qual=error percentage)")(
        "recordBaseIdentifier", value<string>(&record_base_identifier)->default_value("ReseqRead"),
        "Base Identifier for the simulated fastq records, followed by a count and other information about the read")(
        "refBias", value<string>(),
        "Way to select the reference biases for simulation (keep [from refIn]/no [biases]/draw [with replacement from "
        "original biases]/file) [keep/no]")(
        "refBiasFile", value<string>(), "File to read reference biases from: One sequence per file (identifier bias)")(
        "refSim,R", value<string>(), "Reference sequences in fasta format to simulate from [<refIn>]")(
        "seed", value<uintSeed>(), "Seed used for simulation, if none is given random seed will be used")(
        "vcfSim,V", value<string>(), "Defines genotypes to simulate alleles or populations")(
        "writeSysError", value<string>(),
        "Write the randomly drawn systematic errors to file in fastq format (seq=dominant error, qual=error "
        "percentage)");
    opt_desc.add(opt_desc_sim);
    opt_desc_full.add(opt_desc);

    string usage_str = "Usage:  reseq illuminaPE -b <file.bam> -r <ref.fa> -1 <file1.fq> -2 <file2.fq> [options]\n";
    variables_map opts_map;
    try {
        store(command_line_parser(args).options(opt_desc).run(), opts_map);
        notify(opts_map);
    } catch (const exception& e) {
        printErr << "Could not parse illuminaPE command line arguments: " << e.what() << std::endl;
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
    } else if (0 != num_read_pairs && 0.0 != coverage) {
        printErr << "numReads and coverage set the same value. Use either the one or the other." << std::endl;
        if (0 < kVerbosityLevel) {
            cerr << usage_str;
            cerr << opt_desc_full << std::endl;
        }
        return 1;
    } else if (!AutoDetectThreads(num_threads, opt_desc_full, usage_str)) {
        return 1;
    } else {
        auto it_ref_in = opts_map.find("refIn");
        auto it_ref_out = opts_map.find("refSim");
        bool stop_after_estimation = opts_map.count("stopAfterEstimation");
        bool text_format = opts_map.count("textFormat");
        bool both_formats = opts_map.count("bothFormats");
        if (opts_map.end() == it_ref_in && opts_map.end() == it_ref_out &&
            (!stop_after_estimation || !opts_map.count("statsIn") || opts_map.count("writeSysError"))) {
            printErr << "refIn or refSim option mandatory." << std::endl;
            if (0 < kVerbosityLevel) {
                cerr << usage_str;
                cerr << opt_desc_full << std::endl;
            }
            return 1;
        } else {
            string ref_input, ref_output;
            if (opts_map.end() != it_ref_in) {
                ref_input = it_ref_in->second.as<string>();
                printInfo << "Reading reference from " << ref_input << std::endl;
            }
            if (opts_map.end() != it_ref_out) {
                ref_output = it_ref_out->second.as<string>();
                printInfo << "Simulating reads from reference " << ref_output << std::endl;
            }

            Reference species_reference;
            if (opts_map.end() != it_ref_in && !species_reference.ReadFasta(ref_input.c_str())) {
                return 1;
            } else {
                DataStats real_data_stats((opts_map.end() == it_ref_in ? nullptr : &species_reference),
                                          maximum_insert_length, minimum_mapping_quality);
                bool stats_only = opts_map.count("statsOnly");
                string stats_file;
                bool loaded_stats = false; // defensive init (see design spec)

                GetDataStats(real_data_stats, stats_file, loaded_stats, stats_only, max_ref_seq_bin_size, opts_map,
                             opt_desc_full, usage_str, num_threads, text_format, both_formats);

                if (0 == real_data_stats.TotalNumberReads()) {
                    return 1;
                }

                if (!stats_only) { // Data stats have been loaded or computed
                    string probs_in, probs_out;

                    PrepareProbabilityEstimation(probs_in, probs_out, stats_file + ".ipf", loaded_stats, opts_map);

                    ProbabilityEstimates probabilities;

                    if (!probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads,
                                                probs_out.c_str(), probs_in.c_str(), text_format)) {
                        return 1;
                    }
                    if (both_formats) {
                        string alt_probs = probs_out + (text_format ? ".bin" : ".text");
                        probabilities.Save(alt_probs.c_str(), !text_format);
                    }

                    if (!stop_after_estimation || opts_map.count("writeSysError")) {
                        if ((opts_map.end() != it_ref_out) &&
                            (opts_map.end() == it_ref_in || ref_input != ref_output) &&
                            !species_reference.ReadFasta(ref_output.c_str())) {
                            return 1;
                        } else {
                            probabilities.PrepareResult();

                            string sys_error_file;
                            uintSeed seed;
                            if (!WriteSysError(sys_error_file, seed, stop_after_estimation, opts_map, opt_desc_full,
                                               usage_str, species_reference, real_data_stats, probabilities)) {
                                return 1;
                            } else {
                                if (!stop_after_estimation) {
                                    string sim_output_first, sim_output_second;

                                    PrepareSimulation(sim_output_first, sim_output_second, opts_map);

                                    RefSeqBiasSimulation ref_bias_model;
                                    auto it_ref_bias = opts_map.find("refBias");
                                    auto it_ref_bias_file = opts_map.find("refBiasFile");
                                    if (opts_map.end() == it_ref_bias) {
                                        if (opts_map.end() == it_ref_bias_file) {
                                            if (opts_map.end() == it_ref_out) {
                                                printInfo << "Keeping reference sequence biases from read in[default]."
                                                          << std::endl;
                                                ref_bias_model = RefSeqBiasSimulation::kKeep;
                                            } else {
                                                printInfo << "Removing all reference sequence biases[default]."
                                                          << std::endl;
                                                ref_bias_model = RefSeqBiasSimulation::kNo;
                                            }
                                        } else {
                                            printInfo << "Reading reference sequence biases from file." << std::endl;
                                            ref_bias_model = RefSeqBiasSimulation::kFile;
                                        }
                                    } else {
                                        auto ref_bmodel = it_ref_bias->second.as<string>();
                                        if ("keep" == ref_bmodel) {
                                            printInfo << "Keeping reference sequence biases from read in." << std::endl;
                                            ref_bias_model = RefSeqBiasSimulation::kKeep;
                                        } else if ("no" == ref_bmodel) {
                                            printInfo << "Removing all reference sequence biases." << std::endl;
                                            ref_bias_model = RefSeqBiasSimulation::kNo;
                                        } else if ("draw" == ref_bmodel) {
                                            printInfo << "Drawing with replacement from reference sequence biases from "
                                                         "read in."
                                                      << std::endl;
                                            ref_bias_model = RefSeqBiasSimulation::kDraw;
                                        } else if ("file" == ref_bmodel) {
                                            printInfo << "Reading reference sequence biases from file." << std::endl;
                                            ref_bias_model = RefSeqBiasSimulation::kFile;
                                        } else {
                                            printErr << "Unknown option for refBias: " << ref_bmodel << std::endl;
                                            if (0 < kVerbosityLevel) {
                                                cerr << usage_str;
                                                cerr << opt_desc_full << std::endl;
                                            }
                                            ref_bias_model = RefSeqBiasSimulation::kError;
                                        }
                                    }

                                    if (RefSeqBiasSimulation::kError == ref_bias_model) {
                                        return 1;
                                    }

                                    string ref_bias_file;
                                    if (RefSeqBiasSimulation::kFile == ref_bias_model) {
                                        if (opts_map.end() == it_ref_bias_file) {
                                            printErr << "refBiasFile option mandatory if for refBias option file was "
                                                        "chosen"
                                                     << std::endl;
                                            if (0 < kVerbosityLevel) {
                                                cerr << usage_str;
                                                cerr << opt_desc_full << std::endl;
                                            }
                                            ref_bias_model = RefSeqBiasSimulation::kError;
                                        } else {
                                            ref_bias_file = it_ref_bias_file->second.as<string>();
                                        }
                                    } else {
                                        if (opts_map.end() != it_ref_bias_file) {
                                            printErr
                                                << "refBiasFile option only allowed if for refBias option file was "
                                                   "chosen"
                                                << std::endl;
                                            if (0 < kVerbosityLevel) {
                                                cerr << usage_str;
                                                cerr << opt_desc_full << std::endl;
                                            }
                                            ref_bias_model = RefSeqBiasSimulation::kError;
                                        }
                                    }

                                    if (RefSeqBiasSimulation::kError == ref_bias_model) {
                                        return 1;
                                    } else {
                                        // Load variation
                                        auto it_var_file = opts_map.find("vcfSim");

                                        if (opts_map.end() == it_var_file && opts_map.count("vcfIn") &&
                                            opts_map.count("statsIn")) {
                                            printErr << "vcfIn specified but not used as stats were loaded. Did you "
                                                        "mean vcfSim?"
                                                     << std::endl;
                                            if (0 < kVerbosityLevel) {
                                                cerr << usage_str;
                                                cerr << opt_desc_full << std::endl;
                                            }
                                            return 1;
                                        } else {
                                            string var_file("");
                                            if (opts_map.end() == it_var_file) {
                                                printInfo << "Simulating reference allele/population." << std::endl;
                                            } else {
                                                var_file = it_var_file->second.as<string>();
                                                printInfo << "Simulating variance from file: '" << var_file << "'"
                                                          << std::endl;
                                            }

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
                                                    printErr << "--error_multiplier must be a positive value."
                                                             << std::endl;
                                                    return 1;
                                                } else {
                                                    probabilities.ChangeErrorRate(error_multiplier);
                                                }
                                            }

                                            Simulator sim;
                                            if (!sim.Simulate(sim_output_first.c_str(), sim_output_second.c_str(),
                                                              species_reference, real_data_stats, probabilities,
                                                              num_threads, seed, num_read_pairs, coverage,
                                                              ref_bias_model, ref_bias_file, sys_error_file,
                                                              record_base_identifier, var_file, meth_file)) {
                                                return 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

} // namespace reseq::cli
