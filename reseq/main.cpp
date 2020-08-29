#include <algorithm>
using std::count;
using std::max;
#include <exception>
using std::exception;
#include <iostream>
using std::cout;
using std::cerr;
#include <stdint.h>
#include <stdio.h>
#include <string>
using std::string;
using std::to_string;
#include <sys/stat.h>
#include <thread>
#include <vector>
using std::vector;

#include "gtest/gtest.h"

namespace reseq{
	uint16_t kVerbosityLevel = 99;
	bool kNoDebugOutput = false;
}
#include "reportingUtils.hpp"
using reseq::kVerbosityLevel;

#include <boost/program_options.hpp>
using boost::program_options::command_line_parser;
using boost::program_options::include_positional;
using boost::program_options::notify;
using boost::program_options::options_description;
using boost::program_options::parsed_options;
using boost::program_options::store;
using boost::program_options::value;
using boost::program_options::variable_value;
using boost::program_options::variables_map;

#include "CMakeConfig.h"
#include "DataStats.h"
using reseq::DataStats;
#include "FragmentDistributionStats.h"
using reseq::RefSeqBiasSimulation;
#include "ProbabilityEstimates.h"
using reseq::ProbabilityEstimates;
#include "Reference.h"
using reseq::Reference;
#include "Simulator.h"
using reseq::Simulator;
#include "utilities.hpp"
using reseq::uintQualPrint;
using reseq::uintNumThreads;
using reseq::uintNumFits;
using reseq::uintSeqLen;
using reseq::uintFragCount;
using reseq::uintSeed;
using reseq::utilities::DeleteFile;
using reseq::utilities::TrueRandom;

#include "AdapterStatsTest.h"
#include "DataStatsTest.h"
#include "FragmentDistributionStatsTest.h"
#include "FragmentDuplicationStatsTest.h"
#include "ProbabilityEstimatesTest.h"
#include "ReferenceTest.h"
#include "SeqQualityStatsTest.h"
#include "SimulatorTest.h"
#include "SurroundingTest.h"
#include "TileStatsTest.h"
#include "VectTest.h"
#include "utilitiesTest.h"

// Helper functions
bool AutoDetectThreads(uintNumThreads &num_threads, const options_description &opt_desc, const string &usage_str ){
	if(0 == num_threads){
		num_threads = std::thread::hardware_concurrency();
		if(0 == num_threads){
			printErr << "Automatic detection of available cores failed." << std::endl;
			if(0 < kVerbosityLevel){
				cerr << usage_str;
				cerr << opt_desc << std::endl;
			}
			return false;
		}
		else{
			printInfo << "Detected " << num_threads << " cores to be used." << std::endl;
		}
	}

	return true;
}

void DefaultExtensionFile(string &file_name, const string folder, const string extension){
	if (FILE *file = fopen(file_name.c_str(), "r")){
		// File exists without modification
		fclose(file);
		return;
	}

	if( 0 == count(file_name.begin(), file_name.end(), '.' ) ){
		file_name += extension;
		if (FILE *file = fopen(file_name.c_str(), "r")){
			// File exists with adding extension
			fclose(file);
			return;
		}
	}

	if( 0 == count(file_name.begin(), file_name.end(), '/' ) ){
		file_name = folder + file_name;
		if (FILE *file = fopen(file_name.c_str(), "r")){
			// File exists with adding folder
			fclose(file);
			return;
		}
	}

	printErr << "Automatic filename deduction failed for: '" << file_name << "'. Please provide valid path." << std::endl;
	file_name = "";
	return;
}

void GetDataStats(DataStats &real_data_stats, string &stats_file, bool &loaded_stats, bool stats_only, uintSeqLen max_ref_seq_bin_size, const variables_map &opts_map, const options_description &opt_desc, const string &usage_str, uintNumThreads num_threads ){
	auto it_bam_in = opts_map.find("bamIn");
	auto it_stats_in = opts_map.find("statsIn");
	auto it_stats_out = opts_map.find("statsOut");
	bool no_tiles = opts_map.count("noTiles");
	bool no_bias_calculation = opts_map.count("noBias");
	bool tiles = opts_map.count("tiles");
	if(no_tiles && tiles){
		printErr << "noTiles and tiles options are exclusive." << std::endl;
		if(0 < kVerbosityLevel){
			cerr << usage_str;
			cerr << opt_desc << std::endl;
		}
	}

	if(opts_map.end() == it_bam_in){ // "bamIn" hasn't been found
		if( stats_only ){
			printErr << "With statsOnly option bamIn option is mandatory." << std::endl;
			if(0 < kVerbosityLevel){
				cerr << usage_str;
				cerr << opt_desc << std::endl;
			}
		}
		else if( no_tiles || tiles ){
			printErr << (no_tiles ? "noTiles" : "tiles") << " option can only be specified for generation of statistics. After that it depends on the statistics loaded." << std::endl;
			if(0 < kVerbosityLevel){
				cerr << usage_str;
				cerr << opt_desc << std::endl;
			}
		}
		else if(opts_map.end() == it_stats_in){ // "statsIn" hasn't been found
			printErr << "Either bamIn or statsIn option are mandatory." << std::endl;
			if(0 < kVerbosityLevel){
				cerr << usage_str;
				cerr << opt_desc << std::endl;
			}
		}
		else{
			if(opts_map.end() != it_stats_out){
				printWarn << "statsOut option ignored as statsIn is also given." << std::endl;
			}

			stats_file = it_stats_in->second.as<string>();
			printInfo << "Reading real data statistics from " << stats_file << std::endl;

			if( real_data_stats.Load( stats_file.c_str() ) ){
				if(no_bias_calculation){
					real_data_stats.SetUniformBias();
				}
				real_data_stats.PrepareProcessing();

				loaded_stats = true;
			}
		}
	}
	else{
		if(opts_map.end() != it_stats_in){
			printErr << "statsIn and bamIn options are exclusive." << std::endl;
		}
		else if(!real_data_stats.HasReference()){
			printErr << "Reading in statistics requires refIn option." << std::endl;
		}
		else{
			string bam_input( it_bam_in->second.as<string>() );
			printInfo << "Reading mapping from " << bam_input << std::endl;

			if(opts_map.end() == it_stats_out){
				stats_file = bam_input + ".reseq";
			}
			else{
				stats_file = it_stats_out->second.as<string>();
			}
			printInfo << "Storing real data statistics in " << stats_file << std::endl;

			// Set adapter files to input or default
			string adapter_file, adapter_matrix;
			auto it_adapter_file = opts_map.find("adapterFile");
			auto it_adapter_matrix = opts_map.find("adapterMatrix");
			if( opts_map.end() != it_adapter_file || opts_map.end() != it_adapter_matrix ){
				if( opts_map.end() != it_adapter_file ){
					adapter_file = it_adapter_file->second.as<string>();
					DefaultExtensionFile(adapter_file, string(PROJECT_SOURCE_DIR)+"/adapters/", ".fa");
				}
				if( opts_map.end() != it_adapter_matrix ){
					adapter_matrix = it_adapter_matrix->second.as<string>();
					DefaultExtensionFile(adapter_matrix, string(PROJECT_SOURCE_DIR)+"/adapters/", ".mat");
				}

				if( adapter_file.size() ){
					if( 0 == adapter_matrix.size() ){
						// Only sequences
						bool failed = true;
						if(2 < adapter_file.size()){
							adapter_matrix = adapter_file;
							adapter_matrix.replace(adapter_matrix.size()-3,3,".mat");
							if (FILE *file = fopen(adapter_matrix.c_str(), "r")){
								fclose(file);
								failed = false;
							}
						}
						if(failed){
							printErr << "Switching extension of adapter-sequence file did not result in valid matrix file. Please specify the matrix file." << std::endl;
							if(0 < kVerbosityLevel){
								cerr << usage_str;
								cerr << opt_desc << std::endl;
							}

							adapter_file = "";
							adapter_matrix = "error";
						}
					}
				}
				else{
					// Only matrix
					bool failed = true;
					if(3 < adapter_matrix.size()){
						adapter_file = adapter_matrix;
						adapter_file.replace(adapter_file.size()-4,4,".fa");
						if (FILE *file = fopen(adapter_file.c_str(), "r")){
							fclose(file);
							failed = false;
						}
					}
					if(failed){
						printErr << "Switching extension of adapter-matrix file did not result in valid sequence file. Please specify the sequence file." << std::endl;
						if(0 < kVerbosityLevel){
							cerr << usage_str;
							cerr << opt_desc << std::endl;
						}
						adapter_file = "";
						adapter_matrix = "error";
					}
				}
			}

			if( !tiles ){
				printInfo << "Tiles will be ignored and all statistics will be generated like all reads are from the same tile." << std::endl;
				real_data_stats.IgnoreTiles();
			}
			else{
				printInfo << "Statistics will be split into tiles if possible." << std::endl;
			}

			if( !adapter_file.empty() || adapter_matrix.empty() ){
				string variant_file = "";
				auto it_variant_file = opts_map.find("vcfIn");
				if( opts_map.end() != it_variant_file ){
					variant_file = it_variant_file->second.as<string>();
					printInfo << "Ignoring all variant positions listed in '" << variant_file << "' for error statistics." << std::endl;
				}

				if( real_data_stats.ReadBam( bam_input.c_str(), adapter_file.c_str(), adapter_matrix.c_str(), variant_file, max_ref_seq_bin_size, num_threads, !no_bias_calculation ) ){
					real_data_stats.Save( stats_file.c_str() );
					real_data_stats.PrepareProcessing();
					real_data_stats.ClearReference(); // It shouldn't be used after the read in to guarantee that we can remove or change the reference
				}
				else{
					DeleteFile( stats_file.c_str() ); // Remove the stats file if one existed previously so it is clear that we encountered an error and do not accidentally continue with the old file
				}
				loaded_stats = false;
			}
		}
	}
}

void GetProbsOut( string &probs_out, const string &fallback_out, const variables_map &opts_map){
	auto it_probs_out = opts_map.find("probabilitiesOut");
	if(opts_map.end() == it_probs_out){
		probs_out = fallback_out;
	}
	else{
		probs_out = it_probs_out->second.as<string>();
	}
}

void PrepareProbabilityEstimation( string &probs_in, string &probs_out, const string &standard_probs_out, bool loaded_stats, const variables_map &opts_map ){
	auto it_probs_in = opts_map.find("probabilitiesIn");
	if(opts_map.end() == it_probs_in){
		GetProbsOut( probs_out, standard_probs_out, opts_map);

		probs_in = "";
		if( loaded_stats ){
			// Check if standard probabilities output file exists and in case it does use it as input for the probabilities
			if (FILE *file = fopen(standard_probs_out.c_str(), "r")){
				probs_in = standard_probs_out;
				fclose(file);
			}
		}
	}
	else{
		probs_in = it_probs_in->second.as<string>();
		GetProbsOut( probs_out, probs_in, opts_map);
	}
}

uintSeed GetSeed(const variables_map &opts_map){
	uintSeed seed;
	auto it = opts_map.find("seed");
	if(opts_map.end() == it){
		seed = TrueRandom();
		printInfo << "Randomly generated seed is " << seed << std::endl;
	}
	else{
		seed = it->second.as<uintSeed>();
		printInfo << "Using seed " << seed << std::endl;
	}

	return seed;
}

bool WriteSysError(string &sys_error_file, uintSeed &seed, bool stop_after_estimation, const variables_map &opts_map, const options_description &opt_desc, const string &usage_str, const Reference &ref, const DataStats &stats, const ProbabilityEstimates &estimates){
	auto it_read = opts_map.find("readSysError");
	auto it_write = opts_map.find("writeSysError");

	if(opts_map.end() == it_write){
		if(stop_after_estimation){
			if(opts_map.end() != it_read ){
				printWarn << "readSysError option is only for simulation, so it's useless with the stopAfterEstimation option." << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc << std::endl;
				}
				return false;
			}
		}
		else{
			if(opts_map.end() != it_read ){
				sys_error_file = it_read->second.as<string>();
			}
			else{
				sys_error_file = "";
			}

			seed = GetSeed(opts_map);
		}
	}
	else{
		if(opts_map.end() != it_read){
			printErr << "writeSysError and readSysError option are mutually exclusive. Specify the one or the other." << std::endl;
			if(0 < kVerbosityLevel){
				cerr << usage_str;
				cerr << opt_desc << std::endl;
			}
			return false;
		}
		else{
			seed = GetSeed(opts_map);

			sys_error_file = it_write->second.as<string>();
			Simulator sim;
			if( !sim.CreateSystematicErrorProfile(sys_error_file.c_str(), ref, stats, estimates, seed) ){
				DeleteFile(sys_error_file.c_str()); // Remove the sys error file if one existed previously so it is clear that we encountered an error and do not accidentally continue with the old file
				return false;
			}
		}
	}

	return true;
}

void PrepareSimulation( string &sim_output_first, string &sim_output_second, const variables_map &opts_map ){
	auto it = opts_map.find("firstReadsOut");
	if(opts_map.end() == it){
		sim_output_first = "reseq-R1.fq";
	}
	else{
		sim_output_first = it->second.as<string>();
	}

	it = opts_map.find("secondReadsOut");
	if(opts_map.end() == it){
		sim_output_second = "reseq-R2.fq";
	}
	else{
		sim_output_second = it->second.as<string>();
	}
	printInfo << "Storing simulated data in " << sim_output_first << " and " << sim_output_second << std::endl;
}

int RunGoogleTests(const string &new_test, string &tests_already_run){
	if( !tests_already_run.empty() ){
		tests_already_run += ":";
	}
	tests_already_run += new_test;

	::testing::GTEST_FLAG(filter) = new_test;
	return RUN_ALL_TESTS();
}

// Main
int main(int argc, char *argv[]) {
	uintNumThreads num_threads;
	options_description opt_desc_full("General");
	opt_desc_full.add_options() // Returns a special object with defined operator ()
		("help,h", "Prints help information and exits")
		("threads,j", value<uintNumThreads>(&num_threads)->default_value(0), "Number of threads used (0=auto)")
		("verbosity", value<uint16_t>(&reseq::kVerbosityLevel)->default_value(4), "Sets the level of verbosity (4=everything, 0=nothing)")
		("version", "Prints version info and exits");

	vector<string> unrecognized_opts;
	variables_map general_opts_map;
	try{
		parsed_options general_opts = command_line_parser(argc, argv).options(opt_desc_full).allow_unregistered().run();
		unrecognized_opts = collect_unrecognized(general_opts.options, include_positional);

		store(general_opts, general_opts_map);
		notify(general_opts_map);
	}
	catch(const exception& e) {
		printErr << "Could not parse general command line arguments: " << e.what() << std::endl;
		if(0 < kVerbosityLevel){
			cerr << opt_desc_full << std::endl;
		}
		return 1;
	}

	if(general_opts_map.count("version")) { // Check if user only wants to know version
		cerr << "ReSeq version " << RESEQ_VERSION_MAJOR << '.' << RESEQ_VERSION_MINOR << std::endl;
		return 0;
	}

	string general_usage =
		string("\nProgram: reseq (REal SEQuence replicator)\n")+
		"Version: "+to_string(RESEQ_VERSION_MAJOR)+'.'+to_string(RESEQ_VERSION_MINOR)+'\n'+
		"Contact: Stephan Schmeing <stephan.schmeing@uzh.ch>\n\n"+
		"Usage:  reseq <command> [options]\n"+
		"Commands:\n"+
		"  illuminaPE\t\t"+"simulates illumina paired-end data\n"+
		"  queryProfile\t\t"+"queries reseq statistic files for information\n"+
		"  replaceN\t\t"+"replaces N's in reference\n"+
		"  seqToIllumina\t\t"+"applies illumina quality and error model to input sequences\n"+
		"  test\t\t\t"+"tests the program\n";

	int return_code = 0;
	if(0 == unrecognized_opts.size()){
		cerr << general_usage << std::endl;
	}
	else{
		printInfo << "Running ReSeq version " << RESEQ_VERSION_MAJOR << '.' << RESEQ_VERSION_MINOR; // Always show version

		if ("queryProfile" == unrecognized_opts.at(0)) {
			if(2 < kVerbosityLevel){
				cerr << " in queryProfile mode" << std::endl;
			}
			unrecognized_opts.erase(unrecognized_opts.begin());

			options_description opt_desc("queryProfile");
			opt_desc.add_options() // Returns a special object with defined operator ()
				("maxLenDeletion", "Output lengths of longest detected deletion to stdout")
				("maxReadLength", "Output lengths of longest detected read to stdout")
				("ref,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")
				("refSeqBias", value<string>(), "Output reference sequence bias to file (tsv format; - for stdout)")
				("stats,s", value<string>(), "Reseq statistics file to extract reference sequence bias");
			opt_desc_full.add(opt_desc);

			string usage_str = "Usage:  reseq queryProfile -r <ref.fa> -s <stats.reseq>\n";
			variables_map opts_map;
			try{
				store( command_line_parser(unrecognized_opts).options(opt_desc).run(), opts_map );
				notify(opts_map);
			}
			catch (const exception& e) {
				printErr << "Could not parse queryProfile command line arguments: " << e.what() << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}

			if ( general_opts_map.count("help") ) {
				cerr << usage_str;
				cerr << opt_desc_full << std::endl;
			}
			else{
				auto it_ref = opts_map.find("ref");
				auto it_refseq_bias = opts_map.find("refSeqBias");
				string ref_file = "";
				if(opts_map.end() == it_ref){
					if (opts_map.end() != it_refseq_bias){
						printErr << "ref option is required for refSeqBias output." << std::endl;
						if(0 < kVerbosityLevel){
							cerr << usage_str;
							cerr << opt_desc_full << std::endl;
						}
						return 1;
					}
				}
				else{
					ref_file = it_ref->second.as<string>();
					printInfo << "Reading reference from " << ref_file << std::endl;
				}

				auto it_stats = opts_map.find("stats");
				if(opts_map.end() == it_stats){
					printErr << "stats option is mandatory." << std::endl;
					if(0 < kVerbosityLevel){
						cerr << usage_str;
						cerr << opt_desc_full << std::endl;
					}
					return 1;
				}
				else{
					auto stats_file = it_stats->second.as<string>();
					printInfo << "Reading reference sequence biases from " << stats_file << std::endl;

					string refseq_bias_file = "";
					if(opts_map.end() != it_refseq_bias){
						refseq_bias_file = it_refseq_bias->second.as<string>();
						if("-" == refseq_bias_file){
							refseq_bias_file = "";
							printInfo << "Writing reference sequence biases to stdout" << std::endl;
						}
						else{
							printInfo << "Writing reference sequence biases to " << refseq_bias_file << std::endl;
						}
					}

					DataStats real_data_stats(NULL);
					Reference species_reference;
					if(opts_map.end() != it_ref){
						if( species_reference.ReadFasta( ref_file.c_str() ) ){
							real_data_stats.SetReference(&species_reference);
						}
						else{
							return 1;
						}
					}

					if( real_data_stats.Load( stats_file.c_str() ) ){
						bool error = false;
						bool no_output = true;

						if( opts_map.count("maxLenDeletion") ){
							real_data_stats.CalculateMaxLenDeletion();
							cout << "maxLenDeletion: " << real_data_stats.Errors().MaxLenDeletion() << std::endl;
							no_output = false;
						}

						if( opts_map.count("maxReadLength") ){
							cout << "maxReadLength: " << max(real_data_stats.ReadLengths(0).to(), real_data_stats.ReadLengths(1).to())-1 << std::endl;
							no_output = false;
						}

						if( opts_map.end() != it_refseq_bias ){
							if( !real_data_stats.FragmentDistribution().WriteRefSeqBias(refseq_bias_file, species_reference) ){
								error = true;
							}
							no_output = false;
						}

						if(no_output){
							printErr << "No output option was selected." << std::endl;
							if(0 < kVerbosityLevel){
								cerr << usage_str;
								cerr << opt_desc_full << std::endl;
							}
							return 1;
						}

						if(error){
							return 1;
						}
					}
					else{
						return 1;
					}
				}
			}
		}
		else if ("replaceN" == unrecognized_opts.at(0)) {
			if(2 < kVerbosityLevel){
				cerr << " in replaceN mode" << std::endl;
			}
			unrecognized_opts.erase(unrecognized_opts.begin());

			options_description opt_desc("ReplaceN");
			opt_desc.add_options() // Returns a special object with defined operator ()
				("refIn,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")
				("refSim,R", value<string>(), "File to where reference sequences in fasta format with N's randomly replace should be written to")
				("seed", value<uintSeed>(), "Seed used for replacing N, if none is given random seed will be used");
			opt_desc_full.add(opt_desc);

			string usage_str = "Usage:  reseq replaceN -r <refIn.fa> -R <refSim.fa> [options]\n";
			variables_map opts_map;
			try{
				store( command_line_parser(unrecognized_opts).options(opt_desc).run(), opts_map );
				notify(opts_map);
			}
			catch (const exception& e) {
				printErr << "Could not parse replaceN command line arguments: " << e.what() << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}

			if ( general_opts_map.count("help") ) {
				cerr << usage_str;
				cerr << opt_desc_full << std::endl;
			}
			else if(!AutoDetectThreads(num_threads, opt_desc_full, usage_str)){
				return 1;
			}
			else{
				auto it_ref_in = opts_map.find("refIn");
				auto it_ref_out = opts_map.find("refSim");
				string ref_input, ref_output;
				if(opts_map.end() == it_ref_in){
					printErr << "refIn option is mandatory." << std::endl;
					if(0 < kVerbosityLevel){
						cerr << usage_str;
						cerr << opt_desc_full << std::endl;
					}
					return 1;
				}
				else{
					ref_input =  it_ref_in->second.as<string>();
					printInfo << "Reading reference from " << ref_input << std::endl;

					if(opts_map.end() == it_ref_out){
						printErr << "refSim option is mandatory." << std::endl;
						if(0 < kVerbosityLevel){
							cerr << usage_str;
							cerr << opt_desc_full << std::endl;
						}
						return 1;
					}
					else{
						ref_output = it_ref_out->second.as<string>();
						printInfo << "Writing reference without N to " << ref_output << std::endl;

						Reference species_reference;
						if( species_reference.ReadFasta( ref_input.c_str() ) ){
							auto seed = GetSeed(opts_map);

							species_reference.ReplaceN( seed );

							if( species_reference.WriteFasta( ref_output.c_str() ) ){
								printInfo << "Finished replacing N's." << std::endl;
							}
							else{
								return 1;
							}
						}
						else{
							return 1;
						}
					}
				}
			}
		}
		else if("illuminaPE" == unrecognized_opts.at(0)) {
			if(2 < kVerbosityLevel){
				cerr << " in illuminaPE mode" << std::endl;
			}
			unrecognized_opts.erase(unrecognized_opts.begin());

			uintNumFits ipf_iterations;
			double ipf_precision;
			uintFragCount num_read_pairs;
			double coverage;
			double error_mutliplier;
			uintSeqLen maximum_insert_length;
			uintQualPrint minimum_mapping_quality;
			uintSeqLen max_ref_seq_bin_size;
			std::string record_base_identifier, meth_file;

			options_description opt_desc("Stats");
			opt_desc.add_options() // Returns a special object with defined operator ()
				("adapterFile", value<string>(), "Fasta file with adapter sequences [(AutoDetect)]")
				("adapterMatrix", value<string>(), "0/1 matrix with valid adapter pairing (first read in rows, second read in columns) [(AutoDetect)]")
				("bamIn,b", value<string>(), "Position sorted bam/sam file with reads mapped to refIn")
				("binSizeBiasFit", value<uintSeqLen>(&max_ref_seq_bin_size)->default_value(100000000), "Reference sequences large then this are split for bias fitting to limit memory consumption")
				("maxFragLen", value<uintSeqLen>(&maximum_insert_length)->default_value(2000), "Maximum fragment length to include pairs into statistics")
				("minMapQ", value<uintQualPrint>(&minimum_mapping_quality)->default_value(10), "Minimum mapping quality to include pairs into statistics")
				("noBias", "Do not perform bias fit. Results in uniform coverage if simulated from")
				("noTiles", "Ignore tiles for the statistics [default]")
				("refIn,r", value<string>(), "Reference sequences in fasta format (gz and bz2 supported)")
				("statsOnly", "Only generate the statistics")
				("statsIn,s", value<string>(), "Skips statistics generation and reads directly from stats file")
				("statsOut,S", value<string>(), "Stores the real data statistics for reuse in given file [<bamIn>.reseq]")
				("tiles", "Use tiles for the statistics")
				("vcfIn,v", value<string>(), "Ignore all positions with a listed variant for stats generation");

			options_description opt_desc_ipf("Probabilities");
			opt_desc_ipf.add_options()
				("ipfIterations", value<uintNumFits>(&ipf_iterations)->default_value(200), "Maximum number of iterations for iterative proportional fitting")
				("ipfPrecision", value<double>(&ipf_precision)->default_value(5), "Iterative proportional fitting procedure stops after reaching this precision (%)")
				("probabilitiesIn,p", value<string>(), "Loads last estimated probabilities and continues from there if precision is not met [<statsIn>.ipf]")
				("probabilitiesOut,P", value<string>(), "Stores the probabilities estimated by iterative proportional fitting [<probabilitiesIn>]")
				("stopAfterEstimation", "Stop after estimating the probabilities");
			opt_desc.add(opt_desc_ipf);

			options_description opt_desc_sim("Simulation");
			opt_desc_sim.add_options()
				("firstReadsOut,1", value<string>(), "Writes the simulated first reads into this file [reseq-R1.fq]")
				("secondReadsOut,2", value<string>(), "Writes the simulated second reads into this file [reseq-R2.fq]")
				("coverage,c", value<double>(&coverage)->default_value(0.0), "Approximate average read depth simulated (0 = Corrected original coverage)")
				("errorMutliplier", value<double>(&error_mutliplier)->default_value(1.0), "Divides the original probability of correct base calls(no substitution error) by this value and renormalizes")
				("methylation", value<string>(&meth_file)->default_value(""), "Extended bed graph file specifying methylation for regions. Multiple score columns for individual alleles are possible, but must match with vcfSim. C->T conversions for 1-specified value in region.")
				("numReads", value<uintFragCount>(&num_read_pairs)->default_value(0), "Approximate number of read pairs simulated (0 = Use <coverage>)")
				("readSysError", value<string>(), "Read systematic errors from file in fastq format (seq=dominant error, qual=error percentage)")
				("recordBaseIdentifier", value<string>(&record_base_identifier)->default_value("ReseqRead"), "Base Identifier for the simulated fastq records, followed by a count and other information about the read")
				("refBias", value<string>(), "Way to select the reference biases for simulation (keep [from refIn]/no [biases]/draw [with replacement from original biases]/file) [keep/no]")
				("refBiasFile", value<string>(), "File to read reference biases from: One sequence per file (identifier bias)")
				("refSim,R", value<string>(), "Reference sequences in fasta format to simulate from [<refIn>]")
				("seed", value<uintSeed>(), "Seed used for simulation, if none is given random seed will be used")
				("vcfSim,V", value<string>(), "Defines genotypes to simulate alleles or populations")
				("writeSysError", value<string>(), "Write the randomly drawn systematic errors to file in fastq format (seq=dominant error, qual=error percentage)");
			opt_desc.add(opt_desc_sim);
			opt_desc_full.add(opt_desc);

			string usage_str = "Usage:  reseq illuminaPE -b <file.bam> -r <ref.fa> -1 <file1.fq> -2 <file2.fq> [options]\n";
			variables_map opts_map;
			try{
				store( command_line_parser(unrecognized_opts).options(opt_desc).run(), opts_map );
				notify(opts_map);
			}
			catch(const exception& e) {
				printErr << "Could not parse illuminaPE command line arguments: " << e.what() << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}

			if( general_opts_map.count("help") ) {
				cerr << usage_str;
				cerr << opt_desc_full << std::endl;
			}
			else if( ipf_precision < 0.0 ){
				printErr << "ipfPrecision must be positive." << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}
			else if(0 != num_read_pairs && 0.0 != coverage){
				printErr << "numReads and coverage set the same value. Use either the one or the other." << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}
			else if(!AutoDetectThreads(num_threads, opt_desc_full, usage_str)){
				return 1;
			}
			else{
				auto it_ref_in = opts_map.find("refIn");
				auto it_ref_out = opts_map.find("refSim");
				bool stop_after_estimation = opts_map.count("stopAfterEstimation");
				if(opts_map.end() == it_ref_in && opts_map.end() == it_ref_out && (!stop_after_estimation || !opts_map.count("statsIn") || opts_map.count("writeSysError"))){
					printErr << "refIn or refSim option mandatory." << std::endl;
					if(0 < kVerbosityLevel){
						cerr << usage_str;
						cerr << opt_desc_full << std::endl;
					}
					return 1;
				}
				else{
					string ref_input, ref_output;
					if(opts_map.end() != it_ref_in){
						ref_input =  it_ref_in->second.as<string>();
						printInfo << "Reading reference from " << ref_input << std::endl;
					}
					if(opts_map.end() != it_ref_out){
						ref_output = it_ref_out->second.as<string>();
						printInfo << "Simulating reads from reference " << ref_output << std::endl;
					}

					Reference species_reference;
					if( opts_map.end() != it_ref_in && !species_reference.ReadFasta( ref_input.c_str() ) ){
						return 1;
					}
					else{
						DataStats real_data_stats( (opts_map.end() == it_ref_in ? NULL : &species_reference), maximum_insert_length, minimum_mapping_quality);
						bool stats_only = opts_map.count("statsOnly");
						string stats_file;
						bool loaded_stats;

						GetDataStats(real_data_stats, stats_file, loaded_stats, stats_only, max_ref_seq_bin_size, opts_map, opt_desc_full, usage_str, num_threads);

						if( 0 == real_data_stats.TotalNumberReads() ){
							return 1;
						}

						if( !stats_only ){ // Data stats have been loaded or computed
							string probs_in, probs_out;

							PrepareProbabilityEstimation( probs_in, probs_out, stats_file+".ipf", loaded_stats, opts_map );

							ProbabilityEstimates probabilities;

							if( !probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads, probs_out.c_str(), probs_in.c_str()) ){
								return 1;
							}

							if( !stop_after_estimation || opts_map.count("writeSysError") ){
								if((opts_map.end() != it_ref_out) && (opts_map.end() == it_ref_in || ref_input != ref_output) && !species_reference.ReadFasta( ref_output.c_str() )){
									return 1;
								}
								else{
									probabilities.PrepareResult();

									string sys_error_file;
									uintSeed seed;
									if( !WriteSysError(sys_error_file, seed, stop_after_estimation, opts_map, opt_desc_full, usage_str, species_reference, real_data_stats, probabilities) ){
										return 1;
									}
									else{
										if( !stop_after_estimation ){
											string sim_output_first, sim_output_second;

											PrepareSimulation( sim_output_first, sim_output_second, opts_map);

											RefSeqBiasSimulation ref_bias_model;
											auto it_ref_bias = opts_map.find("refBias");
											if(opts_map.end() == it_ref_bias){
												if(opts_map.end() == it_ref_out){
													printInfo << "Keeping reference sequence biases from read in[default]." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kKeep;
												}
												else{
													printInfo << "Removing all reference sequence biases[default]." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kNo;
												}
											}
											else{
												auto ref_bmodel = it_ref_bias->second.as<string>();
												if("keep" == ref_bmodel){
													printInfo << "Keeping reference sequence biases from read in." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kKeep;
												}
												else if("no" == ref_bmodel){
													printInfo << "Removing all reference sequence biases." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kNo;
												}
												else if("draw" == ref_bmodel){
													printInfo << "Drawing with replacement from reference sequence biases from read in." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kDraw;
												}
												else if("file" == ref_bmodel){
													printInfo << "Reading reference sequence biases from file." << std::endl;
													ref_bias_model = RefSeqBiasSimulation::kFile;
												}
												else{
													printErr << "Unknown option for refBias: " << ref_bmodel << std::endl;
													if(0 < kVerbosityLevel){
														cerr << usage_str;
														cerr << opt_desc_full << std::endl;
													}
													ref_bias_model = RefSeqBiasSimulation::kError;
												}
											}

											if(RefSeqBiasSimulation::kError == ref_bias_model){
												return 1;
											}

											string ref_bias_file;
											auto it_ref_bias_file = opts_map.find("refBiasFile");
											if(RefSeqBiasSimulation::kFile == ref_bias_model){
												if(opts_map.end() == it_ref_bias_file){
													printErr << "refBiasFile option mandatory if for refBias option file was chosen" << std::endl;
													if(0 < kVerbosityLevel){
														cerr << usage_str;
														cerr << opt_desc_full << std::endl;
													}
													ref_bias_model = RefSeqBiasSimulation::kError;
												}
												else{
													ref_bias_file = it_ref_bias_file->second.as<string>();
												}
											}
											else{
												if(opts_map.end() != it_ref_bias_file){
													printErr << "refBiasFile option only allowed if for refBias option file was chosen" << std::endl;
													if(0 < kVerbosityLevel){
														cerr << usage_str;
														cerr << opt_desc_full << std::endl;
													}
													ref_bias_model = RefSeqBiasSimulation::kError;
												}
											}

											if(RefSeqBiasSimulation::kError == ref_bias_model){
												return 1;
											}
											else{
												// Load variation
												auto it_var_file = opts_map.find("vcfSim");

												if(opts_map.end() == it_var_file && opts_map.count("vcfIn") && opts_map.count("statsIn")){
													printErr << "vcfIn specified but not used as stats were loaded. Did you mean vcfSim?" << std::endl;
													if(0 < kVerbosityLevel){
														cerr << usage_str;
														cerr << opt_desc_full << std::endl;
													}
													return 1;
												}
												else{
													string var_file("");
													if(opts_map.end() == it_var_file){
														printInfo << "Simulating reference allele/population." << std::endl;
													}
													else{
														var_file = it_var_file->second.as<string>();
														printInfo << "Simulating variance from file: '" << var_file << "'" << std::endl;
													}

													if(1.0 != error_mutliplier){
														probabilities.ChangeErrorRate( error_mutliplier );
													}

													Simulator sim;
													if( !sim.Simulate(sim_output_first.c_str(), sim_output_second.c_str(), species_reference, real_data_stats, probabilities, num_threads, seed, num_read_pairs, coverage, ref_bias_model, ref_bias_file, sys_error_file, record_base_identifier, var_file, meth_file) ){
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
		}
		else if("seqToIllumina" == unrecognized_opts.at(0)){
			if(2 < kVerbosityLevel){
				cerr << " in seqToIllumina mode" << std::endl;
			}
			unrecognized_opts.erase(unrecognized_opts.begin());

			double error_mutliplier;
			uintNumFits ipf_iterations;
			double ipf_precision;

			options_description opt_desc("seqToIllumina");
			opt_desc.add_options() // Returns a special object with defined operator ()
				("errorMutliplier", value<double>(&error_mutliplier)->default_value(1.0), "Divides the original probability of correct base calls(no substitution error) by this value and renormalizes")
				("input,i", value<string>(), "Input file (fasta format, gz and bz2 supported) [stdin]")
				("ipfIterations", value<uintNumFits>(&ipf_iterations)->default_value(200), "Maximum number of iterations for iterative proportional fitting")
				("ipfPrecision", value<double>(&ipf_precision)->default_value(5), "Iterative proportional fitting procedure stops after reaching this precision (%)")
				("output,o", value<string>(), "Output file (fastq format, gz and bz2 supported) [stdout]")
				("probabilitiesIn,p", value<string>(), "Loads last estimated probabilities and continues from there if precision is not met [<statsIn>.ipf]")
				("probabilitiesOut,P", value<string>(), "Stores the probabilities estimated by iterative proportional fitting [<probabilitiesIn>]")
				("seed", value<uintSeed>(), "Seed used for simulation, if none is given random seed will be used")
				("statsIn,s", value<string>(), "Profile file that contains the statistics used for simulation");
			opt_desc_full.add(opt_desc);

			string usage_str = "Usage:  reseq seqToIllumina -i <input.fa> -o <output.fq> -s <stats.reseq> [options]\n";
			variables_map opts_map;
			try{
				store( command_line_parser(unrecognized_opts).options(opt_desc).run(), opts_map );
				notify(opts_map);
			}
			catch(const exception& e) {
				printErr << "Could not parse seqToIllumina command line arguments: " << e.what() << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}

			if( general_opts_map.count("help") ) {
				cerr << usage_str;
				cerr << opt_desc_full << std::endl;
			}
			else if( ipf_precision < 0.0 ){
				printErr << "ipfPrecision must be positive." << std::endl;
				if(0 < kVerbosityLevel){
					cerr << usage_str;
					cerr << opt_desc_full << std::endl;
				}
				return 1;
			}
			else if(!AutoDetectThreads(num_threads, opt_desc_full, usage_str)){
				return 1;
			}
			else{
				DataStats real_data_stats(NULL);
				string probs_in, probs_out;

				auto it_stats_in = opts_map.find("statsIn");
				if(opts_map.end() == it_stats_in){ // "statsIn" hasn't been found
					printErr << "statsIn option is mandatory." << std::endl;
					if(0 < kVerbosityLevel){
						cerr << usage_str;
						cerr << opt_desc << std::endl;
					}
				}
				else{
					auto stats_file = it_stats_in->second.as<string>();
					printInfo << "Reading real data statistics from " << stats_file << std::endl;

					if( real_data_stats.Load( stats_file.c_str() ) ){
						real_data_stats.PrepareProcessing();

						PrepareProbabilityEstimation( probs_in, probs_out, stats_file+".ipf", true, opts_map );
					}
				}

				if( 0 == real_data_stats.TotalNumberReads() ){
					return 1;
				}
				else{
					ProbabilityEstimates probabilities;
					if( !probabilities.Estimate(real_data_stats, ipf_iterations, ipf_precision, num_threads, probs_out.c_str(), probs_in.c_str()) ){
						return 1;
					}
					else{
						probabilities.PrepareResult();

						string org_seq_file;
						auto it = opts_map.find("input");
						if(opts_map.end() == it){
							printInfo << "Reading original sequences from stdin" << std::endl;
						}
						else{
							org_seq_file = it->second.as<string>();
							printInfo << "Reading original sequences from " << org_seq_file << std::endl;
						}

						string destination_file;
						it = opts_map.find("output");
						if(opts_map.end() == it){
							printInfo << "Writing simulated data to stdout" << std::endl;
						}
						else{
							destination_file = it->second.as<string>();
							printInfo << "Writing simulated data to " << destination_file << std::endl;
						}

						auto seed = GetSeed(opts_map);

						if(1.0 != error_mutliplier){
							probabilities.ChangeErrorRate( error_mutliplier );
						}

						Simulator sim;
						if( !sim.SimulateErrorModelOnly( destination_file, org_seq_file, real_data_stats, probabilities, num_threads, seed) ){
							return 1;
						}
					}
				}
			}
		}
		else if("test" == unrecognized_opts.at(0)) {
			if(2 < kVerbosityLevel){
				cerr << " in test mode" << std::endl;
			}

			string usage_str = "Usage: reseq test [options]\n";

			if( general_opts_map.count("help") ) {
				cerr << usage_str;
				cerr << opt_desc_full << std::endl;
			}
			else if(AutoDetectThreads(num_threads, opt_desc_full, usage_str)){
				reseq::AdapterStatsTest::Register();
				reseq::DataStatsTest::Register();
				reseq::FragmentDistributionStatsTest::Register(num_threads);
				reseq::FragmentDuplicationStatsTest::Register();
				reseq::ProbabilityEstimatesTest::Register();
				reseq::ReferenceTest::Register();
				reseq::SeqQualityStatsTest::Register();
				reseq::SimulatorTest::Register();
				reseq::TileStatsTest::Register();
				reseq::VectTest::Register();
				reseq::utilitiesTest::Register();

				bool all_tests_run = false;

				// Setup google test
				unrecognized_opts.at(0) = (*argv)[0]; // Restore proper argv array with program call as first parameter
				int new_argc = unrecognized_opts.size();
				::testing::InitGoogleTest(&new_argc, reinterpret_cast<char**>(&unrecognized_opts.at(0)));

				// Run tests in dependency levels
				string tests_already_run;
				return_code = RunGoogleTests("utilitiesTest.*:VectTest.*", tests_already_run);

				if(0 == return_code){
					return_code = RunGoogleTests("SeqQualityStatsTest.*:ReferenceTest.*:SurroundingTest.*", tests_already_run);

					if(0 == return_code){
						 // Test all remaining classes
						::testing::GTEST_FLAG(filter) = string("-") + tests_already_run;
						return_code = RUN_ALL_TESTS();
						all_tests_run = true;

						if(0 == return_code){
							printSucc << "All tests have succeeded" << std::endl;
						}
					}
				}

				if(0 != return_code && !all_tests_run){
					printErr << "Previous class tests have failed and later classes depend on them in a way that they will definitively fail their test, so the test is aborted here." << std::endl;
				}
			}
		}
		else{
			if(2 < kVerbosityLevel){
				cerr << std::endl;
			}
			printErr << "Unrecognized command: '" << unrecognized_opts.at(0) << "'" << std::endl;
			if(0 < kVerbosityLevel){
				cerr << general_usage << std::endl;
			}
		}
	}

	return return_code;
}
