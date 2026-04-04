#ifndef CLI_COMMON_H
#define CLI_COMMON_H

#include <string>

#include <boost/program_options.hpp>

#include "utilities.hpp"

namespace reseq::cli {

bool AutoDetectThreads(uintNumThreads& num_threads, const boost::program_options::options_description& opt_desc,
                       const std::string& usage_str);

uintSeed GetSeed(const boost::program_options::variables_map& opts_map);

void GetProbsOut(std::string& probs_out, const std::string& fallback_out,
                 const boost::program_options::variables_map& opts_map);

void PrepareProbabilityEstimation(std::string& probs_in, std::string& probs_out, const std::string& standard_probs_out,
                                  bool loaded_stats, const boost::program_options::variables_map& opts_map);

} // namespace reseq::cli

#endif // CLI_COMMON_H
