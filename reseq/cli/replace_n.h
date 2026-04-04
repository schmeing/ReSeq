#ifndef CLI_REPLACE_N_H
#define CLI_REPLACE_N_H

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "utilities.hpp"

namespace reseq::cli {

int RunReplaceN(const std::vector<std::string>& args, uintNumThreads num_threads,
                const boost::program_options::variables_map& general_opts,
                boost::program_options::options_description& opt_desc_full);

} // namespace reseq::cli

#endif // CLI_REPLACE_N_H
