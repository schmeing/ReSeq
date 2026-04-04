#ifndef CLI_CONVERT_PROFILE_H
#define CLI_CONVERT_PROFILE_H

#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace reseq::cli {

int RunConvertProfile(const std::vector<std::string>& args, const boost::program_options::variables_map& general_opts,
                      boost::program_options::options_description& opt_desc_full);

} // namespace reseq::cli

#endif // CLI_CONVERT_PROFILE_H
