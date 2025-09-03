#ifndef DELPHY_CMDLINE_H_
#define DELPHY_CMDLINE_H_

#include <memory>
#include <random>

#include "ctpl_stl.h"

#include "io.h"
#include "phylo_tree.h"
#include "run.h"
#include "sequence_utils.h"

namespace delphy {

// TODO: Should live elsewhere
// Repackage a set of equal-length sequences as a reference + diffs
auto fasta_to_maple(
    const std::vector<Fasta_entry>& in_fasta,
    const std::function<void(int,int)>& progress_hook = [](int,int){},
    const std::function<void(const std::string&, const Sequence_warning&)>& warning_hook
    = default_sequence_warning_hook)
    -> Maple_file;

// TODO: Should live elsewhere
auto build_rough_initial_tree_from_maple(
    Maple_file&& in_maple,
    bool random,
    absl::BitGenRef bitgen,
    const std::function<void(int,int)>& progress_hook = [](int,int){})
    -> Phylo_tree;

struct Processed_cmd_line {
  std::unique_ptr<ctpl::thread_pool> thread_pool;
  std::shared_ptr<Run> run;
  int64_t steps;
  std::optional<std::string> log_filename;
  int64_t log_every;
  std::optional<std::string> trees_filename;
  int64_t tree_every;
  std::optional<std::string> delphy_output_filename;
  std::optional<std::string> delphy_output_metadata;
  int64_t delphy_snapshot_every;
};
auto process_args(int argc, char** argv) -> Processed_cmd_line;

extern bool delphy_invoked_via_cli;
extern std::vector<std::string> delphy_cli_args;

}  // namespace delphy

#endif // DELPHY_CMDLINE_H_
