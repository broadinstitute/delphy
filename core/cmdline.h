#ifndef DELPHY_CMDLINE_H_
#define DELPHY_CMDLINE_H_

#include <memory>
#include <random>

#include "ctpl_stl.h"

#include "io.h"
#include "phylo_tree.h"
#include "run.h"

namespace delphy {

// TODO: Should live elsewhere
auto extract_date_from_fasta_id(const std::string_view id) -> std::optional<double>;
auto build_rough_initial_tree_from_fasta(
    const std::vector<Fasta_entry>& in_fasta,
    bool random,
    absl::BitGenRef bitgen)
    -> Phylo_tree;

struct Processed_cmd_line {
  std::unique_ptr<std::mt19937> prng;
  std::unique_ptr<ctpl::thread_pool> thread_pool;
  std::shared_ptr<Run> run;
  int64_t steps;
  std::optional<std::string> log_filename;
  int64_t log_every;
  std::optional<std::string> trees_filename;
  int64_t tree_every;
  bool alpha_move_enabled;
  bool mu_move_enabled;
  double init_mu;
  bool mpox_hack_enabled;
};
auto process_args(int argc, char** argv) -> Processed_cmd_line;

extern bool delphy_invoked_via_cli;
extern std::vector<std::string> delphy_cli_args;

}  // namespace delphy

#endif // DELPHY_CMDLINE_H_
