#ifndef DELPHY_BEASTY_INPUT_H_
#define DELPHY_BEASTY_INPUT_H_

#include <fstream>

#include "phylo_tree.h"
#include "run.h"

namespace delphy {

auto read_beasty_trees(
    std::istream& is,
    std::int64_t min_state,
    std::int64_t every = 1)
    -> std::map<std::int64_t, Phylo_tree>;

auto export_beast_input(
    const Run& run,
    std::string_view beast_version,
    std::ostream& os,
    int64_t chain_length = 10000000,
    int64_t log_every = 1000,
    int64_t tree_every = 1000)
    -> void;

}  // namespace delphy

#endif // DELPHY_BEASTY_INPUT_H_
