#ifndef DELPHY_SITE_STATES_TREE_PROBER_H_
#define DELPHY_SITE_STATES_TREE_PROBER_H_

#include <span>

#include "phylo_tree.h"
#include "pop_model.h"

namespace delphy {

auto probe_site_states_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    Site_index site,
    std::span<const double> probe_times,
    std::span<double> out_values)
    -> void;

}  // namespace delphy

#endif // DELPHY_SITE_STATES_TREE_PROBER_H_
