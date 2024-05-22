#ifndef DELPHY_SITE_STATES_TREE_PROBER_H_
#define DELPHY_SITE_STATES_TREE_PROBER_H_

#include "tree_prober.h"
#include "phylo_tree.h"

namespace delphy {

auto probe_site_states_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    Site_index site,
    double t_start,
    double t_end,
    int num_t_cells)
    -> Staircase_family;

}  // namespace delphy

#endif // DELPHY_SITE_STATES_TREE_PROBER_H_
