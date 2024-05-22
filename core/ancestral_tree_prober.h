#ifndef DELPHY_ANCESTRAL_TREE_PROBER_H_
#define DELPHY_ANCESTRAL_TREE_PROBER_H_

#include "tree_prober.h"
#include "phylo_tree.h"

namespace delphy {

// Idea: "mark" a couple of nodes (m_0, ..., m_{k-1}) on a tree, then ask for a probe sample at time t
// the probability p_i that the closest marked ancestor is m_i.  Finally, p_k is the probability that
// none of the marked ancestors are ancestral to the probe sample.
auto probe_ancestors_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    std::span<const Node_index> marked_ancestors,
    double t_start,
    double t_end,
    int num_t_cells)
    -> Staircase_family;

}  // namespace delphy

#endif // DELPHY_ANCESTRAL_TREE_PROBER_H_
