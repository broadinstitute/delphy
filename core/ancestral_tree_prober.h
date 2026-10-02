#ifndef DELPHY_ANCESTRAL_TREE_PROBER_H_
#define DELPHY_ANCESTRAL_TREE_PROBER_H_

#include <span>

#include "phylo_tree.h"
#include "pop_model.h"

namespace delphy {

// Idea: "mark" a couple of nodes (m_0, ..., m_{k-1}) on a tree, then ask for a probe
// sample at time t the probability p_i that the closest marked ancestor is m_i.  Finally,
// p_k is the probability that none of the marked ancestors are ancestral to the probe
// sample.  Output for marked ancestor index `i` at probe `j` with time `probe_times[j]`
// is stored in `out_values[i * num_t_cells + j]`
auto probe_ancestors_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    std::span<const Node_index> marked_ancestors,
    std::span<const double> probe_times,
    std::span<double> out_values)
    -> void;

}  // namespace delphy

#endif // DELPHY_ANCESTRAL_TREE_PROBER_H_
