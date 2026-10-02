#ifndef DELPHY_WHOLE_TREE_PROBER_H_
#define DELPHY_WHOLE_TREE_PROBER_H_

#include <span>

#include "phylo_tree.h"
#include "pop_model.h"

namespace delphy {

// Calculate the probability that a probe sample at `probe_times[j]` coalesces into the
// tree at branch `i`, i.e., just above every node i (0 <= i < tree.size()).  That value
// is stored in `out_values[i * std::ssize(probe_times) + j]`.  If
// `include_indirect_descendants` is true, then we return instead the related probability
// that the probe coalesces at _or below_ branch `i`.
auto probe_whole_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    std::span<const double> probe_times,
    bool include_indirect_descendants,
    std::span<double> out_values)
    -> void;

}  // namespace delphy

#endif // DELPHY_WHOLE_TREE_PROBER_H_
