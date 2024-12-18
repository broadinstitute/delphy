#include "ancestral_tree_prober.h"

namespace delphy {

// In this file, "CMA" = "closest marked ancestor"

static auto probe_ancestors_on_tree_helper(const Phylo_tree& tree,
                                           Node_index node,
                                           std::span<const Node_index> marked_ancestors,
                                           int cma_index,
                                           Staircase_family& branch_counts_by_cma) -> void {
  // Accumulate counts for branch from parent to node
  if (node != tree.root && cma_index >= 0) {
    add_boxcar(branch_counts_by_cma[cma_index], tree.at_parent_of(node).t, tree.at(node).t, 1.0);
  }
  
  // If `node` is a marked ancestor, adjust cma_index for children
  auto it = std::ranges::find(marked_ancestors, node);
  if (it != marked_ancestors.end()) {
    cma_index = static_cast<int>(std::distance(marked_ancestors.begin(), it));
  }
  
  // Recurse through children
  for (const auto& child : tree.at(node).children) {
    probe_ancestors_on_tree_helper(tree, child, marked_ancestors, cma_index, branch_counts_by_cma);
  }
}

auto probe_ancestors_on_tree(const Phylo_tree& tree,
                             const Pop_model& pop_model,
                             std::span<const Node_index> marked_ancestors,
                             double t_start,
                             double t_end,
                             int num_t_cells) -> Staircase_family {
  for (const auto& node : marked_ancestors) {
    // Including k_no_node in marked_ancestors is useful in case
    // we're iterating over all base trees in an MCC, and a particular
    // base tree doesn't have a particular clade.
    auto node_valid = (node == k_no_node) || (node >= 0 && node < std::ssize(tree));
    if (not node_valid) {
      throw std::out_of_range(absl::StrFormat(
          "Node %d is neither `none` (%d) nor inside the valid range [0, %d)",
          node, k_no_node, std::ssize(tree)));
    }
  }

  auto k = static_cast<int>(std::ssize(marked_ancestors));
  
  // Initial probabilities are all 0.0 if t_start <= t_root
  // If t_start > t_root, we add cells at the beginning until we reach past the root time,
  // but mark those for ignoring (there's a much better way to do this properly, but
  // this crude scheme will do for now)
  auto real_t_start = t_start;
  auto cells_to_skip = 0;
  if (tree.root != k_no_node && t_start > tree.at_root().t) {
    auto cell_size = (t_end - t_start) / num_t_cells;
    while (real_t_start > tree.at_root().t) {
      real_t_start -= cell_size;
      ++num_t_cells;
      ++cells_to_skip;
    }
  }
  
  auto p_initial = std::vector<double>(k+1, 0.0);
  p_initial[k] = 1.0;
  
  // Accumulate counts of branches where CMA is `i`
  auto branch_counts_by_cma = Staircase_family{k+1, real_t_start, t_end, num_t_cells};
  if (tree.root != k_no_node) {
    probe_ancestors_on_tree_helper(tree, tree.root, marked_ancestors, k, branch_counts_by_cma);
  }
  
  auto prober = Tree_prober{branch_counts_by_cma, cells_to_skip, pop_model, std::move(p_initial)};
  
  return prober.p();
}

}  // namespace delphy
