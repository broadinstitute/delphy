#include "ancestral_tree_prober.h"

#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <absl/strings/str_format.h>

#include "generic_tree_prober.h"

namespace delphy {

// In this file, "CMA" = "closest marked ancestor"

auto probe_ancestors_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    std::span<const Node_index> marked_ancestors,
    std::span<const double> probe_times,
    std::span<double> out_values)
    -> void {

  if (tree.root == k_no_node) {
    // If we ever wanted to properly handle this pathological case, either the logic of
    // generating bundle events would have to change to always add the -infty event unconditionally,
    // or we'd have to special-case the output for this.  However, there's no universe where
    // calling probe_ancestors_on_tree on an empty tree makes any sense...
    throw std::invalid_argument("tree lacking a root?");
  }

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

  // The following makes it O(1) to find the CMA index of a particular marked node
  auto cma_index_of = Node_map<int>{};
  for (auto i = 0; i != k; ++i) {
    cma_index_of.try_emplace(marked_ancestors[i], i);  // Silently resolves duplicate marked ancestors to the lowest cma_index
  }

  // Traverse tree and record bundle events to associate bundle `i` with CMA `i`.
  // Bundle `k` is everything above all CMAs
  auto bundle_events = std::vector<Bundle_event>{};
  auto cma_stack = std::vector<std::pair<int, Node_index>>{};
  cma_stack.push_back({k, k_no_node});
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Enter `node`.
      // At this point, top of `cma_stack` reflects the state at `parent`.
      auto parent = tree.at(node).parent;

      // Account for the `parent` -> `node` branch
      auto [cur_bundle, _] = cma_stack.back();
      auto t_parent = (node == tree.root) ? -std::numeric_limits<double>::infinity() : tree.at(parent).t;
      bundle_events.push_back({ .t = t_parent, .adding = true, .target_bundle = cur_bundle });
      bundle_events.push_back({ .t = tree.at(node).t, .adding = false, .target_bundle = cur_bundle });

      // Now moving downstream of node.  If it's a marked ancestor, push current bundle and switch to a new bundle
      if (auto it = cma_index_of.find(node); it != cma_index_of.end()) {
        auto new_bundle = it->second;
        cma_stack.push_back({new_bundle, node});
      }
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`.
      // At this point, top of `cma_stack` reflects the state at `node`.
      auto [_, bundle_start_node] = cma_stack.back();
      if (bundle_start_node == node) {
        // Done with this bundle
        cma_stack.pop_back();
      }
    }
  }

  // Run the tree probe
  generic_probe_tree(bundle_events, k+1, pop_model, probe_times, out_values);
}

}  // namespace delphy
