#include "whole_tree_prober.h"

#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <absl/strings/str_format.h>

#include "generic_tree_prober.h"

namespace delphy {

auto probe_whole_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    std::span<const double> probe_times,
    bool include_indirect_descendants,
    std::span<double> out_values)
    -> void {

  if (tree.root == k_no_node) {
    // If we ever wanted to properly handle this pathological case, either the logic of
    // generating bundle events would have to change to always add the -infty event unconditionally,
    // or we'd have to special-case the output for this.  However, there's no universe where
    // calling probe_whole_tree on an empty tree makes any sense...
    throw std::invalid_argument("tree lacking a root?");
  }

  // Traverse tree and record bundle events to associate bundle `i` with branch `i`.
  auto num_bundles = tree.size();
  auto bundle_events = std::vector<Bundle_event>{};
  for (const auto& node : index_order_traversal(tree)) {
    auto t_parent = (node == tree.root) ? -std::numeric_limits<double>::infinity() : tree.at_parent_of(node).t;
    bundle_events.push_back({ .t = t_parent, .adding = true, .target_bundle = node });
    bundle_events.push_back({ .t = tree.at(node).t, .adding = false, .target_bundle = node });
  }

  // Run the tree probe
  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, out_values);

  // Aggregate indirect descendants if requested
  if (include_indirect_descendants) {
    auto num_probe_times = static_cast<int>(std::ssize(probe_times));
    auto results = estd::View_2d{out_values, num_bundles, num_probe_times};
    for (const auto& node : post_order_traversal(tree)) {
      for (const auto& child : tree.at(node).children) {  // Note: tips have no children!
        for (auto i = 0; i != num_probe_times; ++i) {
          results(node, i) += results(child, i);
        }
      }
    }
  }
}

}  // namespace delphy
