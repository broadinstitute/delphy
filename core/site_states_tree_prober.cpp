#include "site_states_tree_prober.h"

#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <absl/log/check.h>
#include <absl/strings/str_format.h>

#include "generic_tree_prober.h"

namespace delphy {

auto probe_site_states_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    Site_index site,
    std::span<const double> probe_times,
    std::span<double> out_values)
    -> void {

  if (tree.root == k_no_node) {
    // If we ever wanted to properly handle this pathological case, either the logic of
    // generating bundle events would have to change to always add the -infty event unconditionally,
    // or we'd have to special-case the output for this.  However, there's no universe where
    // calling probe_site_states_on_tree on an empty tree makes any sense...
    throw std::invalid_argument("tree lacking a root?");
  }

  auto L = tree.num_sites();
  if (site < 0 || site >= L) {
    throw std::out_of_range(absl::StrFormat("Site %d is outside the valid range [1, %d]", site+1, L));
  }

  // NOTE: we really should check if `site` is missing even at the root, in which case, p_initial should
  // match pi_a and the probabilities will not change across time.  It's not worth handling this edge
  // case properly at the moment.
  auto state_at_root = tree.ref_sequence[site];
  for (const auto& m : tree.at_root().mutations) {
    if (m.site == site) {
      CHECK_EQ(state_at_root, m.from);
      state_at_root = m.to;
    }
  }
  
  // Traverse tree and record bundle events to associate bundle `a` with state `a`.
  auto num_bundles = k_num_real_seq_letters;
  auto bundle_events = std::vector<Bundle_event>{};
  auto cur_state = state_at_root;
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Enter `node`.
      // At this point, cur_state reflects the state at `parent`.
      auto parent = tree.at(node).parent;

      // Account for the `parent` -> `node` branch
      auto t_parent = (node == tree.root) ? -std::numeric_limits<double>::infinity() : tree.at(parent).t;
      bundle_events.push_back({ .t = t_parent, .adding = true, .target_bundle = index_of(cur_state) });
      if (node != tree.root) {  // Mutations above the root are just deltas from the ref_sequence
        for (const auto& m : tree.at(node).mutations) {
          if (m.site == site) {
            CHECK_EQ(cur_state, m.from);
            bundle_events.push_back({ .t = m.t, .adding = false, .target_bundle = index_of(m.from) });
            bundle_events.push_back({ .t = m.t, .adding = true, .target_bundle = index_of(m.to) });
            cur_state = m.to;
          }
        }
      }
      bundle_events.push_back({ .t = tree.at(node).t, .adding = false, .target_bundle = index_of(cur_state) });
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`.
      // At this point, cur_state reflects the state at `node`.
      if (node != tree.root) {  // Mutations above the root are just deltas from the ref_sequence
        for (const auto& m : tree.at(node).mutations | std::views::reverse) {
          if (m.site == site) {
            CHECK_EQ(m.to, cur_state);
            cur_state = m.from;
          }
        }
      }
    }
  }

  // Run the tree probe
  generic_probe_tree(bundle_events, num_bundles, pop_model, probe_times, out_values);
}

}  // namespace delphy
