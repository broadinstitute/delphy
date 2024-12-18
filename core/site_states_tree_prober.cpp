#include "site_states_tree_prober.h"

namespace delphy {

static auto probe_site_states_on_tree_helper(
    const Phylo_tree& tree,
    Node_index node,
    Site_index site,
    Real_seq_letter state,
    Staircase_family& branch_counts_by_state)
    -> void {
  
  // Accumulate counts for branch from parent to node
  // If there's a mutation on that branch, adjust `state` to match the endpoint of that branch
  if (node != tree.root) {
    auto parent = tree.at(node).parent;
    if (auto it = std::ranges::find_if(tree.at(node).mutations, [site](const auto& mut) { return mut.site == site; });
        it != tree.at(node).mutations.end()) {
      
      // A mutation along this branch on the site we're looking at
      const auto& mut = *it;
      
      assert(state == mut.from);
      add_trapezoid(branch_counts_by_state[index_of(state)], tree.at(parent).t, tree.at(node).t, 1.0, 0.0);
      
      state = mut.to;
      add_trapezoid(branch_counts_by_state[index_of(state)], tree.at(parent).t, tree.at(node).t, 0.0, 1.0);
      
    } else {
      add_boxcar(branch_counts_by_state[index_of(state)], tree.at(parent).t, tree.at(node).t, 1.0);
    }
  }

  // Recurse through children
  for (const auto& child : tree.at(node).children) {
    probe_site_states_on_tree_helper(tree, child, site, state, branch_counts_by_state);
  }
}

auto probe_site_states_on_tree(
    const Phylo_tree& tree,
    const Pop_model& pop_model,
    Site_index site,
    double t_start,
    double t_end,
    int num_t_cells)
    -> Staircase_family {

  auto L = tree.num_sites();
  if (site < 0 || site >= L) {
    throw std::out_of_range(absl::StrFormat("Site %d is outside the valid range [1, %d]", site+1, L));
  }

  // Initial probabilities match root sequence if t_start <= t_root
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
  
  auto state_at_root = tree.ref_sequence[site];
  for (const auto& m : tree.at_root().mutations) {
    if (m.site == site) {
      CHECK_EQ(state_at_root, m.from);
      state_at_root = m.to;
    }
  }

  // NOTE: we really should check if `site` is missing even at the root, in which case, p_initial should
  // match pi_a and the probabilities will not change across time.  It's not worth handling this edge
  // case properly at the moment.
  auto p_initial = std::vector<double>(k_num_real_seq_letters, 0.0);
  p_initial[index_of(state_at_root)] = 1.0;
  
  // Accumulate counts of branches where site has states A, C, G or T
  auto branch_counts_by_state = Staircase_family{k_num_real_seq_letters, real_t_start, t_end, num_t_cells};
  if (tree.root != k_no_node) {
    probe_site_states_on_tree_helper(tree, tree.root, site, state_at_root, branch_counts_by_state);
  }

  auto prober = Tree_prober{branch_counts_by_state, cells_to_skip, pop_model, std::move(p_initial)};
  
  return prober.p();
}

}  // namespace delphy
