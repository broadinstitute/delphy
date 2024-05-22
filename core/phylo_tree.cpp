#include "phylo_tree.h"

#include "sequence_overlay.h"
#include "phylo_tree_calc.h"
#include "spr_study.h"
#include "dates.h"

namespace delphy {

auto operator<<(std::ostream& os, const Phylo_tree& tree) -> std::ostream& {
  os << "Ref sequence: '" << absl::StrJoin(tree.ref_sequence, "", absl::StreamFormatter()) << "'\n";
  for (const auto& node : index_order_traversal(tree)) {
    os << (node == tree.root ? '*' : ' ') << absl::StreamFormat("[%3d] ", node) << tree.at(node) << "\n";
  }
  return os;
}

auto assert_mutation_consistency(const Phylo_tree& tree, bool force) -> void {
  if (estd::is_debug_enabled || force) {

    auto cur_seq = Sequence_overlay{tree.ref_sequence};
    for (const auto& [node, children_so_far] : traversal(tree)) {
      if (children_so_far == 0) {
        // Entering `node`

        // Check mutations have reasonable to/from and times
        auto min_t = node != tree.root ? tree.at_parent_of(node).t : -std::numeric_limits<double>::max();
        for (const auto& m : tree.at(node).mutations) {
          CHECK_NE(m.from, m.to);
          CHECK_LE(0, m.site);
          CHECK_LT(m.site, tree.num_sites());
          CHECK_LE(min_t, m.t);
          CHECK_LE(m.t, tree.at(node).t);
          min_t = m.t;  // Mutations should be sorted in ascending order of time
        }
        
        // Apply mutations and check their from states
        for (const auto& m : tree.at(node).mutations) {
          CHECK_EQ(m.from, cur_seq[m.site]);
          cur_seq[m.site] = m.to;
        }
      }

      if (children_so_far == std::ssize(tree.at(node).children)) {
        // Exiting `node`
        // Revert mutations and check their to states
        for (const auto& m : tree.at(node).mutations | std::views::reverse) {
          CHECK_EQ(m.to, cur_seq[m.site])
              << m << ' ' << absl::StrJoin(cur_seq, "", absl::StreamFormatter());
          cur_seq[m.site] = m.from;
        }
      }
    }
  }
}

auto assert_missation_consistency(const Phylo_tree& tree, bool force) -> void {
  if (estd::is_debug_enabled || force) {
    auto cur_seq = Sequence_overlay{tree.ref_sequence};
    auto cur_missing_sites = absl::flat_hash_set<Site_index>{};
    for (const auto& [node, children_so_far] : traversal(tree)) {
      if (children_so_far == 0) {
        // Entering `node`

        // Record sites missing below the beginning of branch ending at `node` and check their from states
        for (const auto& [mi_site, mi_from] : tree.at(node).missations.slow_elements(tree.ref_sequence)) {
          CHECK(not cur_missing_sites.contains(mi_site));
          CHECK_EQ(mi_from, cur_seq[mi_site]);
          cur_missing_sites.insert(mi_site);
        }
        
        // Apply mutations and check their from() states
        for (const auto& m : tree.at(node).mutations) {
          CHECK_EQ(m.from, cur_seq[m.site]);
          cur_seq[m.site] = m.to;
        }
      }

      if (children_so_far == std::ssize(tree.at(node).children)) {
        // Exiting `node`

        // Flag missations present on all child branches (should be recursively merged into an upstream missation)
        if (not tree.at(node).children.empty()) {
          CHECK_EQ(std::ssize(tree.at(node).children), 2);
          auto common_sites = intersect_interval_sets(tree.at(tree.at(node).left_child()).missations.intervals,
                                                      tree.at(tree.at(node).right_child()).missations.intervals);
          CHECK(common_sites.empty()) << absl::StreamFormat(
              "Missations on both left and right children of node %d: %s",
              node, absl::FormatStreamed(common_sites));
        }
        
        // Flag mutations on missing sites
        for (const auto& m : tree.at(node).mutations) {
          CHECK(not cur_missing_sites.contains(m.site));
        }
        
        // Revert mutations and check their to() states
        for (const auto& m : tree.at(node).mutations | std::views::reverse) {
          CHECK_EQ(m.to, cur_seq[m.site]);
          cur_seq[m.site] = m.from;
        }
        
        // Revert missations
        for (const auto& mi_site : tree.at(node).missations.intervals.slow_elements()) {
          CHECK(cur_missing_sites.contains(mi_site));
          cur_missing_sites.erase(mi_site);
        }
      }
    }
  }
}

auto assert_phylo_tree_integrity(const Phylo_tree& tree, bool force) -> void {
  if (estd::is_debug_enabled || force) {
    assert_tree_integrity(tree, force);

    assert_mutation_consistency(tree, force);
    assert_missation_consistency(tree, force);
    // FIXME: Commented out for node since some tests still have zero-length branches in them
    // for (const auto& node : pre_order_traversal(tree)) {
    //   if (tree.at(node).is_inner_node()) {
    //     for (const auto& child : tree.at(node).children) {
    //       CHECK_LT(tree.at(node).t, tree.at(child).t);  // Note: Forbid zero-length branches!
    //     }
    //   }
    // }
  }
}

auto find_MRCA_of(const Phylo_tree& tree, Node_index P, Node_index Q) -> Node_index {
  if (P == k_no_node) { return P; }
  if (Q == k_no_node) { return Q; }

  while (P != Q) {
    auto t_P = tree.at(P).t;
    auto t_Q = tree.at(Q).t;
    
    if (t_P > t_Q) {
      P = tree.at(P).parent;
      CHECK_NE(P, k_no_node);
      
    } else if (t_P < t_Q) {
      Q = tree.at(Q).parent;
      CHECK_NE(Q, k_no_node);
      
    } else if (tree.at(P).is_tip()) {
      // P is certainly not the MRCA
      P = tree.at(P).parent;
      CHECK_NE(P, k_no_node);
      
    } else if (tree.at(Q).is_tip()) {
      // Q is certainly not the MRCA
      Q = tree.at(Q).parent;
      CHECK_NE(Q, k_no_node);
      
    } else {
      // Times are equal and nodes are distinct inner nodes: this case is very rare, so we can
      // afford to do expensive things here (inner node times are random numbers almost always,
      // so it's nearly impossible for two of them to be equal just by chance).
      //
      // We'll enumerate all ancestors of both P and Q until the root.  These lists share a common prefix
      // of all common ancestors of P and Q.  By successively popping nodes off these lists until we hit
      // a difference, we'll end up finding the most recent common ancestor.
      //
      // If this case weren't so rare, we could just enumerate ancestors of P and Q until we hit one with
      // a time distinct from t_P == t_Q (or we hit the root), and then decide if we found an MRCA among those
      // smaller ancestor lists or we need to loop back.  But since this case is rare, I think it's better
      // to use the following straightforwardly correct but inefficient approach.

      auto ancestors_P = std::vector<Node_index>{};
      for (auto cur = P; cur != k_no_node; cur = tree.at(cur).parent) {
        ancestors_P.push_back(cur);
      }
      CHECK_EQ(ancestors_P.front(), P);
      CHECK_EQ(ancestors_P.back(), tree.root);
      
      auto ancestors_Q = std::vector<Node_index>{};
      for (auto cur = Q; cur != k_no_node; cur = tree.at(cur).parent) {
        ancestors_Q.push_back(cur);
      }
      CHECK_EQ(ancestors_Q.front(), Q);
      CHECK_EQ(ancestors_Q.back(), tree.root);

      auto candidate_mrca = tree.root;
      while (not ancestors_P.empty() && not ancestors_Q.empty() && ancestors_P.back() == ancestors_Q.back()) {
        candidate_mrca = ancestors_P.back();
        ancestors_P.pop_back();
        ancestors_Q.pop_back();
      }

      return candidate_mrca;
    }
  }

  CHECK_NE(P, k_no_node);
  return P;  // == Q
}

auto find_MRCA_of(const Phylo_tree& tree, Phylo_tree_loc p, Phylo_tree_loc q) -> Phylo_tree_loc {
  if (p.branch == q.branch) {
    return p.t < q.t ? p : q;
  }

  auto A = find_MRCA_of(tree, p.branch, q.branch);
  return {A, std::min(p.t, std::min(q.t, tree.at(A).t))};
}

auto descends_from(const Phylo_tree& tree, Node_index X, Node_index A) -> bool {
  if (A == k_no_node) { return true; }
  for (auto cur = X; cur != k_no_node; cur = tree.at(cur).parent) {
    if (cur == A) { return true; }
    if (tree.at(cur).t < tree.at(A).t) { return false; }
  }
  return false;  // Went through the root without hitting A!
}

auto descends_from(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc a) -> bool {
  if (x.branch == a.branch) {
    return a.t <= x.t;
  } else {
    return descends_from(tree, x.branch, a.branch);
  }
}

auto rereference_to_root_sequence(Phylo_tree& tree) -> void {
 // Change ref sequence to match sequence at root
  for (const auto& mut : tree.at_root().mutations) {
    tree.ref_sequence[mut.site] = mut.to;
  }
  for (const auto& node : index_order_traversal(tree)) {
    if (not tree.at(node).missations.empty()) {
      for (const auto& mut : tree.at_root().mutations) {
        tree.at(node).missations.ref_seq_changed(mut.site, mut.from, mut.to);
      }
    }
  }
  tree.at_root().mutations = {};
}

auto push_all_mutations_and_missations_to_tips(Phylo_tree& tree) -> void {
  // Build up delta from ref sequence to each tip, and replace mutations and missations at tip with the whole shebang
  const auto& ref_seq = tree.ref_sequence;
  auto cur_seq = Sequence_overlay{ref_seq};
  auto cur_missing_sites = Interval_set<>{};
  auto scratch_sites = Interval_set<>{};
  auto new_mutations = Mutation_list<>{};
  auto new_missations = Missation_map<>{};
  assert_phylo_tree_integrity(tree, true);
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Apply mutations
      for (const auto& m : tree.at(node).mutations) {
        CHECK_EQ(m.from, cur_seq[m.site]);
        cur_seq[m.site] = m.to;
      }
      
      // Record sites missing below the beginning of branch ending at `node`
      CHECK(not interval_sets_intersect(tree.at(node).missations.intervals, cur_missing_sites));
      merge_interval_sets(scratch_sites, tree.at(node).missations.intervals, cur_missing_sites);
      swap(cur_missing_sites, scratch_sites);
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`

      // Materialize accumulated mutations at tips (with times equal to the tip time) 
      if (tree.at(node).is_tip()) {
        new_mutations.clear();
        for (const auto [l, b] : cur_seq.deltas()) {
          auto a = ref_seq[l];
          new_mutations.push_back(Mutation{a, l, b, tree.at(node).t});
        }
        sort_mutations(new_mutations);
      }

      // Materialize accumulated missations at tips
      if (tree.at(node).is_tip()) {
        new_missations.clear();
        new_missations.intervals = cur_missing_sites;
        new_missations.from_states.clear();  // they're all equal to the reference sequence
      }
      
      // Revert mutations
      for (const auto& m : tree.at(node).mutations | std::views::reverse) {
        CHECK_EQ(m.to, cur_seq[m.site]);
        cur_seq[m.site] = m.from;
      }

      // Revert missations
      subtract_interval_sets(scratch_sites, cur_missing_sites, tree.at(node).missations.intervals);
      swap(cur_missing_sites, scratch_sites);
      
      // Write over mutations and missations now that we no longer need them
      if (tree.at(node).is_tip()) {
        tree.at(node).mutations = std::move(new_mutations);
        tree.at(node).missations = std::move(new_missations);
      } else {
        tree.at(node).mutations = {};
        tree.at(node).missations = {};
      }
    }
  }
}

auto fix_up_missations(Phylo_tree& tree) -> void {

  const auto& ref_seq = tree.ref_sequence;
  
  // First, bubble up any common missations from the leaves to the root
  for (const auto& node : post_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      auto L = tree.at(node).left_child();
      auto R = tree.at(node).right_child();
      if (interval_sets_intersect(tree.at(L).missations.intervals, tree.at(R).missations.intervals)) {
        if (tree.at(node).missations.empty()) {
          factor_out_common_missations(tree.at(L).missations, tree.at(R).missations, tree.at(node).missations);
        } else {
          auto extra_node_missations = Missation_map<>{};
          factor_out_common_missations(tree.at(L).missations, tree.at(R).missations, extra_node_missations);
          tree.at(node).missations = merge_missations_nondestructively(tree.at(node).missations, extra_node_missations);
        }
      }
    }
  }

  // Remove any redundant downstream mutations
  auto cur_missing_sites = absl::flat_hash_set<Site_index>{};
  auto cur_seq = Sequence_overlay{ref_seq};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    auto needs_fixing = false;
    if (children_so_far == 0) {
      // Apply mutations
      for (const auto& m : tree.at(node).mutations) {
        CHECK_EQ(m.from, cur_seq[m.site]);
        cur_seq[m.site] = m.to;
      }
      
      // Record sites missing below the beginning of branch ending at `node`
      for (const auto& mi_site : tree.at(node).missations.intervals.slow_elements()) {
        if (cur_missing_sites.contains(mi_site)) { needs_fixing = true; }
        cur_missing_sites.insert(mi_site);
      }
      
      if (needs_fixing) {
        auto new_missations = Missation_map<>{};
        for (auto l : cur_missing_sites) {
          if (node == tree.root || not is_site_missing_at(tree, tree.at(node).parent, l)) {
            new_missations.insert(Missation{l, cur_seq[l]}, ref_seq);
          }
        }
        tree.at(node).missations = std::move(new_missations);
      }
    }
        
    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`

      // Revert mutations
      for (const auto& m : tree.at(node).mutations | std::views::reverse) {
        CHECK_EQ(m.to, cur_seq[m.site]);
        cur_seq[m.site] = m.from;
      }

      // Revert missations
      for (const auto& mi_site : tree.at(node).missations.intervals.slow_elements()) {
        CHECK(cur_missing_sites.contains(mi_site));
        cur_missing_sites.erase(mi_site);
      }
    }
  }

  // Now reconstruct `from_states`
  cur_seq = Sequence_overlay{ref_seq};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      
      // Reconstruct from_states given cur_seq
      tree.at(node).missations.from_states.clear();
      for (const auto& [l, s] : cur_seq.deltas()) {
        if (tree.at(node).missations.contains(l)) {
          tree.at(node).missations.from_states.insert({l, s});
        }
      }

      // Remove mutations on missing sites
      std::erase_if(tree.at(node).mutations, [&](const auto& m) { return is_site_missing_at(tree, node, m.site); });
      
      // Apply mutations
      for (const auto& m : tree.at(node).mutations) {
        CHECK_EQ(m.from, cur_seq[m.site]);
        cur_seq[m.site] = m.to;
      }
    }
    
    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      
      // Revert mutations
      for (const auto& m : tree.at(node).mutations | std::views::reverse) {
        CHECK_EQ(m.to, cur_seq[m.site]);
        cur_seq[m.site] = m.from;
      }
    }
  }
}

auto extract_nodes(const Phylo_tree& tree) -> extract_nodes_results {
  auto result = extract_nodes_results{};
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      result.tips.push_back(node);
    } else {
      result.inner_nodes.push_back(node);
    }
  }
  return result;
}

auto rewire_tree_through_sequential_accretion(
    Phylo_tree& tree,
    const std::vector<Node_index>& tips,
    double inner_node_offset)
    -> void {
  
  // Validation
  auto num_tips = std::ssize(tree) == 0 ? 0 : (std::ssize(tree) + 1) / 2;
  if (std::ssize(tips) != num_tips) {
    throw std::invalid_argument(absl::StrFormat(
        "Mismatch between actual number of tips (%d) and provided tips (%d)", num_tips, std::ssize(tips)));
  } else {
    auto seen_tips = Node_set{};
    for (const auto& tip : tips) {
      if (tip < 0 || tip >= std::ssize(tree)) {
        throw std::out_of_range(absl::StrFormat(
            "Tip index %d out of valid range [0, %d)", tip, std::ssize(tree)));
      }
      if (seen_tips.contains(tip)) {
        throw std::invalid_argument(absl::StrFormat("Duplicate tip index %d", tip));
      } else {
        seen_tips.insert(tip);
      }
      if (not tree.at(tip).is_tip()) {
        throw std::invalid_argument(absl::StrFormat("Node %d is not a tip", tip));
      }
    }
  }
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      if (not tree.at(node).mutations.empty()) {
        throw std::invalid_argument(absl::StrFormat(
            "Inner node unexpectedly has mutations: %s", absl::FormatStreamed(node)));
      }
      if (not tree.at(node).missations.empty()) {
        throw std::invalid_argument(absl::StrFormat(
            "Inner node unexpectedly has missations: %s", absl::FormatStreamed(node)));
      }
    }
  }

  // Off to the races!
  auto inner_nodes = Node_set{};
  for (auto index = 0; index != std::ssize(tree); ++index) {
    inner_nodes.insert(index);
  }
  for (auto tip : tips) {
    inner_nodes.erase(tip);
  }
  
  auto inner_node_it = std::ranges::begin(inner_nodes);
  [[maybe_unused]] auto inner_node_end = std::ranges::end(inner_nodes);
  tree.root = k_no_node;
  for (auto next_tip : tips) {
    if (tree.root == k_no_node) {
      tree.root = next_tip;
      tree.at(next_tip).parent = k_no_node;
    } else {
      assert(inner_node_it != inner_node_end);
      auto new_root = *inner_node_it;
      ++inner_node_it;
      auto new_left = tree.root;
      auto new_right = next_tip;

      tree.root = new_root;
      tree.at(new_root).name = "";
      tree.at(new_root).t = std::min(tree.at(new_left).t, tree.at(new_right).t) - inner_node_offset;
      tree.at(new_root).parent = k_no_node;
      tree.at(new_root).children = {new_left, new_right};
      tree.at(new_left).parent = new_root;
      tree.at(new_right).parent = new_root;
    }
  }
}

auto randomize_mutation_times(Phylo_tree& tree, absl::BitGenRef bitgen, Scratch_space& scratch) -> void {
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root) {
      auto new_mutations = randomize_branch_mutation_times(tree, node, bitgen, scratch);
      tree.at(node).mutations.assign(new_mutations.begin(), new_mutations.end());
    }
  }
}

auto randomize_branch_mutation_times(
    Phylo_tree& tree,
    Branch_index X,
    absl::BitGenRef bitgen,
    Scratch_space& scratch)
    -> Scratch_vector<Mutation> {

  const auto& old_mutations = tree.at(X).mutations;
  
  if (X == tree.root) {
    // Mutations above the root are not really mutations
    return Scratch_vector<Mutation>(old_mutations.begin(), old_mutations.end(), scratch);
  }
  
  auto P = tree.at(X).parent;
  auto t_X = tree.at(X).t;
  auto t_P = tree.at(P).t;
  
  // See if branch has multiply mutated sites (it rarely does)
  auto mutated_sites = Scratch_flat_hash_set<Site_index>{scratch};
  mutated_sites.reserve(std::ssize(old_mutations));
  auto complicated = false;
  for (auto& m : old_mutations) {
    auto [it, inserted] = mutated_sites.insert(m.site);
    if (not inserted) {
      complicated = true;  // This site has 2+ mutations
    }
  }
  
  auto new_mutations = Scratch_vector<Mutation>{scratch};
  new_mutations.reserve(std::ssize(old_mutations));
  if (not complicated) {
    // Pick new mutation times uniformly
    for (const auto& m : old_mutations) {
      new_mutations.push_back(Mutation{m.from, m.site, m.to,
          absl::Uniform(absl::IntervalOpenClosed, bitgen, t_P, t_X)});
    }
  } else {
    // We have to treat each site separately, and propose new mutation times
    // uniformly but in the same order.
    // This path is intentionally written for clarity instead of speed because it is rare.
    auto new_site_mut_times = Scratch_vector<double>{scratch};
    new_site_mut_times.reserve(std::ssize(old_mutations));
    for (auto l : mutated_sites) {
      new_site_mut_times.clear();
      for (const auto& m : old_mutations) {
        if (m.site == l) {
          new_site_mut_times.push_back(absl::Uniform(absl::IntervalOpenClosed, bitgen, t_P, t_X));
        }
      }
      std::ranges::sort(new_site_mut_times);
      auto it = new_site_mut_times.begin();
      for (const auto& m : old_mutations) {
        if (m.site == l) {
          CHECK(it != new_site_mut_times.end());
          new_mutations.push_back(Mutation{m.from, m.site, m.to, *it});
          ++it;
        }
      }
      CHECK(it == new_site_mut_times.end());
    }
  }
  
  sort_mutations(new_mutations);
  
  return new_mutations;
}

auto randomize_tree(Phylo_tree& tree, absl::BitGenRef bitgen) -> void {
  
  push_all_mutations_and_missations_to_tips(tree);

  auto [tips, inner_nodes] = extract_nodes(tree);
  std::ranges::shuffle(tips, bitgen);

  auto num_inner_nodes = (std::ssize(tree) - 1) / 2;
  auto inner_node_offset = 50.0 / num_inner_nodes;
  rewire_tree_through_sequential_accretion(tree, tips, inner_node_offset);

  fix_up_missations(tree);

  auto scratch = Scratch_space{};
  randomize_mutation_times(tree, bitgen, scratch);

  rereference_to_root_sequence(tree);

  assert_phylo_tree_integrity(tree);
}

auto build_random_tree(
    Real_sequence ref_sequence,
    std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen)
    -> Phylo_tree {
  
  // Sanity-check inputs
  for (auto tip = 0; tip != std::ssize(tip_descs); ++tip) {
    const auto& tip_desc = tip_descs[tip];

    // Mutations
    for (const auto& m : tip_desc.seq_deltas) {
      if (m.site < 0 || m.site >= std::ssize(ref_sequence)) {
        throw std::out_of_range(absl::StrFormat(
            "In tip %d (%s), site of mutation %s is outside valid range [0,%d)",
            tip, tip_desc.name, absl::FormatStreamed(m), std::ssize(ref_sequence)));
      }
      if (auto ref_state = ref_sequence[m.site]; m.from != ref_state) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), mutation %s cannot be applied because ref sequence has state %c there",
            tip, tip_desc.name, absl::FormatStreamed(m), to_char(ref_state)));
      }
      if (m.from == m.to) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), mutation %s has equal 'from' and 'to' states",
            tip, tip_desc.name, absl::FormatStreamed(m)));
      }
    }

    // Missations
    for (const auto& [mi_start, mi_end] : tip_desc.missations.intervals) {
      if (mi_start < 0 || mi_start >= std::ssize(ref_sequence) || mi_end < 0 || mi_end >= (std::ssize(ref_sequence)+1)) {
        throw std::out_of_range(absl::StrFormat(
            "In tip %d (%s), missation interval [%d,%d) applies to a site outside valid range [0,%d)",
            tip, tip_desc.name, mi_start+1, mi_end+1, std::ssize(ref_sequence)));
      }
    }

    // Mutations vs missations
    for (const auto& m : tip_desc.seq_deltas) {
      if (tip_desc.missations.contains(m.site)) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), site %d can't both have a missation and a mutation (%s)",
            tip, tip_desc.name, m.site+1, absl::FormatStreamed(m)));
      }
    }
  }

  if (tip_descs.empty()) {
    auto tree = Phylo_tree{};
    tree.ref_sequence = std::move(ref_sequence);
    return tree;
  }
  
  // Pick root time
  auto num_tips = static_cast<int>(std::ssize(tip_descs));
  auto max_root_t = std::ranges::min(tip_descs, {}, [](const auto& tip_desc) { return tip_desc.t; });

  // Set up tree with tips and ensure absence of links
  // Initial mutations are all at tip end
  auto tree = Phylo_tree{2*num_tips - 1};
  tree.ref_sequence = std::move(ref_sequence);
  for (auto tip = 0; tip != std::ssize(tip_descs); ++tip) {
    auto& tip_desc = tip_descs[tip];
    tree.at(tip).parent = k_no_node;
    tree.at(tip).children = {};
    tree.at(tip).name = std::move(tip_desc.name);
    tree.at(tip).t = tip_desc.t;
    for (const auto& m : tip_desc.seq_deltas) {
      tree.at(tip).mutations.push_back(m.place_at(tip_desc.t));
    }
    sort_mutations(tree.at(tip).mutations);
    tree.at(tip).missations = tip_desc.missations;
  }

  // Join up sequentially
  tree.root = 0;
  for (auto i = Node_index{1}; i != num_tips; ++i) {
    // Make node k the parent of tip i and the current root (j), then make k the new root
    auto j = tree.root;
    auto k = i + num_tips - 1;
    tree.at(k).parent = k_no_node;
    tree.at(k).children = {i, j};
    tree.at(i).parent = k;
    tree.at(j).parent = k;
    tree.root = k;
  }

  // Bubble up missations (and remove redundant mutations along the way)
  fix_up_missations(tree);

  // Sanitize inner node times
  for (const auto& node : post_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      auto L = tree.at(node).left_child();
      auto R = tree.at(node).right_child();
      tree.at(node).t = std::min(tree.at(L).t, tree.at(R).t) - 1.0;
    }
  }

  auto scratch = Scratch_space{};
  randomize_mutation_times(tree, bitgen, scratch);

  rereference_to_root_sequence(tree);

  assert_phylo_tree_integrity(tree, true);

  // Randomize tree
  randomize_tree(tree, bitgen);
  
  return tree;
}

auto build_usher_like_tree(
    Real_sequence ref_sequence,
    std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen)
    -> Phylo_tree {
  
  // Sanity-check inputs
  for (auto tip = 0; tip != std::ssize(tip_descs); ++tip) {
    const auto& tip_desc = tip_descs[tip];

    // Mutations
    for (const auto& m : tip_desc.seq_deltas) {
      if (m.site < 0 || m.site >= std::ssize(ref_sequence)) {
        throw std::out_of_range(absl::StrFormat(
            "In tip %d (%s), site of mutation %s is outside valid range [0,%d)",
            tip, tip_desc.name, absl::FormatStreamed(m), std::ssize(ref_sequence)));
      }
      if (auto ref_state = ref_sequence[m.site]; m.from != ref_state) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), mutation %s cannot be applied because ref sequence has state %c there",
            tip, tip_desc.name, absl::FormatStreamed(m), to_char(ref_state)));
      }
      if (m.from == m.to) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), mutation %s has equal 'from' and 'to' states",
            tip, tip_desc.name, absl::FormatStreamed(m)));
      }
    }

    // Missations
    for (const auto& [mi_start, mi_end] : tip_desc.missations.intervals) {
      if (mi_start < 0 || mi_start >= std::ssize(ref_sequence) || mi_end < 0 || mi_end >= (std::ssize(ref_sequence)+1)) {
        throw std::out_of_range(absl::StrFormat(
            "In tip %d (%s), missation interval [%d,%d) applies to a site outside valid range [0,%d)",
            tip, tip_desc.name, mi_start+1, mi_end+1, std::ssize(ref_sequence)));
      }
    }

    // Mutations vs missations
    for (const auto& m : tip_desc.seq_deltas) {
      if (tip_desc.missations.contains(m.site)) {
        throw std::invalid_argument(absl::StrFormat(
            "In tip %d (%s), site %d can't both have a missation and a mutation (%s)",
            tip, tip_desc.name, m.site+1, absl::FormatStreamed(m)));
      }
    }
  }

  if (tip_descs.empty()) {
    auto tree = Phylo_tree{};
    tree.ref_sequence = std::move(ref_sequence);
    return tree;
  }
  
  // Pick root time
  auto num_tips = static_cast<int>(std::ssize(tip_descs));
  auto max_root_t = std::ranges::min(tip_descs, {}, [](const auto& tip_desc) { return tip_desc.t; });

  // Set up tree with tips and ensure absence of links
  // Initial mutations are all at tip end, and missations are all piled up near the tips
  // (a call to fix_up_missations at the end fixes this)
  auto tree = Phylo_tree{2*num_tips - 1};
  tree.ref_sequence = std::move(ref_sequence);
  for (auto tip = 0; tip != std::ssize(tip_descs); ++tip) {
    auto& tip_desc = tip_descs[tip];
    tree.at(tip).name = std::move(tip_desc.name);
    tree.at(tip).t = tip_desc.t;
  }

  // Setting up the first two sequences is tricky (since we want the root sequence to be the ref_sequence)
  CHECK_GE(num_tips, 2);
  {
    auto P = num_tips;  // index of first inner node
    auto A = 0;
    auto B = 1;
    
    tree.root = P;
    tree.at(P).parent = k_no_node;
    
    tree.at(P).children = {A, B};
    tree.at(A).parent = P;
    tree.at(B).parent = P;

    // Pick a reasonable t_P
    // (assume a COVID-like mutation rate of one mutation every 13 days, we get some not-atrocious guess for t_P)
    auto t_A = tree.at(A).t;
    auto t_B = tree.at(B).t;

    auto t_P = std::min(
        t_A - std::ssize(tip_descs[A].seq_deltas) * 13.0,
        t_B - std::ssize(tip_descs[B].seq_deltas) * 13.0)
        - 1.0;

    tree.at(P).t = t_P;

    for (const auto& m : tip_descs[A].seq_deltas) {
      tree.at(A).mutations.push_back(Mutation{
          m.from, m.site, m.to,
          absl::Uniform(absl::IntervalOpenClosed, bitgen, t_P, t_A)});
    }
    sort_mutations(tree.at(A).mutations);
    tree.at(A).missations = tip_descs[A].missations;
    
    for (const auto& m : tip_descs[B].seq_deltas) {
      tree.at(B).mutations.push_back(Mutation{
          m.from, m.site, m.to,
          absl::Uniform(absl::IntervalOpenClosed, bitgen, t_P, t_B)});
    }
    sort_mutations(tree.at(B).mutations);
    tree.at(B).missations = tip_descs[B].missations;
  }
  
  // Sequentially graft every tip where it implies the fewest additional mutations
  auto scratch = Scratch_space{};
  for (auto X = Node_index{2}; X != num_tips; ++X) {
    
    // Study SPR possibilities for attaching X to the budding tree
    auto& tip_desc = tip_descs[X];
    auto t_X = tip_desc.t;
    auto deltas_root_to_X = Site_deltas{scratch};
    for (const auto& delta : tip_desc.seq_deltas) {
      push_back_site_deltas(delta, deltas_root_to_X);
    }
    auto missing_at_X = Scratch_interval_set{scratch};
    missing_at_X = tip_desc.missations.intervals;

    tree.at(X).parent = k_no_node;
    auto builder = Spr_study_builder{tree, k_no_node, tree.at(X).t, missing_at_X, scratch};
    builder.seed_fill_from(tree.root, std::ssize(tree.at_root().mutations), deltas_root_to_X, true);

    // Find least number of mutations required on the P-X branch without changing the rest of the tree
    auto all_min_muts = std::numeric_limits<int>::max();
    for (const auto& region : builder.result) {
      all_min_muts = std::min(all_min_muts, region.min_muts);
    }

    // Pick a random point on one of the regions where grafting implies the minimum number of mutations
    auto chosen_region_idx = -1;
    auto tot_min_T = 0.0;
    for (auto i = 0; i != std::ssize(builder.result); ++i) {
      const auto& region = builder.result[i];
      if (region.min_muts == all_min_muts) {
        if (region.branch == tree.root) {
          // Always pick above the root if that's a possibility
          chosen_region_idx = i;
          break;
        } else {
          tot_min_T += region.t_max - region.t_min;
        }
      }
    }
    if (chosen_region_idx == -1) {  // If inserting above the root isn't optimal...
      auto so_far_min_T = 0.0;
      auto insertion_cum_t = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0.0, tot_min_T);
      for (auto i = 0; i != std::ssize(builder.result); ++i) {
        const auto& region = builder.result[i];
        if (region.min_muts == all_min_muts) {
          so_far_min_T += region.t_max - region.t_min;
          if (insertion_cum_t <= so_far_min_T) {
            chosen_region_idx = i;
            break;
          }
        }
      }
    }
    CHECK_NE(chosen_region_idx, -1);
    
    const auto& chosen_region = builder.result[chosen_region_idx];
    std::cerr << "Attaching tip " << X << " to " << chosen_region << "\n";

    // `chosen_region` is a branch from `G` to `S`.  We'll insert a new node P in the middle of chosen_region,
    // and arrange for P's children to be S and X.
    auto S = chosen_region.branch;
    auto G = tree.at(S).parent;
    auto P = X + num_tips - 1;
    auto t_P = double{};
    auto deltas_P_to_X = Site_deltas{scratch};
    if (S == tree.root) {
      // We're inserting the tip above the root

      // To not change the root sequence, put all the mutations on the P-X branch
      deltas_P_to_X = deltas_root_to_X;
      
      // If we assume a COVID-like mutation rate of one mutation every 13 days, we get some not-atrocious guess for t_P
      auto t_P_guess = t_X - std::ssize(deltas_P_to_X) * 13.0;
      
      // t_P also needs to be above t_S, and we add a bit of slack just in case
      auto t_S = tree.at(S).t;
      t_P = std::min(t_P_guess, t_S) - 1.0;
      
      tree.root = P;
      tree.at(P).mutations.clear();
      swap(tree.at(P).mutations, tree.at(S).mutations);
      
    } else {
      t_P = absl::Uniform(absl::IntervalOpenOpen, bitgen, chosen_region.t_min, chosen_region.t_max);
      
      deltas_P_to_X = deltas_root_to_X;
      displace_site_deltas_start_downwards(tree, deltas_P_to_X, tree.node_loc(tree.root), {S, t_P}, scratch);
      
      auto U = tree.at(G).sibling_of(S);
      tree.at(G).children = {P, U};
      tree.at(P).parent = G;
      CHECK_EQ(tree.at(U).parent, G);

      auto G_S_muts_G_P_end = std::ranges::find_if(tree.at(S).mutations, [&](const auto& m) { return m.t > t_P; });
      auto G_S_muts_end = tree.at(S).mutations.end();

      tree.at(P).mutations.clear();
      swap(tree.at(P).mutations, tree.at(S).mutations);
      tree.at(S).mutations.assign(G_S_muts_G_P_end, G_S_muts_end);
      tree.at(P).mutations.erase(G_S_muts_G_P_end, G_S_muts_end);
    }
    tree.at(P).t = t_P;
    tree.at(P).children = {X, S};
    tree.at(X).parent = P;
    tree.at(S).parent = P;

    // Distribute needed mutations randomly
    for (const auto& [l, delta] : deltas_P_to_X) {
      tree.at(X).mutations.push_back(Mutation{
          delta.from, l, delta.to,
          absl::Uniform(absl::IntervalOpenClosed, bitgen, t_P, t_X)});
    }
    sort_mutations(tree.at(X).mutations);
    tree.at(X).missations = tip_desc.missations;

    // Reuse scratch memory for next tip
    scratch.reset();
  }

  // Bubble up missations (and remove redundant mutations along the way)
  fix_up_missations(tree);

  // Improved guesses for inner node times now that a basic tree is in place
  pseudo_date(tree, bitgen);
  randomize_mutation_times(tree, bitgen, scratch);
  
  assert_phylo_tree_integrity(tree, true);
  
  return tree;
}

}  // namespace delphy
