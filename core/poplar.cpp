#include "poplar.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/strings/str_format.h"

#include "phylo_tree_calc.h"

namespace delphy {

auto map_site_tree_branch(
    const Phylo_tree& tree,
    Branch_index full_branch,
    Site_index l)
    -> std::pair<Node_index, Node_index> {

  // Optional check that full branch exists in pruned site tree
  constexpr auto full_missing_site_check = false;
  if (full_missing_site_check && is_site_missing_at(tree, full_branch, l)) {
    std::cerr << absl::StreamFormat(
        "POSSIBLE BUG! Call of map_site_tree_branch(tree, full_branch=%d, l=%d) when site is missing there\n",
        full_branch, l);
    return {k_no_node, k_no_node};
  }
  
  // Start with the trivial path of the branch in the full tree
  auto end_node = Node_index{full_branch};
  auto start_node = tree.at(end_node).parent;

  // The branch above the root doesn't really exist in the full tree or in any pruned site tree
  if (end_node == tree.root) {
    return {k_no_node, k_no_node};
  }

  // If the `full_branch` contains a missation at l, then `full_branch` isn't even in the
  // pruned site tree
  if (tree.at(end_node).missations.contains(l)) {
    return {k_no_node, k_no_node};
  }
  
  // Push the end node as far downstream as it can go, as long as there are always
  // missations on the offshoots
  while (not tree.at(end_node).is_tip()) {
    auto next_node = k_no_node;
    auto num_children_including_site = 0;
    for (const auto& child : tree.at(end_node).children) {
      if (not tree.at(child).missations.contains(l)) {
        ++num_children_including_site;
        next_node = child;
      }
    }

    if (num_children_including_site > 1) { break; }  // Found one endpoint of the branch in the pruned site tree

    // Missation factoring guarantees that a node that isn't missing a site has at least
    // one child that isn't missing it, either.  Crash rather than index the tree with
    // `k_no_node` if that is somehow violated
    CHECK_NE(next_node, k_no_node);
    end_node = next_node;
  }

  // Push the start node as far upstream as it can go, as long as there are always
  // missations on the offshoots
  auto prev_node = Node_index{full_branch};
  while (start_node != tree.root) {
    auto sibling = tree.at(start_node).sibling_of(prev_node);
    if (tree.at(sibling).missations.contains(l)) {
      // If `sibling` has a missation at `l`, then by the consistency condition on
      // missations, there can't be a missation for `l` upstream of `sibling`, too
      CHECK(not tree.at(prev_node).missations.contains(l));

      prev_node = start_node;
      start_node = tree.at(start_node).parent;
    } else {
      break;
    }
  }

  // If we reach the root, funny things start happening...
  if (start_node == tree.root) {
    // ...the site might be missing from the entire tree
    if (tree.at_root().missations.contains(l)) {
      return {k_no_node, k_no_node};
    }

    // ...the site may be missing from the other side of the root, in which case the
    // branch only exists on the site tree because we insist on the site trees all having
    // the same root
    auto sibling = tree.at_root().sibling_of(prev_node);
    if (tree.at(sibling).missations.contains(l)) {
      return {k_no_node, end_node};
    }
    
    // ...or we may have to go round the root node and push *downstream*, just as with the end_node
    start_node = sibling;

    // Push the start node as far downstream as it can go, as long as there are always
    // missations on the offshoots
    while (not tree.at(start_node).is_tip()) {
      auto next_node = k_no_node;
      auto num_children_including_site = 0;
      for (const auto& child : tree.at(start_node).children) {
        if (not tree.at(child).missations.contains(l)) {
          ++num_children_including_site;
          next_node = child;
        }
      }

      if (num_children_including_site > 1) { break; }  // Found the other endpoint of the branch in the pruned site tree

      CHECK_NE(next_node, k_no_node);  // See above
      start_node = next_node;
    }
  }

  return {start_node, end_node};  // Common case
}

auto map_mutation_extent(
    const Phylo_tree& tree,
    Branch_index branch,
    Site_index l)
    -> std::optional<Mutation_extent> {

  auto [start_node, end_node] = map_site_tree_branch(tree, branch, l);
  
  // Standardize to start_node <= end_node, even if that means we lose info on which is
  // the MRCA(start, end) whenever the start<->end path doesn't include the root
  if (start_node > end_node) {
    std::swap(start_node, end_node);
  }

  // Defend against pruned site tree branches the don't exist (see above)
  static_assert(k_no_node < 0);  // So if there's only one k_no_node, it's the start node
  if (start_node == k_no_node && end_node == k_no_node) {
    return {};
  }

  if (start_node == k_no_node) {
    // The dangling branch from the root to the first node with informative tips
    // downstream (edge case (c)).  Handle it specially, since this "extent" isn't one at
    // all: it only forces a state on the root blob (see `map_tree_mutation_extents`) and
    // never decides where the blob boundaries are
    CHECK_NE(end_node, k_no_node);

    auto maybe_state_at_end = std::optional<Real_seq_letter>{};
    auto maybe_state_at_root = std::optional<Real_seq_letter>{};
    for (auto cur = end_node; cur != tree.root; cur = tree.at(cur).parent) {
      for (const auto& m : tree.at(cur).mutations | std::views::reverse) {
        if (m.site == l) {
          if (not maybe_state_at_end.has_value()) {
            maybe_state_at_end = m.to;
          }
          maybe_state_at_root = m.from;
        }
      }
    }

    if (not maybe_state_at_root.has_value()) {
      // No mutation here, no need to report an extent
      return {};
    }
    
    CHECK(maybe_state_at_end.has_value());
    CHECK(maybe_state_at_root.has_value());

    return Mutation_extent{l, k_no_node, *maybe_state_at_root, end_node, *maybe_state_at_end, k_no_node};
  }

  // We could organize the whole blob calculation more efficiently to avoid recalculating
  // the MRCA below and to avoid mapping the same extent multiple times in case of
  // multiple mutations on the same site tree branch.  We explicitly sacrifice that
  // efficiency for obvious correctness here

  auto mrca_node = find_MRCA_of(tree, start_node, end_node);

  // Walk upwards from each endpoint to the MRCA finding the initial and final states of
  // those walks (if there is any mutation on them).
  auto maybe_state_at_mrca_node_from_end = std::optional<Real_seq_letter>{};
  auto maybe_state_at_end_node = std::optional<Real_seq_letter>{};
  for (auto cur = end_node; cur != mrca_node; cur = tree.at(cur).parent) {
    for (const auto& m : tree.at(cur).mutations | std::views::reverse) {
      if (m.site == l) {
        if (not maybe_state_at_end_node.has_value()) {
          maybe_state_at_end_node = m.to;
        }
        maybe_state_at_mrca_node_from_end = m.from;
      }
    }
  }
  
  auto maybe_state_at_mrca_node_from_start = std::optional<Real_seq_letter>{};
  auto maybe_state_at_start_node = std::optional<Real_seq_letter>{};
  for (auto cur = start_node; cur != mrca_node; cur = tree.at(cur).parent) {
    for (const auto& m : tree.at(cur).mutations | std::views::reverse) {
      if (m.site == l) {
        if (not maybe_state_at_start_node.has_value()) {
          maybe_state_at_start_node = m.to;
        }
        maybe_state_at_mrca_node_from_start = m.from;
      }
    }
  }

  // Now reconcile all the states
  if (maybe_state_at_mrca_node_from_start.has_value() && maybe_state_at_mrca_node_from_end.has_value()) {
    CHECK_EQ(*maybe_state_at_mrca_node_from_start, *maybe_state_at_mrca_node_from_end);
  }
  if (not maybe_state_at_mrca_node_from_start.has_value() && not maybe_state_at_mrca_node_from_end.has_value()) {
    // No mutation on this site tree branch
    return {};
  }
  auto state_at_mrca_node =
      maybe_state_at_mrca_node_from_start
      ? maybe_state_at_mrca_node_from_start.value()
      : maybe_state_at_mrca_node_from_end.value();
  auto state_at_start_node = maybe_state_at_start_node.value_or(state_at_mrca_node);
  auto state_at_end_node = maybe_state_at_end_node.value_or(state_at_mrca_node);

  // Edge case (b): 2+ mutations with no net effect.  We still report these as "neutral"
  // extents (equal states at both end-points) and leave it to the caller to decide what
  // to do with them.

  return Mutation_extent{l, start_node, state_at_start_node, end_node, state_at_end_node, mrca_node};
}

auto map_tree_mutation_extents(const Phylo_tree& tree) -> Mapped_mutation_extents {
  auto result = Mapped_mutation_extents{};

  // Do this in the most straightforward possible way.  As in map_mutation_extent, we
  // could organize the calculation more efficiently, but prefer obvious correctness
  // instead.
  for (const auto& node : index_order_traversal(tree)) {
    if (node == tree.root) { continue; }
    for (const auto& m : tree.at(node).mutations) {
      if (auto maybe_e = map_mutation_extent(tree, node, m.site); maybe_e) {
        const auto& e = *maybe_e;
        
        if (e.start_node == k_no_node) {
          // "Half" extent (edge case (c)): the root blob has the state at the other end
          // of the extent
          result.root_state_overrides.push_back({e.site, e.state_at_end});
          
        } else if (e.state_at_start == e.state_at_end) {
          // "Neutral" extent (edge case (b)): nothing to record, since every blob along
          // it has the common state anyway, unless it goes through the root (edge case (d)),
          // in which case the root blob has that common state, instead of the root's
          // actual state
          if (e.mrca_node == tree.root) {
            result.root_state_overrides.push_back({e.site, e.state_at_start});
          }
          
        } else {
          result.extents.push_back(e);
        }
      }
    }
  }

  // Sort and deduplicate
  std::ranges::sort(result.extents);
  result.extents.erase(std::unique(result.extents.begin(), result.extents.end()), result.extents.end());
  std::ranges::sort(result.root_state_overrides);
  result.root_state_overrides.erase(
      std::unique(result.root_state_overrides.begin(), result.root_state_overrides.end()),
      result.root_state_overrides.end());

  return result;
}

auto annotate_branches_with_extents(
    const Phylo_tree& tree,
    const std::vector<Mutation_extent>& extents)
    -> Branch_vector<Branch_extent_info> {

  // Empty markers for every branch
  auto result = Branch_vector<Branch_extent_info>(std::ssize(tree), Branch_extent_info{{}, {}});

  // Add extent id to the branches it touches
  for (auto ie = Mutation_extent_index{0}; ie != std::ssize(extents); ++ie) {
    const auto& e = extents[ie];

    if (e.mrca_node == tree.root) {
      // If the mrca is the root, then both extent end-points are at the end of their respective branches
      result[e.start_node].mutation_extents_at_branch_end.push_back(ie);
      result[e.end_node].mutation_extents_at_branch_end.push_back(ie);
    } else {
      // Otherwise, the mrca end lies at the start of the branch, and the other lies at the end
      // We have to walk up from the non-mrca end to find the correct child of the mrca node
      auto later_node = (e.mrca_node == e.start_node) ? e.end_node : e.start_node;
      auto prev = later_node;
      auto cur = later_node;
      while (cur != e.mrca_node) {
        prev = cur;
        cur = tree.at(cur).parent;
      }

      result[prev].mutation_extents_at_branch_start.push_back(ie);
      result[later_node].mutation_extents_at_branch_end.push_back(ie);
    }
  }

  // Sort results for to make output deterministic
  for (auto& info : result) {
    std::ranges::sort(info.mutation_extents_at_branch_start);
    std::ranges::sort(info.mutation_extents_at_branch_end);
  }
  
  return result;
}

// The ambiguous state at a site covered by a (real, i.e., non-neutral) mutation extent
static auto ambiguous_state_of(const Mutation_extent& e) -> Seq_letter {
  return to_seq_letter(e.state_at_start) | to_seq_letter(e.state_at_end);
}

auto derive_blob_tree(const Phylo_tree& tree) -> Blob_tree {

  auto blob_tree = Blob_tree{};

  // Map mutation extents and annotate branches where they end
  auto mapped_extents = map_tree_mutation_extents(tree);
  blob_tree.extents = std::move(mapped_extents.extents);
  blob_tree.root_state_overrides = std::move(mapped_extents.root_state_overrides);
  blob_tree.branch_extents = annotate_branches_with_extents(tree, blob_tree.extents);
  
  // Initially, nothing is matched to a blob
  blob_tree.node_2_blob.assign(std::ssize(tree), k_no_blob);
  blob_tree.branch_2_blob.assign(std::ssize(tree), k_no_blob);

  // Kickstart traversal with a stump of a root blob (and record the extents active at the root)
  auto root_blob = blob_tree.add_node();
  blob_tree.root = root_blob;
  blob_tree.at(root_blob).parent = k_no_blob;
  blob_tree.node_2_blob[tree.root] = root_blob;
  blob_tree.at(root_blob).seed_node = tree.root;
  blob_tree.at(root_blob).seed_branch = k_no_branch;
  blob_tree.at(root_blob).nodes.push_back(tree.root);
  for (auto ie = Mutation_extent_index{0}; ie != std::ssize(blob_tree.extents); ++ie) {
    const auto& e = blob_tree.extents[ie];
    if (e.mrca_node == tree.root) {
      blob_tree.at(root_blob).active_mutation_extents.push_back(ie);
    }
  }

  // Calculate the root blob's sequence deltas & fingerprint from scratch, without materializing its sequence.
  // The root blob's sequence is the reference sequence, except at the sites missing in the whole tree (kept
  // aside, see `Blob_tree::missing_sites`, so they never show up in any sequence delta) and at the sites
  // where the root's actual sequence differs from the reference sequence, where there is a
  // `Root_state_override`, or that are covered by an extent going through the root (in increasing order of
  // precedence).
  blob_tree.missing_sites = tree.at_root().missations.intervals;
  auto root_states = Site_map<Seq_letter>{};
  for (const auto& m : tree.at_root().mutations) {
    root_states[m.site] = to_seq_letter(m.to);
  }
  for (const auto& o : blob_tree.root_state_overrides) {
    root_states[o.site] = to_seq_letter(o.state);
  }
  for (const auto& ie : blob_tree.at(root_blob).active_mutation_extents) {
    const auto& e = blob_tree.extents[ie];
    root_states[e.site] = ambiguous_state_of(e);
  }

  // The root blob's fingerprint starts out as the reference sequence's, with the sites missing in the whole tree as N...
  auto root_blob_fp = k_empty_blob_fingerprint;
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    root_blob_fp ^= site_state_fingerprint(l, tree.ref_sequence[l]);
  }
  for (const auto& [start, end] : blob_tree.missing_sites) {
    for (auto l = Site_index{start}; l != end; ++l) {
      root_blob_fp ^=
          site_state_fingerprint(l, tree.ref_sequence[l]) ^
          site_state_fingerprint(l, Seq_letters::N);
    }
  }

  // ...and then changes with every site where the root blob differs from the reference sequence
  for (const auto& [site, state] : root_states) {
    CHECK(not blob_tree.missing_sites.contains(site)) << "Mutation below a missation at site " << site;
    auto from_state = to_seq_letter(tree.ref_sequence[site]);
    if (state != from_state) {
      root_blob_fp ^= site_state_fingerprint(site, from_state) ^ site_state_fingerprint(site, state);
      blob_tree.at(root_blob).sequence_delta.push_back({site, from_state, state});
    }
  }
  blob_tree.at(root_blob).fingerprint = root_blob_fp;
  
  // In the following traversal, the invariant is that the parent node of the node being
  // visited has been assigned to a blob.  We start at the root
  auto cur_blob = root_blob;
  auto cur_fp = root_blob_fp;
  auto cur_es = absl::flat_hash_set<Mutation_extent_index>{};
  cur_es.insert(blob_tree.at(root_blob).active_mutation_extents.begin(),
                blob_tree.at(root_blob).active_mutation_extents.end());
  
  // Record a change in the state of a site in `blob` relative to its parent blob, keeping `cur_fp` in sync
  auto record_delta = [&](Blob_index blob, Site_index site, Seq_letter from_state, Seq_letter to_state) {
    cur_fp ^= site_state_fingerprint(site, from_state) ^ site_state_fingerprint(site, to_state);
    blob_tree.at(blob).sequence_delta.push_back({site, from_state, to_state});
  };

  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (node == tree.root) { continue; }  // We dealt with the root separately above
    
    if (children_so_far == 0) {
      // Entering `node`
      // At this point, cur_* reflects the state at `parent`
      auto parent = tree.at(node).parent;
      auto branch = node;  // parent->node branch
      
      // If there are extents at the beginning of the parent->node branch, that kicks off a new blob
      if (not blob_tree.branch_extents[branch].mutation_extents_at_branch_start.empty()) {
        
        auto new_blob = blob_tree.add_node();
        blob_tree.at(cur_blob).children.push_back(new_blob);
        blob_tree.at(new_blob).parent = cur_blob;
        cur_blob = new_blob;

        blob_tree.at(cur_blob).seed_node = k_no_node;
        blob_tree.at(cur_blob).seed_branch = branch;
        
        for (const auto& ie : blob_tree.branch_extents[branch].mutation_extents_at_branch_start) {
          const auto& e = blob_tree.extents[ie];
          
          cur_es.insert(ie);
          
          // Remove the "certain" state at this site...
          CHECK(parent == e.start_node || parent == e.end_node) << parent << "-" << node << "; " << e;
          auto old_certain_state = (e.start_node == parent) ? e.state_at_start : e.state_at_end;
          
          // ...by the "uncertain" state encoded in the mutation extent
          record_delta(cur_blob, e.site, to_seq_letter(old_certain_state), ambiguous_state_of(e));
        }
        
        blob_tree.at(cur_blob).fingerprint = cur_fp;
        blob_tree.at(cur_blob).active_mutation_extents.assign(cur_es.begin(), cur_es.end());
      }

      // cur_blob is now accurate for the branch, so add the branch to the current blob
      blob_tree.at(cur_blob).branches.push_back(branch);
      blob_tree.branch_2_blob[branch] = cur_blob;

      // If there are extents at the end of the parent->node branch, that also kicks off a new blob
      if (not blob_tree.branch_extents[branch].mutation_extents_at_branch_end.empty()) {
        
        auto new_blob = blob_tree.add_node();
        blob_tree.at(cur_blob).children.push_back(new_blob);
        blob_tree.at(new_blob).parent = cur_blob;
        cur_blob = new_blob;
        
        blob_tree.at(cur_blob).seed_node = node;
        blob_tree.at(cur_blob).seed_branch = k_no_branch;
        
        for (const auto& ie : blob_tree.branch_extents[branch].mutation_extents_at_branch_end) {
          const auto& e = blob_tree.extents[ie];
          
          cur_es.erase(ie);
          
          // Remove the "uncertain" state encoded in the mutation extent...
          // ...by the "certain" state past the endpoint
          CHECK(node == e.start_node || node == e.end_node) << parent << "-" << node << "; " << e;
          auto new_certain_state = (e.start_node == node) ? e.state_at_start : e.state_at_end;
          record_delta(cur_blob, e.site, ambiguous_state_of(e), to_seq_letter(new_certain_state));
        }
        
        blob_tree.at(cur_blob).fingerprint = cur_fp;
        blob_tree.at(cur_blob).active_mutation_extents.assign(cur_es.begin(), cur_es.end());
      }

      // cur_blob is now accurate for the node, so add the node to the current blob
      blob_tree.at(cur_blob).nodes.push_back(node);
      blob_tree.node_2_blob[node] = cur_blob;
      
      // At this point, cur_* reflect the state at `node`
    }

    if (children_so_far == std::ssize(tree.at(node).children)) {
      // Exiting `node`
      // At this point, cur_* reflect the state at `node`
      auto parent = tree.at(node).parent;
      //auto branch = node;  // parent->node branch

      // Don't bother first moving back to `branch` and then to `parent`.  Instead,
      // move back to `parent` directly (can do this directly since we rebuild cur_es from scratch)
      DCHECK_NE(blob_tree.node_2_blob[parent], k_no_blob);
      if (cur_blob != blob_tree.node_2_blob[parent]) {
        cur_blob = blob_tree.node_2_blob[parent];
        cur_fp = blob_tree.at(cur_blob).fingerprint;
        cur_es.clear();
        cur_es.insert(blob_tree.at(cur_blob).active_mutation_extents.begin(),
                      blob_tree.at(cur_blob).active_mutation_extents.end());
      }
      
      // At this point, cur_* reflect the state at `parent`
    }
  }

  // Sort all vectors for determinism
  for (const auto& blob : index_order_traversal(blob_tree)) {
    std::ranges::sort(blob_tree.at(blob).nodes);
    std::ranges::sort(blob_tree.at(blob).branches);
    std::ranges::sort(blob_tree.at(blob).active_mutation_extents);
    std::ranges::sort(blob_tree.at(blob).sequence_delta, {}, &Blob_site_delta::site);
  }

  return blob_tree;
}

// Split a state that is ambiguous between exactly 2 bases (e.g., `{A,G}`) into those 2 bases, in A < C < G < T order
static auto make_blob_sequence_ambiguity(Site_index site, Seq_letter ambiguous_state) -> Blob_sequence_ambiguity {
  auto bases = std::array<Real_seq_letter, 2>{};
  auto num_bases = 0;
  for (const auto& a : {Real_seq_letter::A, Real_seq_letter::C, Real_seq_letter::G, Real_seq_letter::T}) {
    if (to_seq_letter(a) & ambiguous_state) {
      CHECK_LT(num_bases, 2) << "Site " << site;
      bases[num_bases++] = a;
    }
  }
  CHECK_EQ(num_bases, 2) << "Site " << site;
  return {site, bases[0], bases[1]};
}

auto query_blob_sequence(
    const Blob_tree& blob_tree,
    Blob_index blob)
    -> Blob_sequence {

  auto result = Blob_sequence{};

  // The sites missing in the whole tree are the same for every blob
  result.missing_sites = blob_tree.missing_sites;

  // Walk up from `blob` to the root blob, combining the sequence deltas along the way.  For every site touched by
  // any of them, `states` tracks the state at the bottom of the walk (the to_state of the delta closest to `blob`,
  // i.e., the state at `blob`) and the state at the top (the from_state of the delta closest to the root blob,
  // which is the reference sequence's state, since the root blob's deltas are relative to it).
  struct Site_states {
    Seq_letter topmost_from_state;
    Seq_letter bottommost_to_state;
  };
  auto states = Site_map<Site_states>{};
  for (auto b = blob; b != k_no_blob; b = blob_tree.at(b).parent) {
    for (const auto& [site, from_state, to_state] : blob_tree.at(b).sequence_delta) {
      auto [it, inserted] = states.try_emplace(site, Site_states{from_state, to_state});
      if (not inserted) {
        // This delta is further up than the ones we've seen so far at this site, so it has to
        // end where the delta right below it started
        CHECK_EQ(to_state, it->second.topmost_from_state) << "Blob " << b << ", site " << site;
        it->second.topmost_from_state = from_state;
      }
    }
  }

  // Flatten, leaving out sites that ended up back where they started
  for (const auto& [site, site_states] : states) {
    auto [topmost_from_state, bottommost_to_state] = site_states;

    // The topmost delta at a site starts from the reference sequence's state, which is always a certain base
    CHECK(not is_ambiguous(topmost_from_state)) << "Site " << site;

    if (bottommost_to_state == topmost_from_state) {
      continue;
    } else if (is_ambiguous(bottommost_to_state)) {
      result.ambiguities.push_back(make_blob_sequence_ambiguity(site, bottommost_to_state));
    } else {
      result.mutations.push_back(
          {site, to_real_seq_letter(topmost_from_state), to_real_seq_letter(bottommost_to_state)});
    }
  }
  std::ranges::sort(result.mutations, {}, &Blob_sequence_mutation::site);
  std::ranges::sort(result.ambiguities, {}, &Blob_sequence_ambiguity::site);

  return result;
}

auto calc_blob_fingerprint(
    const Real_sequence& ref_sequence,
    const Blob_sequence& blob_sequence)
    -> Blob_fingerprint {

  auto result = k_empty_blob_fingerprint;

  // First, fingerprint the reference sequence
  for (auto l = Site_index{0}; l != std::ssize(ref_sequence); ++l) {
    result ^= site_state_fingerprint(l, ref_sequence[l]);
  }

  // Flag all the missing sites
  for (const auto& [start, end] : blob_sequence.missing_sites) {
    for (auto l = Site_index{start}; l != end; ++l) {
      result ^=
          site_state_fingerprint(l, ref_sequence[l]) ^
          site_state_fingerprint(l, Seq_letters::N);
    }
  }

  // Flip all definitely mutated sites
  for (const auto& [l, from_state, to_state] : blob_sequence.mutations) {
    DCHECK(not blob_sequence.missing_sites.contains(l));
    CHECK_EQ(from_state, ref_sequence[l]);
    result ^=
        site_state_fingerprint(l, ref_sequence[l]) ^
        site_state_fingerprint(l, to_state);
  }

  // Mark all ambiguous sites
  for (const auto& [l, state1, state2] : blob_sequence.ambiguities) {
    DCHECK(not blob_sequence.missing_sites.contains(l));
    result ^=
        site_state_fingerprint(l, ref_sequence[l]) ^
        site_state_fingerprint(l, state1) ^
        site_state_fingerprint(l, state2);
  }
  
  return result;
}

}  // namespace delphy
