#ifndef DELPHY_PHYLO_TREE_H_
#define DELPHY_PHYLO_TREE_H_

#include "tree.h"
#include "sequence.h"
#include "mutations.h"

namespace delphy {

struct Phylo_node : public Binary_node {
  std::string name;
  double t;
  Mutation_list<> mutations;
  Missation_map<> missations;

  auto operator==(const Phylo_node& that) const -> bool = default;
};
inline auto operator<<(std::ostream& os, const Phylo_node& node) -> std::ostream& {
  return os << absl::StreamFormat(
      "Phylo_node{name=\"%s\", t=%f, mutations=[%s], missations=%s; parent=%d, children=[%s]}",
      node.name,
      node.t,
      absl::StrJoin(node.mutations, ", ", absl::StreamFormatter()),
      absl::FormatStreamed(node.missations),
      node.parent,
      absl::StrJoin(node.children, ", "));
}

// Phylo tree locations
// ====================

// A Phylo_tree_loc specifies a point on a tree in terms of a branch and a time on that branch.
struct Phylo_tree_loc {
  Branch_index branch;
  double t;
  auto operator<=>(const Phylo_tree_loc& that) const = default;
};
inline auto operator<<(std::ostream& os, const Phylo_tree_loc& loc) -> std::ostream& {
  return os << absl::StreamFormat("[[%d,%g]]", loc.branch, loc.t);
}

// Phylo trees
// ===========

struct Phylo_tree : public Tree<Phylo_node> {
  using Tree::Tree;

  Real_sequence ref_sequence;

  auto num_sites() const -> Site_index { return std::ssize(ref_sequence); }
  auto node_loc(Node_index node) const -> Phylo_tree_loc { return {node, at(node).t}; }
  auto branch_begin_t(Branch_index branch) const -> double { return at_parent_of(branch).t; }
  auto branch_end_t(Branch_index branch) const -> double { return at(branch).t; }
  auto branch_begin(Branch_index branch) const -> Phylo_tree_loc { return {branch, branch_begin_t(branch)}; }
  auto branch_end(Branch_index branch) const -> Phylo_tree_loc { return {branch, branch_end_t(branch)}; }
};

auto operator<<(std::ostream& os, const Phylo_tree& tree) -> std::ostream&;

// Assertions

auto assert_mutation_consistency(const Phylo_tree& tree, bool force = false) -> void;
auto assert_missation_consistency(const Phylo_tree& tree, bool force = false) -> void;
auto assert_phylo_tree_integrity(const Phylo_tree& tree, bool force = false) -> void;

// Common operations
// =================

auto find_MRCA_of(const Phylo_tree& tree, Node_index P, Node_index Q) -> Node_index;
auto find_MRCA_of(const Phylo_tree& tree, Phylo_tree_loc p, Phylo_tree_loc q) -> Phylo_tree_loc;

// Does X descend from A?  Note: X descends from itself, and everyone descends from the "above-root" empty A
auto descends_from(const Phylo_tree& tree, Node_index X, Node_index A) -> bool;
auto descends_from(const Phylo_tree& tree, Phylo_tree_loc x, Phylo_tree_loc a) -> bool;

// Reference sequence changes
auto rereference_to_root_sequence(Phylo_tree& tree) -> void;

// Whole-tree manipulations

// Aggregates all mutations from root to each tip to live in the branch ending at the tips,
// while removing all mutations upstream of all inner nodes.  Does the same for missations.
// This is a useful transformation for analyzing tip sequences and for randomizing trees.
// NOTE: Almost all MCMC and analysis code expects common missations to be pushed as far
// upstream as possible (e.g., assert_missation_consistency will fail).
// Call fix_up_missations() as soon as possible.
auto push_all_mutations_and_missations_to_tips(Phylo_tree& tree) -> void;

// Recursively merges missations that appear on both child branches of a node and pushes
// them up to the parent.  Any mutations on the parent branch for the site with missing data
// are dropped.  Redundant missations in downstream branches are dropped.
// After fix_up_missations, and assuming that `tree` has a valid topology,
// a call to assert_missation_consistency should succeed.
auto fix_up_missations(Phylo_tree& tree) -> void;

// Extract tips and inner nodes
struct extract_nodes_results {
  std::vector<Node_index> tips;
  std::vector<Node_index> inner_nodes;
};
auto extract_nodes(const Phylo_tree& tree) -> extract_nodes_results;

// Rewire a tree in place so that tips in the provided list are connected one at a time:
//  Rewire [x, y] = (x, y)
//  Rewire x :: xs = (x, Rewire(xs))
//
// The tip indices should all be distinct and completely enumerate all tips.
// There should be no mutations on the inner nodes (e.g., the tree has been through push_all_mutations_to_tips).
//
// Afterwards, inner nodes have times equal to min(child.t) - inner_node_offset.
// Mutations are all at the tips and their times are set to the tip times.
auto rewire_tree_through_sequential_accretion(
    Phylo_tree& tree,
    const std::vector<Node_index>& ordered_tips,
    double inner_node_offset = 0.0)
    -> void;

auto randomize_tree(Phylo_tree& tree, absl::BitGenRef bitgen) -> void;

auto randomize_mutation_times(Phylo_tree& tree, absl::BitGenRef bitgen, Scratch_space& scratch) -> void;
auto randomize_branch_mutation_times(
    Phylo_tree& tree,
    Branch_index X,
    absl::BitGenRef bitgen,
    Scratch_space& scratch)
    -> Scratch_vector<Mutation>;

// Whole-tree construction

struct Tip_desc {
  std::string name;
  double t;
  std::vector<Seq_delta> seq_deltas;
  Missation_map<> missations;
};

// Given a reference sequence and a vector of tip descriptions, build a random tree (equivalent to building
// a tree somehow, and then calling randomize_tree).
auto build_random_tree(
    Real_sequence ref_sequence,
    std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen)
    -> Phylo_tree;

// Given a reference sequence and a vector of tip descriptions, build an UShER-like tree (equivalent to greedily
// attaching each sequence in turn to the place in the growing tree where the fewest additional mutations are needed).
auto build_usher_like_tree(
    Real_sequence ref_sequence,
    std::vector<Tip_desc> tip_descs,
    absl::BitGenRef bitgen)
    -> Phylo_tree;

}  // namespace delphy

#endif // DELPHY_PHYLO_TREE_H_
