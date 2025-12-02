#ifndef DELPHY_MCC_TREE_H_
#define DELPHY_MCC_TREE_H_

#include <memory>
#include <vector>
#include <random>

#include <absl/random/bit_gen_ref.h>

#include "tree.h"
#include "phylo_tree.h"

namespace delphy {

class Mcc_tree;

// An Mcc_tree summarizes the properties of a list of related "base trees"
using Base_tree_index = int;

// Stores a `T` for every base tree in an Mcc_tree, indexed by a base tree's index
template<typename T, typename Alloc = std::allocator<T>>
using Base_tree_vector = std::vector<T, Alloc>;

// Each node in the MCC defines a useful grouping of tips.  The "corresponding node" in each base tree is the MRCA
// of those tips in the base tree.  If those tips are the exact set of descendant tips of the corresponding node,
// we say that the corresponding node is monophyletic in the base tree.  More concisely, we say that the corresponding
// node is an "exact match" to the MCC node in that base tree.
// The properties of the usual MCC trees reflect only the properties the corresponding nodes that are exact matches.
struct Mcc_corresponding_node_info {  // Meant to be the value type in a Base_tree_vector or a Base_tree_map
  Node_index node_in_base_tree;
  bool is_monophyletic_in_base_tree;
};

class Mcc_node : public Binary_node {
 public:
  Mcc_node() : corresponding_node_infos_{} {}
  explicit Mcc_node(Base_tree_vector<Mcc_corresponding_node_info> corresponding_node_infos)
      : corresponding_node_infos_{std::move(corresponding_node_infos)} {}

  auto corresponding_node_infos() -> Base_tree_vector<Mcc_corresponding_node_info>& {
    return corresponding_node_infos_; }
  auto corresponding_node_infos() const -> const Base_tree_vector<Mcc_corresponding_node_info>& {
    return corresponding_node_infos_; }

  // Posterior support = fraction of base trees with an exact match for this node
  auto posterior_support() const -> double { return posterior_support_; }
  auto set_posterior_support(double p) -> void { posterior_support_ = p; }

  // Mean of base tree node times over exact matches only
  auto t() const -> double { return t_; }
  auto set_t(double t) -> void { t_ = t; }

  // Mean of base tree node times over *all* corresponding nodes
  auto t_mrca() const -> double { return t_mrca_; }
  auto set_t_mrca(double t_mrca) -> void { t_mrca_ = t_mrca; }
  
 private:
  Base_tree_vector<Mcc_corresponding_node_info> corresponding_node_infos_;
  double t_{};
  double t_mrca_{};
  double posterior_support_{};
};

class Mcc_tree : public Tree<Mcc_node> {
 public:
  Mcc_tree() {}
  Mcc_tree(Node_index num_nodes, Base_tree_vector<Phylo_tree*> base_trees, int master_base_tree_index)
      : Tree{num_nodes},
        base_trees_{std::move(base_trees)},
        master_base_tree_index_{master_base_tree_index} {}

  auto num_base_trees() const -> Base_tree_index { return std::ssize(base_trees_); }
  auto base_trees() -> std::vector<Phylo_tree*>& { return base_trees_; }
  auto base_trees() const -> const std::vector<Phylo_tree*>& { return base_trees_; }
  auto master_base_tree_index() const -> Base_tree_index { return master_base_tree_index_; }

  auto calculate_derived_quantities() -> void;

 private:
  Base_tree_vector<Phylo_tree*> base_trees_;
  int master_base_tree_index_ = -1;
};

inline
auto corresponding_node_to(
    const Mcc_tree& mcc_tree,
    Node_index mcc_node,
    Base_tree_index base_tree_index)
    -> Node_index {
  return mcc_tree.at(mcc_node).corresponding_node_infos().at(base_tree_index).node_in_base_tree;
}

// The main show
// =============
auto derive_mcc_tree(Base_tree_vector<Phylo_tree*> base_trees, absl::BitGenRef bitgen) -> Mcc_tree;

}  // namespace delphy

#endif // DELPHY_MCC_TREE_H_
