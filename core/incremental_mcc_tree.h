#ifndef DELPHY_INCREMENTAL_MCC_TREE_H_
#define DELPHY_INCREMENTAL_MCC_TREE_H_

#include "phylo_tree.h"
#include "mcc_tree.h"

namespace delphy {

// TODO: Revive this, improve it so that it detects non-monophyletic mappings too (see incremental_mcc_tree.cpp:175)
// class Incremental_mcc_tree {
//  public:
//   Incremental_mcc_tree(int max_trees, int beam_size, Node_index num_nodes, std::mt19937& prng);

//   auto max_trees() const -> int { return max_trees_; }
//   auto beam_size() const -> int { return beam_size_; }
//   auto num_trees() const -> int { return num_trees_; }

//   auto add_base_tree(Phylo_tree base_tree) -> void;

//   auto derive_mcc_tree() -> Mcc_tree;

//  private:
//   std::mt19937& prng_;

//   int max_trees_;
//   int beam_size_;  // # of trees with top likelihoods which we keep fully up to date

//   // We maintain a circular buffer of at most `max_trees_` base trees, with a new index for every incoming
//   // base tree:
//   // - base_trees_[first_base_tree_index_ % max_trees] is the first tree in the buffer
//   // - base_trees_[(first_base_tree_index_ + num_trees_ - 1) % max_trees] is the last tree in the buffer
//   // If num_trees_ == max_trees_, the buffer is full; to add a new base tree, the first tree must be removed
//   int num_trees_;
//   int first_base_tree_index_;
//   std::vector<Phylo_tree> base_trees_;

//   auto base_tree_by_index(int base_tree_index) -> Phylo_tree& {
//     assert(base_tree_index >= first_base_tree_index_);
//     assert(base_tree_index < first_base_tree_index_ + num_trees_);
//     return base_trees_[base_tree_index % max_trees_];
//   }
//   auto pop_first_base_tree() -> void {
//     assert(num_trees_ > 0);
//     base_tree_by_index(first_base_tree_index_) = {};  // Clear out buffers
//     ++first_base_tree_index_;
//     --num_trees_;
//   }
//   auto get_next_base_tree_index() const -> int { return first_base_tree_index_ + num_trees_; }
//   auto push_base_tree(Phylo_tree new_base_tree) -> int {
//     assert(num_trees_ < max_trees_);
//     int new_base_tree_index = first_base_tree_index_ + num_trees_;
//     base_trees_[new_base_tree_index % max_trees_] = std::move(new_base_tree);
//     ++num_trees_;
//     return new_base_tree_index;
//   }

//   std::vector<uint64_t> tip_fingerprints_;  // Random fingerprint of i'th leaf (corresponding to (N_s - 1 + i)'th node)

//   // clade_base_nodes_[clade] is a map of base tree index to node index for each appearance of `clade`
//   std::unordered_map<uint64_t, std::unordered_map<int, Node_index>> clade_base_nodes_;

//   struct Beam_tree {
//     int base_tree_index;
//     double log_clade_credibility;
//   };

//   std::vector<Beam_tree> beam_trees_;

//   void adjust_clade_base_nodes(const Phylo_tree& base_tree, int base_tree_index, bool add_node);
//   double calc_log_clade_credibility(const Phylo_tree& base_tree);
// };

}  // namespace delphy

#endif // DELPHY_INCREMENTAL_MCC_TREE_H_
