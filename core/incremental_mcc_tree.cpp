#include <random>

#include "incremental_mcc_tree.h"

// TODO: Refactor to use Mcc_tree for tracking of base trees, clade frequencies and clade -> base node mappings,
//   then just use beam search to pick the topology of the mcc tree.  (maybe factor the common stuff out
//   to something like a Clade_repository).

// TODO: Switch value clade_base_nodes_ to a btree_map with a compound (clade, base_node_tree_index) key to
//   easily walk

// namespace delphy {

// Incremental_mcc_tree::Incremental_mcc_tree(int max_trees, int beam_size, Node_index num_nodes, std::mt19937& prng)
//     : prng_{prng},
//       max_trees_{max_trees},
//       beam_size_{beam_size},
//       num_trees_{0},
//       first_base_tree_index_{0},
//       base_trees_(max_trees, Phylo_tree{}) {
//   for (auto i = 0; i != num_nodes; ++i) {
//     tip_fingerprints_.push_back(std::uniform_int_distribution<uint64_t>{}(prng_));
//   }
// }

// auto Incremental_mcc_tree::adjust_clade_base_nodes(
//     const Phylo_tree& base_tree,
//     int base_tree_index,
//     bool add_node)
//     -> void {
//   auto node_clade = Node_vector<uint64_t>(base_tree.size());
//   for (const auto& node : post_order_traversal(base_tree)) {
//     auto& clade = node_clade[node.index()];
//     if (node.is_tip()) {
//       clade = tip_fingerprints_[node.index()];
//     } else {
//       clade = node_clade[node.left_child().index()] ^ node_clade[node.right_child().index()];
//     }

//     if (add_node) {
//       clade_base_nodes_[clade].emplace(base_tree_index, node.index());
//     } else {
//       auto& this_clade_base_nodes = clade_base_nodes_[clade];
//       this_clade_base_nodes.erase(base_tree_index);
//       if (this_clade_base_nodes.empty()) {
//         clade_base_nodes_.erase(clade);
//       }
//     }
//   }
// }

// auto Incremental_mcc_tree::add_base_tree(Phylo_tree base_tree) -> void {
//   // If we're already storing `max_trees_` base trees, slide out the first one
//   if (num_trees_ == max_trees_) {
//     auto base_tree_index = first_base_tree_index_;
//     const auto& first_tree = base_tree_by_index(base_tree_index);
//     adjust_clade_base_nodes(first_tree, base_tree_index, false);
//     pop_first_base_tree();

//     auto beam_tree_it_to_remove = std::ranges::find(beam_trees_, base_tree_index,
//                                                     [](const auto& beam_tree) { return beam_tree.base_tree_index; });
//     if (beam_tree_it_to_remove != beam_trees_.end()) {
//       *beam_tree_it_to_remove = beam_trees_.back();
//       beam_trees_.pop_back();
//     }
//   }

//   // Add in new tree
//   auto new_base_tree_index = get_next_base_tree_index();
//   adjust_clade_base_nodes(base_tree, new_base_tree_index, true);
//   push_base_tree(std::move(base_tree));

//   // Recalculate log clade credibility of each beam tree
//   for (auto& beam_tree : beam_trees_) {
//     beam_tree.log_clade_credibility = calc_log_clade_credibility(base_tree_by_index(beam_tree.base_tree_index));
//   }

//   // Pick a random other tree, calculate its log clade credibility and see if it should enter the beam
//   // (possibly displacing another beam tree)
//   if (num_trees_ > std::ssize(beam_trees_)) {
//     auto random_base_tree_index = std::uniform_int_distribution{first_base_tree_index_, first_base_tree_index_ + num_trees_ - 1}(prng_);

//     if (std::ranges::find(beam_trees_, random_base_tree_index,
//                           [](const auto& beam_tree) { return beam_tree.base_tree_index; }) == beam_trees_.end()) {
//       // Random tree is not already in beam
//       auto log_clade_credibility = calc_log_clade_credibility(base_tree_by_index(random_base_tree_index));

//       if (std::ssize(beam_trees_) < beam_size_) {
//         // Add unconditionally
//         beam_trees_.emplace_back(Beam_tree{random_base_tree_index, log_clade_credibility});
//       } else {
//         auto weakest_beam_tree =
//             std::ranges::min_element(beam_trees_.begin(), beam_trees_.end(), {},
//                                      [](const auto& beam_tree) { return beam_tree.log_clade_credibility; });
//         if (log_clade_credibility > weakest_beam_tree->log_clade_credibility) {
//           // Replace
//           *weakest_beam_tree = Beam_tree{random_base_tree_index, log_clade_credibility};
//         }
//       }
//     }
//   }
// }

// auto Incremental_mcc_tree::calc_log_clade_credibility(const Phylo_tree& base_tree) -> double {
//   auto result = 0.0;

//   auto node_clade = Node_vector<uint64_t>(base_tree.size());
//   for (const auto& node : post_order_traversal(base_tree)) {
//     auto& clade = node_clade[node.index()];
//     if (node.is_tip()) {
//       clade = tip_fingerprints_[node.index()];
//     } else {
//       clade = node_clade[node.left_child().index()] ^ node_clade[node.right_child().index()];
//     }

//     result += std::log(static_cast<double>(std::ssize(clade_base_nodes_[clade])) / num_trees_);
//   }

//   return result;
// }

// auto Incremental_mcc_tree::derive_mcc_tree() -> Mcc_tree {
//   if (beam_trees_.empty()) {
//     return Mcc_tree{};
//   }

//   // Pick member of beam with highest clade credibility, and use it to make the Mcc_tree
//   const auto& best_beam_tree =
//       std::ranges::max(beam_trees_,
//                        {}, [](const auto& beam_tree) { return beam_tree.log_clade_credibility; });
//   auto master_base_tree_index = best_beam_tree.base_tree_index;

//   auto base_trees = std::vector<Phylo_tree>{};
//   for (auto i = 0; i != num_trees_; ++i) {
//     base_trees.push_back(base_tree_by_index(first_base_tree_index_ + i));
//   }

//   const auto& master_base_tree = base_tree_by_index(master_base_tree_index);

//   auto mcc_tree = Mcc_tree{master_base_tree.size(), std::move(base_trees), master_base_tree_index};
//   copy_topology(master_base_tree, mcc_tree);

//   for (const auto& [clade, node_mappings] : clade_base_nodes_) {
//     const auto base_tree_node_mapping_it = node_mappings.find(master_base_tree_index);
//     if (base_tree_node_mapping_it != node_mappings.end()) {
//       auto translated_node_mappings = std::unordered_map<int, Node_index>{};
//       if (first_base_tree_index_ == 0) {
//         translated_node_mappings = node_mappings;
//       } else {
//         for (const auto [orig_base_tree_index, node_index] : node_mappings) {
//           assert(orig_base_tree_index >= first_base_tree_index_);
//           assert(orig_base_tree_index < first_base_tree_index_ + num_trees_);
//           translated_node_mappings[orig_base_tree_index - first_base_tree_index_] = node_index;
//         }
//       }
//       const auto [_, base_tree_node_index] = *base_tree_node_mapping_it;
//       auto mcc_node = Mut_node_ref{mcc_tree, base_tree_node_index};
//       mcc_node->base_tree_mappings().clear();
//       for (const auto& [tree_index, node_index] : translated_node_mappings) {
//         mcc_node->base_tree_mappings().push_back({
//             .base_tree_index = tree_index,
//             .base_tree_node_index = node_index,
//             .is_monophyletic_in_base_tree = true // TODO: Detect non-monophyletic mappings in Incremental_mcc_tree
//           });
//       }
//     }
//   }

//   mcc_tree.calculate_derived_quantities();

//   return mcc_tree;
// }

//}  // namespace delphy
