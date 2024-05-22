#include "mcc_tree.h"

#include <cmath>
#include <random>

#include "absl/hash/hash.h"

namespace delphy {

class Clade_fingerprint {
 public:
  Clade_fingerprint() : data_{0} {}
  explicit Clade_fingerprint(absl::BitGenRef bitgen)
      : data_{std::uniform_int_distribution<uint64_t>{}(bitgen)} {}
  Clade_fingerprint(Clade_fingerprint left, Clade_fingerprint right)
      : data_{left.data_ ^ right.data_} {}

  auto operator<=>(const Clade_fingerprint& rhs) const = default;

  template <typename H>
  friend H AbslHashValue(H h, const Clade_fingerprint& f) {
    return H::combine(std::move(h), f.data_);
  }

 private:
  uint64_t data_;
};

// Given clade fingerprints in the tips, propagate them up to inner nodes
auto calc_inner_node_clade_fingerprints(
    const Phylo_tree& tree,
    Node_vector<Clade_fingerprint>& clade_fingerprints)  // in-out
    -> void {
  for (const auto& node : post_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      clade_fingerprints[node] = Clade_fingerprint{
        clade_fingerprints[tree.at(node).left_child()],
        clade_fingerprints[tree.at(node).right_child()]};
    }
  }
}

static auto count_base_trees_for_each_clade(
    const Base_tree_vector<Phylo_tree*>& base_trees,
    Node_vector<Clade_fingerprint>& clade_fingerprints)
    -> absl::flat_hash_map<Clade_fingerprint, int>
{
  auto result = absl::flat_hash_map<Clade_fingerprint, int>{};  // counts implicitly start at 0
  for (const auto* base_tree : base_trees) {
    calc_inner_node_clade_fingerprints(*base_tree, clade_fingerprints);
    for (auto clade_fingerprint : clade_fingerprints) { // loop over clade fingerprints of *all* nodes
      result[clade_fingerprint] += 1;
    }
  }
  return result;
}

auto derive_mcc_tree(Base_tree_vector<Phylo_tree*> base_trees, absl::BitGenRef bitgen) -> Mcc_tree {
  auto M = std::ssize(base_trees);  // M = number of base trees
  CHECK_GT(M, 0);
  
  auto num_nodes = Node_index(std::ssize(*base_trees[0]));
  if (estd::is_debug_enabled) {
    for (const auto* base_tree : base_trees) {
      CHECK_EQ(std::ssize(*base_tree), num_nodes);
    }
  }

  // Generate random clade fingerprints for tips
  auto clade_fingerprints = Node_vector<Clade_fingerprint>(num_nodes, Clade_fingerprint{});
  for (const auto& node : index_order_traversal(*base_trees[0])) {
    if (base_trees[0]->at(node).is_tip()) {
      clade_fingerprints[node] = Clade_fingerprint{bitgen};
    }
  }

  // Pre-calc log(i / M) for i in [0,M]
  auto log_i_over_M = std::vector<double>{};
  log_i_over_M.reserve(M + 1);
  log_i_over_M.push_back(-std::numeric_limits<double>::infinity());  // log(0)
  auto log_M = std::log(M);
  for (auto i = 1; i <= M; ++i) {
    log_i_over_M.push_back(std::log(i) - log_M);
  }

  // Calculate log(clade_credibility) for each base tree
  auto base_tree_log_cc = std::vector<double>(M, 0.0);
  auto clade_base_tree_counts = count_base_trees_for_each_clade(base_trees, clade_fingerprints);
  for (auto i = 0; i != M; ++i) {
    const auto& base_tree = *base_trees[i];
    
    calc_inner_node_clade_fingerprints(base_tree, clade_fingerprints);
    for (const auto& node : index_order_traversal(base_tree)) {
      if (not base_tree.at(node).is_tip()) {  // All tips contribute equally to log(cc) in each base tree
        auto clade_fingerprint = clade_fingerprints[node];
        base_tree_log_cc[i] += log_i_over_M[clade_base_tree_counts[clade_fingerprint]];
      }
    }
  }

  // Identify master base tree as the one with maximum clade credibility
  auto master_base_tree_index = static_cast<int>(std::distance(
      begin(base_tree_log_cc),
      std::max_element(begin(base_tree_log_cc), end(base_tree_log_cc))));
  const auto& master_base_tree = *base_trees[master_base_tree_index];

  // Set up MCC tree using master base tree as template
  // NOTE: `base_trees` is moved here, use mcc_tree.base_trees() from here onwards
  auto mcc_tree = Mcc_tree{num_nodes, std::move(base_trees), master_base_tree_index};
  copy_topology(master_base_tree, mcc_tree);

  auto mcc_tree_clade_fingerprints = Node_vector<Clade_fingerprint>{clade_fingerprints};  // copy
  calc_inner_node_clade_fingerprints(master_base_tree, mcc_tree_clade_fingerprints);

  // Node i in MCC represents the MRCA of all the tips below i in the MCC.
  // In the loop over `base_tree`s below, the MRCA of those tips in `base_tree`
  // is the node mrca_in_base_of_tips_below_mcc_node[i]
  auto corresponding_nodes = Node_vector<Node_index>(num_nodes, k_no_node);
  for (const auto& mcc_node : index_order_traversal(mcc_tree)) {
    if (mcc_tree.at(mcc_node).is_tip()) {
      // For tips: 1-1 correspondence of tips, always
      // For inner nodes: correponding node depends on base tree (see below)
      corresponding_nodes[mcc_node] = mcc_node;
    }
  }
  
  for (auto i = 0; i != M; ++i) {
    const auto& base_tree = *mcc_tree.base_trees()[i];

    calc_inner_node_clade_fingerprints(base_tree, clade_fingerprints);
    
    // Map each node in the MCC to the MRCA of its children in the base tree
    for (const auto& mcc_node : post_order_traversal(mcc_tree)) {
      if (mcc_tree.at(mcc_node).is_tip()) {
        //corresponding_nodes[mcc_node] = mcc_node;   // <-- Set up above!
      } else {
        corresponding_nodes[mcc_node] = find_MRCA_of(base_tree, 
                                                     corresponding_nodes[mcc_tree.at(mcc_node).left_child()],
                                                     corresponding_nodes[mcc_tree.at(mcc_node).right_child()]);
      }

      const auto& clade_fingerprint = clade_fingerprints[corresponding_nodes[mcc_node]];
      mcc_tree.at(mcc_node).corresponding_node_infos().push_back({
          .node_in_base_tree = corresponding_nodes[mcc_node],
          .is_monophyletic_in_base_tree = (clade_fingerprint == mcc_tree_clade_fingerprints[mcc_node])
        });
    }
  }

  assert_tree_integrity(mcc_tree);

  mcc_tree.calculate_derived_quantities();

  return mcc_tree;
}

auto Mcc_tree::calculate_derived_quantities() -> void {
  for (const auto& mcc_node : index_order_traversal(*this)) {
    // Find mean time of mapped nodes (only monophyletic ones, and all of them too)
    auto sum_t = 0.0;
    auto num_exact_matches = 0;
    auto sum_t_mrca = 0.0;
    for (auto i = 0; i != num_base_trees(); ++i) {
      const auto& base_tree = *base_trees().at(i);
      const auto& corresponding_node_info = at(mcc_node).corresponding_node_infos().at(i);
      auto corresponding_node = corresponding_node_info.node_in_base_tree;
      sum_t_mrca += base_tree.at(corresponding_node).t;
      if (corresponding_node_info.is_monophyletic_in_base_tree) {
        sum_t += base_tree.at(corresponding_node).t;
        ++num_exact_matches;
      }
    }

    CHECK_GT(num_exact_matches, 0);  // Corresponding node is always monophyletic in master base tree
    at(mcc_node).set_t(sum_t / num_exact_matches);
    at(mcc_node).set_t_mrca(sum_t_mrca / num_base_trees());
  }
}

}  // namespace delphy
