#ifndef DELPHY_TREE_PARTITIONING_H_
#define DELPHY_TREE_PARTITIONING_H_

#include <random>

#include <absl/random/bit_gen_ref.h>

#include "tree.h"

namespace delphy {

// To partition a tree, we first choose a set of N cut points.  Each cut point P is the root of a
// subtree T_P that forms a member of the partition.  The subtree T_P contains all descendants of the cut point
// that are not themselves a cut point or a descendant thereof.  If P' is a cut point with closest ancestor
// cut point is P, then P' appears as a tip in T_P and as the root of T_{P'}.

// For a partition of a tree T, stores the information about the part of the partition rooted at
// the `cut_point` in T.
struct Partition_part_info {
  Node_index cut_point;
};

// Generates a random partition stencil of tree `tree` into at most `num_parts` parts (possibly fewer).  A best effort
// is made to make the parts roughly of equal size.
auto generate_random_partition_stencil(
    const Tree_like auto& t,
    int num_parts,
    absl::BitGenRef bitgen)
    -> std::vector<Partition_part_info>;

class Partition_part_node : public Binary_node {
 public:
  auto orig_tree_index() const -> Node_index { return orig_tree_index_; }
  auto set_orig_tree_index(Node_index orig_tree_index) -> void { orig_tree_index_ = orig_tree_index; }

 private:
  Node_index orig_tree_index_{k_no_node};
};

// A part of a partition.  Includes a map from node indices in the subtree to the node indices in the (whole) tree
class Partition_part : public Tree<Partition_part_node> {
 public:
  Partition_part() : Tree{0}, info_{} {}
  explicit Partition_part(Partition_part_info info) : Tree{0}, info_{info} {}

  auto info() const -> Partition_part_info { return info_; }

 private:
  Partition_part_info info_;
};
static_assert(Tree_like<Partition_part>);

// A partition.  Consists of many parts, one of which contains the root of the original tree
class Partition {
 public:
  Partition() : Partition(0) {}
  explicit Partition(int num_parts)
      : parts_(num_parts, Partition_part{}),
        root_part_index_{num_parts},
        tree_index_to_partition_index_{} {}

  auto parts() -> std::vector<Partition_part>& { return parts_; }
  auto parts() const -> const std::vector<Partition_part>& { return parts_; }
  auto root_part_index() const -> int { return root_part_index_; }
  auto set_root_part_index(int root_part_index) -> void { root_part_index_ = root_part_index; }
  auto tree_index_to_partition_index() -> Node_vector<int>& { return tree_index_to_partition_index_; }
  auto tree_index_to_partition_index() const -> const Node_vector<int>& { return tree_index_to_partition_index_; }

 private:
  std::vector<Partition_part> parts_;
  int root_part_index_;
  Node_vector<int> tree_index_to_partition_index_;
};

// Given a list of partition part infos describing a set of cut points and the size of the parts of the
// induced partition, creates the corresponding subtrees and maps to the original tree
auto partition_tree(const Tree_like auto& tree, const std::vector<Partition_part_info>& part_infos) -> Partition;

// Reassembles a tree from a given partition
auto reassemble_tree(const Partition& partition) -> Tree<Binary_node>;


// Implementation
// ==============

namespace details {

auto make_partition_part(
    const Tree_like auto& src_tree,
    const Node_index cut_point,
    const Node_set& src_cut_points)
    -> Partition_part {
  
  auto part = Partition_part{Partition_part_info{cut_point}};

  struct Work_item {
    Node_index src_node;
    Node_index dst_node;
  };
  auto work_stack = std::stack<Work_item>{};

  auto root_dst_node = part.add_node();
  part.root = root_dst_node;
  work_stack.push({.src_node = cut_point, .dst_node = root_dst_node});
  
  while (not work_stack.empty()) {
    // Every node in the destination tree will appear exactly once in this loop as cur.dst_node
    auto [src_node, dst_node] = work_stack.top();
    work_stack.pop();

    part.at(dst_node).set_orig_tree_index(src_node);
    
    auto src_node_is_cut_point = src_cut_points.contains(src_node);
    if (src_tree.at(src_node).is_tip() || (src_node_is_cut_point && src_node != cut_point)) {
      // dst_node will be a tip in this subtree
      part.at(dst_node).children = {};
    } else {
      // dst_node is an inner node in this subtree
      DCHECK_EQ(std::ssize(src_tree.at(src_node).children), 2);

      auto dst_left = part.add_node();
      auto dst_right = part.add_node();
      part.at(dst_node).children = {dst_left, dst_right};
      part.at(dst_left).parent = dst_node;
      part.at(dst_right).parent = dst_node;
      
      work_stack.push({.src_node = src_tree.at(src_node).left_child(), .dst_node = dst_left});
      work_stack.push({.src_node = src_tree.at(src_node).right_child(), .dst_node = dst_right});
    }
  }

  assert_tree_integrity(part);

  return part;
}

}  // namespace details

auto generate_random_partition_stencil(
    const Tree_like auto& tree,
    int num_parts,
    absl::BitGenRef bitgen)
    -> std::vector<Partition_part_info> {
  
  auto part_infos = std::vector<Partition_part_info>{};
  part_infos.reserve(num_parts);

  // As we proceed, descendants[i] will count the number of descendants of i that aren't already assigned a part_index.
  auto descendants = Node_vector<int>(std::ssize(tree), 0);

  auto num_branches_left = std::ssize(tree);
  auto num_parts_left = num_parts;
  for (auto node : randomized_post_order_traversal(tree, bitgen)) {  // insert a bit of randomness in traversal

    // The root never goes explicitly into the stencils
    if (node == tree.root) { break; }
    
    // Stop if we've already cut into num_parts-1 parts (the last one will implicitly start at the root)
    if (std::ssize(part_infos) == (num_parts - 1)) { break; }
    
    // Accumulate descendants
    descendants[node] = 1;
    for (auto child : tree.at(node).children) {
      descendants[node] += descendants[child];
    }

    // Should we cut here?
    auto min_subtree_size = std::max(10L, num_branches_left / (num_parts_left + 1));  // +1 => a little bit of slack
    if (descendants[node] >= min_subtree_size) {

      // Can we cut here?
      auto is_allowed = true;
      
      // 1. Don't cut if the remaining stump would be too small
      if (is_allowed && (num_branches_left - (descendants[node] - 1)) < min_subtree_size) { is_allowed = false; }

      // 2. A bit of randomness
      if (is_allowed && std::bernoulli_distribution{0.5}(bitgen)) { is_allowed = false; }

      if (is_allowed) {
        // Cut here!
        auto cut_point = node;
        auto num_branches_in_part = descendants[cut_point] - 1;
        num_branches_left -= num_branches_in_part;
        CHECK_GE(num_branches_left, 0);
        part_infos.push_back({cut_point});
        descendants[cut_point] = 1;
        --num_parts_left;
      }
    }
  }

  return part_infos;
}

auto partition_tree(const Tree_like auto& tree, const std::vector<Partition_part_info>& stencil) -> Partition {
  // The root partition is implicitly the last one in the stencil if it's not explicitly in the stencil
  auto root_part_index = static_cast<int>(std::ssize(stencil));
  auto root_in_stencil = false;
  for (auto i = 0; i != std::ssize(stencil); ++i) {
    const auto& p = stencil[i];
    if (p.cut_point == tree.root) {
      root_in_stencil = true;
      root_part_index = i;
      break;
    }
  }
  
  auto num_partitions = static_cast<int>(std::ssize(stencil) + (root_in_stencil ? 0 : 1));
  
  auto partition = Partition{num_partitions};

  auto partition_cut_points = Node_set{};
  partition_cut_points.reserve(num_partitions);
  for (const auto& p : stencil) {
    partition_cut_points.insert(p.cut_point);
  }

  if (not root_in_stencil) {
    partition_cut_points.insert(tree.root);  // Implicit final partition
    root_part_index = num_partitions - 1;
  }
  partition.set_root_part_index(root_part_index);

  for (auto i = 0; i != num_partitions; ++i) {
    auto cut_point = (i == partition.root_part_index()) ? tree.root : stencil[i].cut_point;
    partition.parts()[i] = details::make_partition_part(tree, cut_point, partition_cut_points);
  }

  partition.tree_index_to_partition_index().assign(std::ssize(tree), -1);
  for (auto i = 0; i != num_partitions; ++i) {
    const auto& partition_part = partition.parts()[i];
    for (const auto& pnode : index_order_traversal(partition_part)) {
      partition.tree_index_to_partition_index()[partition_part.at(pnode).orig_tree_index()] = static_cast<int>(i);
    }
  }

  return partition;
}

}  // namespace delphy

#endif // DELPHY_TREE_PARTITIONING_H_
