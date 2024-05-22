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
  Node_index num_nodes;  // 1 more than the number of branches
};

// Generates a random partition of tree `tree` into at most `num_parts` parts (possibly fewer).  A best effort
// is made to make the parts roughly of equal size.  If `forbidden_cut_points` is specified and not
// empty, then those points will not be chosen as cut points for the partition (see above).  Useful to
// generate a sequence of random partitions that eventually include every inner node of `tree` as an
// inner node of some part.
auto generate_random_partition(
    const Tree_like auto& t,
    int num_parts,
    absl::BitGenRef bitgen,
    const Node_set& forbidden_cut_points)
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

// Picks a random path starting at the root and report the node with number of descendants closest to target_branches.
// The random path is biased towards children with more descendants, to avoid going down rabbit holes.
// Unless we choose the root, the chosen cut point will never be one with an index in forbidden_cut_node_indices
auto choose_cut_point(
    const Tree_like auto& tree,
    const Node_vector<int>& descendants,
    int target_branches,
    absl::BitGenRef bitgen,
    const Node_set& forbidden_cut_points)
    -> Node_index {
  
  auto best = tree.root;
  auto best_score = std::abs(target_branches - descendants[best]);

  auto cur = tree.root;
  while (descendants[cur] != 0) {
    auto cur_descendants = descendants[cur];
    auto cur_score = std::abs(target_branches - cur_descendants);
    auto is_allowed = not forbidden_cut_points.contains(cur);
    if (is_allowed && cur_score < best_score) {
      best = cur;
      best_score = cur_score;
    }
    if (cur_descendants < target_branches && cur_score >= best_score) {
      // Score of descendants can only be worse as we go further down the tree; no point in going further
      break;
    }

    assert(not tree.at(cur).children.empty());
    auto left = tree.at(cur).left_child();
    auto right = tree.at(cur).right_child();

    auto w_left = std::pow(descendants[left], 2.0);
    auto w_right = std::pow(descendants[right], 2.0);
    auto p_left = w_left / (w_left + w_right);

    cur = std::uniform_real_distribution{0.0, 1.0}(bitgen) < p_left ? left : right;
  }

  return best;
}

auto discount_descendants_of(const Tree_like auto& tree, Node_index node, Node_vector<int>& descendants) -> void {
  auto num_gone = descendants[node];
  for (auto cur = node; cur != k_no_node; cur = tree.at(cur).parent) {
    descendants[cur] -= num_gone;
  }
}

auto make_partition_part(
    const Tree_like auto& src_tree,
    const Partition_part_info& part_info,
    const Node_set& src_cut_points)
    -> Partition_part {
  
  auto part = Partition_part{part_info};

  struct Work_item {
    Node_index src_node;
    Node_index dst_node;
  };
  auto work_stack = std::stack<Work_item>{};

  auto root_dst_node = part.add_node();
  part.root = root_dst_node;
  work_stack.push({.src_node = part_info.cut_point, .dst_node = root_dst_node});
  
  while (not work_stack.empty()) {
    // Every node in the destination tree will appear exactly once in this loop as cur.dst_node
    auto [src_node, dst_node] = work_stack.top();
    work_stack.pop();

    part.at(dst_node).set_orig_tree_index(src_node);
    
    auto src_node_is_cut_point = src_cut_points.contains(src_node);
    if (src_tree.at(src_node).is_tip() || (src_node_is_cut_point && src_node != part_info.cut_point)) {
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

auto generate_random_partition(
    const Tree_like auto& tree,
    int num_parts,
    absl::BitGenRef bitgen,
    const Node_set& forbidden_cut_points)
    -> std::vector<Partition_part_info> {
  
  auto part_infos = std::vector<Partition_part_info>{};
  part_infos.reserve(num_parts);

  // As we proceed, descendants[i] is the number of descendants of i that aren't already assigned a part_index.
  // For efficiency, if descendants[i] == 0, then we don't update or use descendants[j] for all descendants j of i
  auto descendants = count_all_descendants(tree);

  auto num_branches_left = descendants[tree.root];
  for (auto num_parts_left = num_parts; num_parts_left > 1 && num_branches_left > 0; --num_parts_left) {

    auto ideal_target_branches = num_branches_left / num_parts_left;
    auto min_target_branches = std::max(static_cast<int>(std::floor(ideal_target_branches * 0.7)), 1);
    auto max_target_branches = std::min(static_cast<int>(std::ceil(ideal_target_branches * 1.3)), num_branches_left);

    auto target_branches =
        (min_target_branches >= max_target_branches)
        ? num_branches_left
        : std::uniform_int_distribution{min_target_branches, max_target_branches}(bitgen);

    auto cut_point = details::choose_cut_point(tree, descendants, target_branches, bitgen, forbidden_cut_points);
    auto num_branches_in_part = descendants[cut_point];
    num_branches_left -= num_branches_in_part;
    details::discount_descendants_of(tree, cut_point, descendants);

    part_infos.push_back({cut_point, 1 + num_branches_in_part});
  }

  // All remaining branches are lumped into a part_index rooted at the tree root
  if (num_branches_left > 0) {
    part_infos.push_back({tree.root, 1 + num_branches_left});
  }

  return part_infos;
}

auto partition_tree(const Tree_like auto& tree, const std::vector<Partition_part_info>& part_infos) -> Partition {
  DCHECK_EQ(std::ranges::count_if(part_infos,
                                  [&](const auto& part_info) { return part_info.cut_point == tree.root; }),
            1);

  auto partition = Partition{static_cast<int>(std::ssize(part_infos))};
  partition.parts().reserve(std::ssize(part_infos));

  auto partition_cut_points = Node_set{};
  partition_cut_points.reserve(std::ssize(part_infos));
  for (const auto& p : part_infos) {
    partition_cut_points.insert(p.cut_point);
  }

  auto root_part_index = std::ssize(part_infos);
  for (auto i = 0; i != std::ssize(part_infos); ++i) {
    const auto& p = part_infos[i];
    partition.parts()[i] = details::make_partition_part(tree, p, partition_cut_points);
    if (p.cut_point == tree.root) {
      root_part_index = i;
    }
  }
  CHECK_NE(root_part_index, std::ssize(part_infos));
  partition.set_root_part_index(root_part_index);

  partition.tree_index_to_partition_index().assign(std::ssize(tree), -1);
  for (auto i = 0; i != std::ssize(part_infos); ++i) {
    const auto& partition_part = partition.parts()[i];
    for (const auto& pnode : index_order_traversal(partition_part)) {
      partition.tree_index_to_partition_index()[partition_part.at(pnode).orig_tree_index()] = static_cast<int>(i);
    }
  }

  return partition;
}

}  // namespace delphy

#endif // DELPHY_TREE_PARTITIONING_H_
