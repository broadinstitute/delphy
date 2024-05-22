#include "tree_partitioning.h"

#include <algorithm>

namespace delphy {

class Partition_reassembler {
  using Dst_tree = Tree<Binary_node>;

 public:
  Partition_reassembler(Dst_tree& dst_t,
                        const Partition& partition,
                        const Node_map<int>& cut_point_to_part)
      : dst_t_{&dst_t}, partition_{&partition}, cut_point_to_part_{&cut_point_to_part} {}

  auto add_part(const Partition_part& part) -> Node_index {
    return add_node(part, part.root);
  }

  auto add_node(const Partition_part& src, Node_index src_node) -> Node_index {
    auto dst_node = src.at(src_node).orig_tree_index();

    if (src.at(src_node).children.empty()) {
      // It's either a real tip, or a placeholder for the root of another partition
      auto it = cut_point_to_part_->find(dst_node);
      if (it == cut_point_to_part_->end()) {
        // The common case: just a normal tip
        dst_t_->at(dst_node).children = {};
      } else {
        // Bring in a whole new tree
        auto [_, new_part] = *it;
        [[maybe_unused]] auto new_tree_root = add_part(partition_->parts()[new_part]);
        assert(new_tree_root == dst_node);
      }

    } else {
      // Inner node
      auto dst_left = add_node(src, src.at(src_node).left_child());
      auto dst_right = add_node(src, src.at(src_node).right_child());

      dst_t_->at(dst_node).children = {dst_left, dst_right};
      dst_t_->at(dst_left).parent = dst_node;
      dst_t_->at(dst_right).parent = dst_node;
    }

    return dst_node;
  }

 private:
  Tree<Binary_node>* dst_t_;
  const Partition* partition_;
  const Node_map<int>* cut_point_to_part_;
};

auto reassemble_tree(const Partition& partition) -> Tree<Binary_node> {

  // Map each cut point to a part
  auto cut_point_to_part = Node_map<int>{};
  cut_point_to_part.reserve(partition.parts().size());
  for (auto i = 0; i != std::ssize(partition.parts()); ++i) {
    const auto& part = partition.parts()[i];
    cut_point_to_part.insert({part.info().cut_point, i});
  }

  // Set up new empty tree
  auto result = Tree<Binary_node>{};
  if (partition.parts().empty()) {
    return result;
  }

  // Reassemble all subtrees back to make original tree, starting with the partition containing the
  // original tree's root
  auto reassembler = Partition_reassembler{result, partition, cut_point_to_part};
  const auto& root_part = partition.parts()[partition.root_part_index()];
  auto new_root = reassembler.add_part(root_part);
  assert(root_part.at_root().orig_tree_index() == new_root);
  result.at(new_root).parent = k_no_node;
  result.root = new_root;

  assert_tree_integrity(result);

  return result;
}

}  // namespace delphy
