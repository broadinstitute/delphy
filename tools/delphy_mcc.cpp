#include <iostream>

#include "absl/log/initialize.h"
#include "absl/log/check.h"
#include "absl/random/random.h"

#include "beasty_input.h"
#include "mcc_tree.h"
#include "version.h"

namespace delphy {

auto output_mcc_tree(const Mcc_tree& mcc_tree, std::ostream& os) -> void {
  auto mcc_node_name = [&](Node_index mcc_node) -> const std::string& {
    CHECK(mcc_tree.at(mcc_node).is_tip());  // Otherwise, "node name" can be different in each base tree
    const auto& base_tree = *mcc_tree.base_trees().front();
    auto base_node = corresponding_node_to(mcc_tree, mcc_node, 0);
    return base_tree.at(base_node).name;
  };
  
  auto num_tips = (std::ssize(mcc_tree) + 1) / 2;
  auto next_tip_index = Node_index{1};
  auto node_to_tip_index = Node_vector<Node_index>(std::ssize(mcc_tree), k_no_node);
  for (const auto& node : index_order_traversal(mcc_tree)) {
    if (mcc_tree.at(node).is_tip()) {
      node_to_tip_index[node] = next_tip_index;
      ++next_tip_index;
    }
  }
  
  os << "#NEXUS\n"
     << "\n"
     << "Begin taxa;\n"
     << "    Dimensions ntax=" << num_tips << ";\n"
     << "        Taxlabels\n";
  
  for (const auto& node : index_order_traversal(mcc_tree)) {
    if (mcc_tree.at(node).is_tip()) {
      os << "            " << mcc_node_name(node) << "\n";
    }
  }

  os << "            ;\n"
     << "End;\n"
     << "Begin trees;\n"
     << "    Translate\n";

  auto first_tip = true;
  for (const auto& node : index_order_traversal(mcc_tree)) {
    if (mcc_tree.at(node).is_tip()) {
      if (first_tip) {
        first_tip = false;
      } else {
        os << ",\n";
      }
      os << absl::StrFormat("        %4d %s", node_to_tip_index[node], mcc_node_name(node));
    }
  }
  os << "\n";

  os << ";\n";

  os << "tree TREE1 = ";

  for (const auto& [node, children_so_far] : traversal(mcc_tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      if (mcc_tree.at(node).is_inner_node()) {
        os << "(";
      }
    }

    if (0 < children_so_far && children_so_far < std::ssize(mcc_tree.at(node).children)) {
      // In between two children
      os << ",";
    }

    if (children_so_far == std::ssize(mcc_tree.at(node).children)) {
      // Exiting `node`
      if (mcc_tree.at(node).is_inner_node()) {
        os << ")";
      }
      if (mcc_tree.at(node).is_tip()) {
        os << node_to_tip_index.at(node);
      }
      
      auto posterior = mcc_tree.at(node).posterior_support();
      os << "[&posterior=" << posterior << "]";
      
      auto days = node == mcc_tree.root
          ? 0.0
          : mcc_tree.at(node).t_mrca() - mcc_tree.at_parent_of(node).t_mrca();
      os << ":" << (days / 365.0);  // Days (internal) -> Years (output)
    }
  }
  
  os << ";\n";
  os << "End;\n";
}

}  // namespace delphy

auto main(int argc, char** argv) -> int {
  using namespace delphy;

  absl::InitializeLog();
  
  auto scope = Local_arena_scope{};

  if (argc != 3) {
    std::cerr << absl::StreamFormat("Delphy MCC Version %s (build %d, commit %s)",
                                    k_delphy_version_string,
                                    k_delphy_build_number,
                                    k_delphy_commit_string) << "\n\n";
    std::cerr << "Usage: delphy_mcc <input-trees-file> <output-mcc-file>\n\n";
    std::cerr << "Reads a BEAST-format posterior trees files and computes a\n"
              << "Maximum Clade Credibility (MCC) tree.\n"
              << "A 30% burn-in is applied automatically.\n";
    return EXIT_FAILURE;
  }

  auto tree_is = std::ifstream{argv[1]};
  CHECK(tree_is);

  auto min_state = 0L;
  auto trees = read_beasty_trees(tree_is, min_state);
  tree_is.close();

  std::cerr << absl::StreamFormat("Read %d trees", std::ssize(trees)) << std::endl;

  min_state = trees.begin()->first;
  auto max_state = trees.rbegin()->first;

  auto max_burnin_state = min_state + ((max_state - min_state) * 3) / 10;  // Hard-coded 30% burn-in
  
  auto num_burnin_trees = 0;
  for (auto i = trees.begin();
       i != trees.end() && i->first <= max_burnin_state;
       i = trees.erase(i)) {
    ++num_burnin_trees;
  }
  std::cerr << absl::StreamFormat("Discarded %d trees as burn-in", num_burnin_trees) << std::endl;
  std::cerr << absl::StreamFormat("Read %d post-burnin trees", std::ssize(trees)) << std::endl;

  auto base_trees = std::vector<Phylo_tree*>{};
  for (auto& [state, tree] : trees) {
    base_trees.push_back(&tree);
  }

  std::cerr << absl::StreamFormat("Building MCC with %d trees", std::ssize(base_trees)) << std::endl;

  auto bitgen = absl::BitGen{};
  auto mcc_tree = derive_mcc_tree(std::move(base_trees), bitgen);

  auto mcc_os = std::ofstream{argv[2]};
  CHECK(mcc_os);
  output_mcc_tree(mcc_tree, mcc_os);
  
  return 0;
}
