#include "tree.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

namespace delphy {

TEST(Tree_test, nary_empty) {
  auto tree = Tree<Nary_node>{};

  EXPECT_THAT(tree, testing::IsEmpty());
  EXPECT_THAT(tree.root, testing::Eq(k_no_node));
}

TEST(Tree_test, nary_simple) {
  auto tree = Tree<Nary_node>{2};
  auto i = 0;
  auto j = 1;
  auto k = tree.add_node();

  tree.at(i).children = {j, k};
  tree.at(j).parent = i;
  tree.at(k).parent = i;
  
  tree.root = i;

  assert_tree_integrity(tree);

  EXPECT_THAT(tree.at(i).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.at(i).children, testing::ElementsAre(1, 2));
}

TEST(Tree_test, binary_empty) {
  auto tree = Tree<Binary_node>{};

  EXPECT_THAT(tree, testing::IsEmpty());
  EXPECT_THAT(tree.root, testing::Eq(k_no_node));
}

TEST(Tree_test, binary_simple) {
  auto tree = Tree<Binary_node>{2};
  auto i = 0;
  auto j = 1;
  auto k = tree.add_node();
  
  tree.at(i).children = {j, k};
  tree.at(j).parent = i;
  tree.at(k).parent = i;
  
  tree.root = i;

  assert_tree_integrity(tree);

  EXPECT_THAT(tree.at(i).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.at(i).children, testing::ElementsAre(1, 2));

  EXPECT_THAT(tree.at(i).left_child(), testing::Eq(j));
  EXPECT_THAT(tree.at(i).right_child(), testing::Eq(k));
  EXPECT_THAT(tree.at(j).parent, testing::Eq(i));
  EXPECT_THAT(tree.at(k).parent, testing::Eq(i));
  EXPECT_THAT(tree.at(i).sibling_of(j), testing::Eq(k));
  EXPECT_THAT(tree.at(i).sibling_of(k), testing::Eq(j));
}

TEST(Tree_test, meaty_simple) {
  struct Test_node : public Nary_node {
    std::string name{};
    double t{};
  };

  auto tree = Tree<Test_node>{2};
  auto i = 0;
  auto j = 1;
  tree.at(i).name = "Hello";
  tree.at(i).t = 5.0;
  tree.at(j).name = "World";
  tree.at(j).t = 10.0;

  tree.root = i;
  tree.at(i).children = {j};
  tree.at(j).parent = i;

  assert_tree_integrity(tree);

  EXPECT_THAT(tree.at(tree.at_root().children[0]).name, testing::StrEq("World"));
}

TEST(Tree_test, empty_traversals) {
  struct Test_node : public Binary_node {
    std::string str{};
  };
  auto tree = Tree<Test_node>{};
  
  for (const auto& node : in_order_traversal(tree)) {
    FAIL() << "Tree is empty" << node;
  }
  
  for (const auto& node : pre_order_traversal(tree)) {
    FAIL() << "Tree is empty" << node;
  }
  
  for (const auto& node : post_order_traversal(tree)) {
    FAIL() << "Tree is empty" << node;
  }
  
  for (const auto& node : index_order_traversal(tree)) {
    FAIL() << "Tree is empty" << node;
  }
  
  for (const auto& [node, children_so_far] : traversal(tree)) {
    FAIL() << "Tree is empty" << node << children_so_far;
  }
}

TEST(Tree_test, typical_traversals) {
  struct Test_node : public Binary_node {
    std::string str{};
  };
  auto tree = Tree<Test_node>{};
  auto c = tree.add_node();
  tree.at(c).str = "c";
  auto b = tree.add_node();
  tree.at(b).str = "b";
  auto a = tree.add_node();
  tree.at(a).str = "a";

  tree.root = a;
  tree.at(a).children = {b, c};
  tree.at(b).parent = a;
  tree.at(c).parent = a;

  assert_tree_integrity(tree);

  auto result = std::vector<std::string>{};

  result.clear();
  for (const auto& node : in_order_traversal(tree)) { result.push_back(tree.at(node).str); }
  EXPECT_THAT(result, testing::ElementsAre("b", "a", "c"));

  result.clear();
  for (const auto& node : pre_order_traversal(tree)) { result.push_back(tree.at(node).str); }
  EXPECT_THAT(result, testing::ElementsAre("a", "b", "c"));

  result.clear();
  for (const auto& node : post_order_traversal(tree)) { result.push_back(tree.at(node).str); }
  EXPECT_THAT(result, testing::ElementsAre("b", "c", "a"));

  result.clear();
  for (const auto& node : index_order_traversal(tree)) { result.push_back(tree.at(node).str); }
  EXPECT_THAT(result, testing::ElementsAre("c", "b", "a"));
  
  result.clear();
  for (const auto& [node, children_so_far] : traversal(tree)) {
    result.push_back(absl::StrFormat("(%s, %d)", tree.at(node).str, children_so_far));
  }
  EXPECT_THAT(result, testing::ElementsAre(
      "(a, 0)",
      "(b, 0)",
      "(a, 1)",
      "(c, 0)",
      "(a, 2)"));
}

}  // namespace delphy
