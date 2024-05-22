#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "api.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

static auto build_simple_tree() -> Phylo_tree {
  // The tree that we build
  //
  // Time:             -1.0          0.0        1.0        2.0        3.0
  //                   
  //                                  +-- T0C -- a (CANN)
  //                                  |
  //                     +A2N- A0T ---+ x (TANN)
  //                     |            |
  // (AACA) A3N- C2A --+ r (AAAN)     +-------- A1G ------- b (TGNN)
  //                     |
  //                     +A1N------------------ A0G ------------------ c (GNAN)

  auto ref_sequence = Real_sequence{rA, rA, rC, rA};
  auto tree = Phylo_tree{5};
  
  auto r = Node_index{0};
  auto x = Node_index{1};
  auto a = Node_index{2};
  auto b = Node_index{3};
  auto c = Node_index{4};
  
  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {x, c};
  tree.at(r).name = "r";
  tree.at(r).t = -1.0;
  tree.at(r).mutations = {Mutation{rC, 2, rA, -std::numeric_limits<double>::max()}};
  tree.at(r).missations = {{Missation{3, rA}}, ref_sequence};
  
  tree.at(x).parent = r;
  tree.at(x).children = {a, b};
  tree.at(x).name = "x";
  tree.at(x).t = 0.0;
  tree.at(x).mutations = {Mutation{rA, 0, rT, -0.5}};
  tree.at(x).missations = {{Missation{2, rA}}, ref_sequence};
  
  tree.at(a).parent = x;
  tree.at(a).children = {};
  tree.at(a).name = "a";
  tree.at(a).t = 1.0;
  tree.at(a).mutations = {Mutation{rT, 0, rC, 0.5}};
  
  tree.at(b).parent = x;
  tree.at(b).children = {};
  tree.at(b).name = "b";
  tree.at(b).t = 2.0;
  tree.at(b).mutations = {Mutation{rA, 1, rG, 1.0}};
  
  tree.at(c).parent = r;
  tree.at(c).children = {};
  tree.at(c).name = "c";
  tree.at(c).t = 3.0;
  tree.at(c).mutations = {Mutation{rA, 0, rG, 1.0}};
  tree.at(c).missations = {{Missation{1, rA}}, ref_sequence};
  
  assert_phylo_tree_integrity(tree);

  return tree;
}

TEST(Api_test, api_export_mcc_tree) {
  auto base_tree_1 = build_simple_tree();
  auto base_tree_2 = build_simple_tree();

  // Add a tiny bit of variation between base trees
  base_tree_1.at_root().t = -1.0;
  base_tree_2.at_root().t = -2.0;
  
  auto bitgen = std::mt19937_64{12345};
  auto mcc_tree = derive_mcc_tree({&base_tree_1, &base_tree_2}, bitgen);

  // Using traditional definition of MCC inner nodes
  EXPECT_THAT(api_export_mcc_tree(mcc_tree, false), testing::StrEq(
      absl::StrFormat("((a:%f,b:%f)NODE_0000001:%f,c:%f)NODE_0000000:%f;",
                      1.0/365.0, 2.0/365.0, (1.0+2.0)/2/365.0, (4.0+5.0)/2/365.0, 0.0)));

  // Using MRCA definition of MCC inner nodes
  EXPECT_THAT(api_export_mcc_tree(mcc_tree, true), testing::StrEq(
      absl::StrFormat("((a:%f,b:%f)NODE_0000001:%f,c:%f)NODE_0000000:%f;",
                      1.0/365.0, 2.0/365.0, (1.0+2.0)/2/365.0, (4.0+5.0)/2/365.0, 0.0)));
}

}  // namespace delphy
