#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <absl/random/mocking_bit_gen.h>
#include <absl/random/mock_distributions.h>

#include "phylo_tree.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Phylo_tree_test, push_all_mutations_and_missations_to_tips_trivial) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rA, rA, rA, rA};
  tree.root = 0;

  push_all_mutations_and_missations_to_tips(tree);

  EXPECT_THAT(tree.at(0).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(0).missations, testing::IsEmpty());
}

TEST(Phylo_tree_test, push_all_mutations_and_missations_to_tips_accept_muts_above_root) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rC, rA, rA, rA};
  tree.root = 0;
  tree.at(0).mutations = {Mutation{rC, 0, rA, -1.0}};

  push_all_mutations_and_missations_to_tips(tree);  // Should not throw
}

TEST(Phylo_tree_test, push_all_mutations_and_missations_to_tips_simple) {
  //
  // The tree that we build:
  //
  // Time:        0.0        1.0        2.0        3.0
  //         
  //               +-- T0C -- a (CT)
  //               |
  // (AT) -- A0T --+ r (TT)
  //               |
  //               +T1N----- T0G ------- b (GN)
  // 


  auto ref_sequence = Real_sequence{rA, rT};
  auto tree = Phylo_tree{3};

  const auto r = Node_index{0};
  const auto a = Node_index{1};
  const auto b = Node_index{2};

  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {a, b};
  tree.at(r).t = 0.0;
  tree.at(r).mutations = {Mutation{rA, 0, rT, -0.5}};

  tree.at(a).parent = r;
  tree.at(a).children = {};
  tree.at(a).t = 1.0;
  tree.at(a).mutations = {Mutation{rT, 0, rC, 0.5}};

  tree.at(b).parent = r;
  tree.at(b).children = {};
  tree.at(b).t = 2.0;
  tree.at(b).mutations = {Mutation{rT, 0, rG, 1.0}};
  tree.at(b).missations = {{Missation{1, rT}}, ref_sequence};

  push_all_mutations_and_missations_to_tips(tree);

  EXPECT_THAT(tree.at(r).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(r).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rA, 0, rC, 1.0}));
  EXPECT_THAT(tree.at(a).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rG, 2.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rT}));
}

class Phylo_tree_complex_test : public testing::Test {
 protected:
  //
  // The tree that we build:
  //
  // Time:            -1.0         0.0        1.0        2.0        3.0
  //                 
  //                                +-- T0C -- a (CANN)
  //                                |
  //                    +A2N A0T ---+ x (TANN)
  //                    |           |
  //                    |           +-------- A1G ------- b (TGNN)
  // (AACA) A3N-- C2A --+ r (AAAN)
  //                    |
  //                    +A1N----------------- A0G ------------------ c (GNAN)
  // 

  Real_sequence ref_sequence{rA, rA, rC, rA};
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Phylo_tree_complex_test() {
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
  }
};

TEST_F(Phylo_tree_complex_test, push_all_mutations_and_missations_to_tips_bigger) {
  push_all_mutations_and_missations_to_tips(tree);

  EXPECT_THAT(tree.at(r).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(r).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(x).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rA, 0, rC, 1.0},
                                                         Mutation{rC, 2, rA, 1.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{2, rC}, Missation{3, rA}));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rT, 2.0},
                                                         Mutation{rA, 1, rG, 2.0},
                                                         Mutation{rC, 2, rA, 2.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{2, rC}, Missation{3, rA}));
  EXPECT_THAT(tree.at(c).mutations, testing::ElementsAre(Mutation{rA, 0, rG, 3.0},
                                                         Mutation{rC, 2, rA, 3.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(c).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rA}, Missation{3, rA}));
}

TEST(Phylo_tree_test, extract_nodes_trivial) {
  auto tree = Phylo_tree{};
  auto [tips, inner_nodes] = extract_nodes(tree);

  EXPECT_THAT(tips, testing::IsEmpty());
  EXPECT_THAT(inner_nodes, testing::IsEmpty());
}

TEST_F(Phylo_tree_complex_test, extract_nodes_typical) {
  auto [tips, inner_nodes] = extract_nodes(tree);

  EXPECT_THAT(tips, testing::UnorderedElementsAre(a, b, c));
  EXPECT_THAT(inner_nodes, testing::UnorderedElementsAre(r, x));
}

TEST(Phylo_tree_test, rewire_tree_through_sequential_accretion_trivial) {
  auto tree = Phylo_tree{};
  auto tips = std::vector<Node_index>{};

  rewire_tree_through_sequential_accretion(tree, tips);

  EXPECT_THAT(std::ssize(tree), testing::Eq(0));
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_wrong_number_of_tip_indices) {
  auto tips = std::vector<Node_index>{};

  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips)), std::invalid_argument);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_duplicate_tip_indices) {
  auto tips = std::vector<Node_index>{a, a, a};

  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips)), std::invalid_argument);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_non_tips_in_tip_indices) {
  auto tips = std::vector<Node_index>{a, r, x};

  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips)), std::invalid_argument);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_out_of_range_indices) {
  auto tips_1 = std::vector<Node_index>{-1, -1, -1};
  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips_1)), std::out_of_range);
  
  auto tips_2 = std::vector<Node_index>{9, 9, 9};
  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips_2)), std::out_of_range);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_mutations_on_inner_nodes) {
  // The tree in the test fixture already has mutations on x !  Ensure missations aren't causing exceptions
  for (const auto& node : index_order_traversal(tree)) {
    tree.at(node).missations = {};
  }
  
  auto tips = std::vector<Node_index>{a, b, c};
  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips)), std::invalid_argument);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_invalid_missations_on_inner_nodes) {
  // The tree in the test fixture already has missations on x !  Ensure mutations aren't causing exceptions
  for (const auto& node : index_order_traversal(tree)) {
    tree.at(node).mutations = {};
  }
  
  auto tips = std::vector<Node_index>{a, b, c};
  EXPECT_THROW((rewire_tree_through_sequential_accretion(tree, tips)), std::invalid_argument);
}

TEST_F(Phylo_tree_complex_test, rewire_tree_through_sequential_accretion_typical) {
  push_all_mutations_and_missations_to_tips(tree);

  rewire_tree_through_sequential_accretion(tree, {c, b, a}, 10.0);

  EXPECT_THAT(tree.at(c).parent, testing::Eq(tree.at(b).parent));
  auto y = tree.at(c).parent;
  EXPECT_THAT(tree.at(y).parent, testing::Eq(tree.at(a).parent));
  EXPECT_THAT(tree.at(y).parent, testing::Eq(tree.root));
  auto z = tree.at(y).parent;

  EXPECT_THAT(tree.at(a).t, testing::DoubleNear(1.0, 1e-6));
  EXPECT_THAT(tree.at(b).t, testing::DoubleNear(2.0, 1e-6));
  EXPECT_THAT(tree.at(c).t, testing::DoubleNear(3.0, 1e-6));
  EXPECT_THAT(tree.at(y).t, testing::DoubleNear(std::min(2.0, 3.0) - 10.0, 1e-6));
  EXPECT_THAT(tree.at(z).t, testing::DoubleNear(std::min(tree.at(y).t, 1.0) - 10.0, 1e-6));

  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rA, 0, rC, 1.0},
                                                         Mutation{rC, 2, rA, 1.0}));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rT, 2.0},
                                                         Mutation{rA, 1, rG, 2.0},
                                                         Mutation{rC, 2, rA, 2.0}));
  EXPECT_THAT(tree.at(c).mutations, testing::ElementsAre(Mutation{rA, 0, rG, 3.0},
                                                         Mutation{rC, 2, rA, 3.0}));
  EXPECT_THAT(tree.at(y).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(z).mutations, testing::IsEmpty());
}

TEST_F(Phylo_tree_complex_test, assert_mutation_consistency_ok) {
  assert_mutation_consistency(tree);  // should be fine
}

#ifndef __EMSCRIPTEN__  // We want tests to be compiled to exercise clang, but EXPECT_DEBUG_DEATH doesn't exist in Emscripten
TEST_F(Phylo_tree_complex_test, assert_mutation_consistency_different_from_root) {
  // Mutation from != to
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rA, 0, rA, -0.5}};
        assert_mutation_consistency(tree);
      }), "A vs. A");

  // Mutation with site < 0
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rA, -15, rT, -0.5}};
        assert_mutation_consistency(tree);
      }), "0 <= m.site");

  // Mutation with site >= genome size
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rA, 15, rT, -0.5}};
        assert_mutation_consistency(tree);
      }), "15 vs. 4");

  // Mutation with time before parent node time
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rA, 0, rT, -1.5}};
        assert_mutation_consistency(tree);
      }), "-1 vs. -1.5");
  
  // Mutation with time after node time
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rA, 0, rT, +1.5}};
        assert_mutation_consistency(tree);
      }), "1.5 vs. 0");

  // Mutation with from not matching current seq
  EXPECT_DEBUG_DEATH(({
        tree.at(x).mutations = {Mutation{rC, 0, rT, -0.5}};  // root is A at site 0!
        assert_mutation_consistency(tree);
      }), "C vs. A");
}
#endif

TEST_F(Phylo_tree_complex_test, assert_missation_consistency_ok) {
  assert_missation_consistency(tree);
}

#ifndef __EMSCRIPTEN__  // We want tests to be compiled to exercise clang, but EXPECT_DEBUG_DEATH doesn't exist in Emscripten
TEST_F(Phylo_tree_complex_test, assert_missation_consistency_bad) {
  // Missation below another missation
  EXPECT_DEBUG_DEATH(({
        tree.at(a).missations = {{Missation{2, rA}}, ref_sequence};  // conflicts with A2N on upstream r->x branch
        assert_missation_consistency(tree);
      }), "not.*cur_missing.*contains");

  // Mutation on a site that's marked missing on the same branch
  EXPECT_DEBUG_DEATH(({
        tree.at(x).missations = {{Missation{0, rA}}, ref_sequence};  // conflicts with A0T mutation on r->x branch
        assert_missation_consistency(tree);
      }), "not.*cur_missing.*contains");

  // Mutation on a site that's marked missing on an ancestor branch
  EXPECT_DEBUG_DEATH(({
        tree.at(x).missations = {{Missation{1, rA}}, ref_sequence};  // conflicts with A1G on x->b branch
        tree.at(c).missations = {};              // ensure A1N doesn't appear on both branches out of r
        assert_missation_consistency(tree);
      }), "not.*cur_missing.*contains");

  // Missations for the same site on all branches leading out of a node
  EXPECT_DEBUG_DEATH(({
        tree.at(x).missations = {{Missation{2, rA}}, ref_sequence};
        tree.at(c).missations = {{Missation{2, rA}}, ref_sequence};  // missation should bubble up to above root
        assert_missation_consistency(tree);
      }), "common_sites.empty()");
}
#endif

TEST_F(Phylo_tree_complex_test, find_MRCA_of_nodes) {
  EXPECT_THAT(find_MRCA_of(tree, a, a), testing::Eq(a));
  EXPECT_THAT(find_MRCA_of(tree, a, b), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, a, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, a, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, a, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, a, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, b, a), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, b, b), testing::Eq(b));
  EXPECT_THAT(find_MRCA_of(tree, b, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, b, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, b, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, b, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, c, a), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, b), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, c), testing::Eq(c));
  EXPECT_THAT(find_MRCA_of(tree, c, x), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, x, a), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, b), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, x, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, x, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, r, a), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, b), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, x), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, a), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, b), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, c), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, x), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, r), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, k_no_node), testing::Eq(k_no_node));
}

TEST_F(Phylo_tree_complex_test, find_MRCA_of_tree_locs) {
  EXPECT_THAT(find_MRCA_of(tree, {a, 1.0}, {a,  1.0}), testing::Eq(tree.node_loc(a)));
  EXPECT_THAT(find_MRCA_of(tree, {a, 1.0}, {a,  0.5}), testing::Eq(Phylo_tree_loc{a, 0.5}));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {b,  1.0}), testing::Eq(tree.node_loc(x)));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {c,  1.0}), testing::Eq(tree.node_loc(r)));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {x,  0.0}), testing::Eq(tree.node_loc(x)));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {x, -0.5}), testing::Eq(Phylo_tree_loc{x, -0.5}));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {r, -1.0}), testing::Eq(tree.node_loc(r)));
  EXPECT_THAT(find_MRCA_of(tree, {a, 0.5}, {r, -1.5}), testing::Eq(Phylo_tree_loc{r, -1.5}));
}

TEST_F(Phylo_tree_complex_test, find_MRCA_of_nodes_all_times_equal) {
  for (const auto& node : index_order_traversal(tree)) {
    tree.at(node).t = 0.0;
  }
  
  EXPECT_THAT(find_MRCA_of(tree, a, a), testing::Eq(a));
  EXPECT_THAT(find_MRCA_of(tree, a, b), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, a, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, a, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, a, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, a, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, b, a), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, b, b), testing::Eq(b));
  EXPECT_THAT(find_MRCA_of(tree, b, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, b, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, b, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, b, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, c, a), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, b), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, c), testing::Eq(c));
  EXPECT_THAT(find_MRCA_of(tree, c, x), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, c, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, x, a), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, b), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, x, x), testing::Eq(x));
  EXPECT_THAT(find_MRCA_of(tree, x, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, x, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, r, a), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, b), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, c), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, x), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, r), testing::Eq(r));
  EXPECT_THAT(find_MRCA_of(tree, r, k_no_node), testing::Eq(k_no_node));
  
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, a), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, b), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, c), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, x), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, r), testing::Eq(k_no_node));
  EXPECT_THAT(find_MRCA_of(tree, k_no_node, k_no_node), testing::Eq(k_no_node));
}

TEST_F(Phylo_tree_complex_test, descends_from_nodes) {
  EXPECT_THAT(descends_from(tree, a, a), true);
  EXPECT_THAT(descends_from(tree, a, b), false);
  EXPECT_THAT(descends_from(tree, a, c), false);
  EXPECT_THAT(descends_from(tree, a, x), true);
  EXPECT_THAT(descends_from(tree, a, r), true);
  EXPECT_THAT(descends_from(tree, a, k_no_node), true);
  
  EXPECT_THAT(descends_from(tree, b, a), false);
  EXPECT_THAT(descends_from(tree, b, b), true);
  EXPECT_THAT(descends_from(tree, b, c), false);
  EXPECT_THAT(descends_from(tree, b, x), true);
  EXPECT_THAT(descends_from(tree, b, r), true);
  EXPECT_THAT(descends_from(tree, b, k_no_node), true);
  
  EXPECT_THAT(descends_from(tree, c, a), false);
  EXPECT_THAT(descends_from(tree, c, b), false);
  EXPECT_THAT(descends_from(tree, c, c), true);
  EXPECT_THAT(descends_from(tree, c, x), false);
  EXPECT_THAT(descends_from(tree, c, r), true);
  EXPECT_THAT(descends_from(tree, c, k_no_node), true);
  
  EXPECT_THAT(descends_from(tree, x, a), false);
  EXPECT_THAT(descends_from(tree, x, b), false);
  EXPECT_THAT(descends_from(tree, x, c), false);
  EXPECT_THAT(descends_from(tree, x, x), true);
  EXPECT_THAT(descends_from(tree, x, r), true);
  EXPECT_THAT(descends_from(tree, x, k_no_node), true);
  
  EXPECT_THAT(descends_from(tree, r, a), false);
  EXPECT_THAT(descends_from(tree, r, b), false);
  EXPECT_THAT(descends_from(tree, r, c), false);
  EXPECT_THAT(descends_from(tree, r, x), false);
  EXPECT_THAT(descends_from(tree, r, r), true);
  EXPECT_THAT(descends_from(tree, r, k_no_node), true);
  
  EXPECT_THAT(descends_from(tree, k_no_node, a), false);
  EXPECT_THAT(descends_from(tree, k_no_node, b), false);
  EXPECT_THAT(descends_from(tree, k_no_node, c), false);
  EXPECT_THAT(descends_from(tree, k_no_node, x), false);
  EXPECT_THAT(descends_from(tree, k_no_node, r), false);
  EXPECT_THAT(descends_from(tree, k_no_node, k_no_node), true);
}

TEST_F(Phylo_tree_complex_test, descends_from_tree_locs) {
  EXPECT_THAT(descends_from(tree, {a, 1.0}, {a,  1.0}), true);
  EXPECT_THAT(descends_from(tree, {a, 1.0}, {a,  0.5}), true);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {a,  1.0}), false);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {b,  1.0}), false);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {c,  1.0}), false);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {x,  0.0}), true);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {x, -0.5}), true);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {r, -1.0}), true);
  EXPECT_THAT(descends_from(tree, {a, 0.5}, {r, -1.5}), true);
}

TEST_F(Phylo_tree_complex_test, multiple_mutations_on_branch_are_ok) {
  //
  // Time:            -1.0             0.0          1.0        2.0        3.0
  //
  //   before           +----- A0T -----+ x (TANN)
  //
  //   after            +- A0C --- C0T -+ x (TANN)
  //
  tree.at(x).mutations = {Mutation{rA, 0, rC, -0.75}, Mutation{rC, 0, rT, -0.25}};
  assert_phylo_tree_integrity(tree);
}

TEST(Phylo_tree_test, fix_up_missations_trivial) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rA, rA, rA, rA};
  tree.root = 0;

  fix_up_missations(tree);

  EXPECT_THAT(tree.at(0).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(0).missations, testing::IsEmpty());

  assert_missation_consistency(tree);
}

TEST(Phylo_tree_test, fix_up_missation_redundant_downstream_missations) {
  //
  // The tree that we build:
  //
  // Time:   0.0        1.0        2.0        3.0
  //         
  //          +--------- a (N)
  //          |
  // (T) T0N--+ r (N)
  //          |
  //          +T0N----------------- b (N)  <-- this missation is redudant
  // 

  auto ref_sequence = Real_sequence{rT};
  auto tree = Phylo_tree{3};

  const auto r = Node_index{0};
  const auto a = Node_index{1};
  const auto b = Node_index{2};

  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {a, b};
  tree.at(r).t = 0.0;
  tree.at(r).missations = {{Missation{0, rT}}, ref_sequence};

  tree.at(a).parent = r;
  tree.at(a).children = {};
  tree.at(a).t = 1.0;

  tree.at(b).parent = r;
  tree.at(b).children = {};
  tree.at(b).t = 2.0;
  tree.at(b).missations = {{Missation{0, rT}}, ref_sequence};

  fix_up_missations(tree);

  EXPECT_THAT(estd::ranges::to_vec(tree.at(r).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rT}));
  EXPECT_THAT(tree.at(a).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).missations, testing::IsEmpty());

  assert_missation_consistency(tree);
}

TEST(Phylo_tree_test, fix_up_missations_mutations_below_missations) {
  //
  // The tree that we build:
  //
  // Time:      0.0        1.0        2.0        3.0
  //         
  //             +-- A0T -- a (NA)  <-- this mutation can be dropped
  //             |
  //  (AA) A0N---+ r (NA)
  //             |
  //             +-------- A1C ------- b (NC)
  // 

  auto ref_sequence = Real_sequence{rA, rA};
  auto tree = Phylo_tree{3};

  const auto r = Node_index{0};
  const auto a = Node_index{1};
  const auto b = Node_index{2};

  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {a, b};
  tree.at(r).t = 0.0;
  tree.at(r).missations = {{Missation{0, rA}}, ref_sequence};

  tree.at(a).parent = r;
  tree.at(a).children = {};
  tree.at(a).t = 1.0;
  tree.at(a).mutations = {Mutation{rA, 0, rT, 0.5}};

  tree.at(b).parent = r;
  tree.at(b).children = {};
  tree.at(b).t = 2.0;
  tree.at(b).mutations = {Mutation{rA, 1, rC, 1.0}};

  fix_up_missations(tree);

  EXPECT_THAT(tree.at(r).mutations, testing::IsEmpty());
  EXPECT_THAT(estd::ranges::to_vec(tree.at(r).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}));
  EXPECT_THAT(tree.at(a).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(a).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 1, rC, 1.0}));
  EXPECT_THAT(tree.at(b).missations, testing::IsEmpty());

  assert_missation_consistency(tree);
}

TEST(Phylo_tree_test, fix_up_missations_merge_up_once) {
  //
  // The tree that we build:
  //
  // Time:   0.0        1.0        2.0        3.0
  //         
  //          +A0N------ a (N)            <-- this missation...
  //          |
  // (A) -----+ r (A)
  //          |
  //          +A0N----------------- b (N) <-- ...and this one can be merged into a missation above r
  // 

  auto ref_sequence = Real_sequence{rA};
  auto tree = Phylo_tree{3};

  const auto r = Node_index{0};
  const auto a = Node_index{1};
  const auto b = Node_index{2};

  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {a, b};
  tree.at(r).t = 0.0;

  tree.at(a).parent = r;
  tree.at(a).children = {};
  tree.at(a).t = 1.0;
  tree.at(a).missations = {{Missation{0, rA}}, ref_sequence};

  tree.at(b).parent = r;
  tree.at(b).children = {};
  tree.at(b).t = 2.0;
  tree.at(b).missations = {{Missation{0, rA}}, ref_sequence};

  fix_up_missations(tree);

  EXPECT_THAT(estd::ranges::to_vec(tree.at(r).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}));
  EXPECT_THAT(tree.at(a).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).missations, testing::IsEmpty());

  assert_missation_consistency(tree);
}

TEST_F(Phylo_tree_complex_test, fix_up_missations_merge_up_twice_and_drop_upstream_mutations) {
  //
  // The tree that we build:
  //
  // Time:            -1.0         0.0           1.0        2.0        3.0
  //                 
  //                                +T0N-- T0C -- a (NANN)  <-- this mutation dropped; this missation and 
  //                                |                           the sibling one bubble up...
  //                    +A2N A0T ---+ x (TANN)              <-- ...and force us to drop this mutations
  //                    |           |
  //                    |           +T0N-------- A1G ------- b (NGNN)
  // (AACA) A3N-- C2A --+ r (AAAN)
  //                    |
  //                    +A1N-------------------- A0G ------------------ c (GNAN)

  tree.at(a).missations = {{Missation{0, rT}}, ref_sequence};
  tree.at(b).missations = {{Missation{0, rT}}, ref_sequence};

  fix_up_missations(tree);

  EXPECT_THAT(tree.at(a).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(a).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 1, rG, 1.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(x).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}, Missation{2, rA}));
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
}

TEST_F(Phylo_tree_complex_test, fix_up_missations_two_mutations_on_same_logical_branch_not_changed) {
  //
  // The tree that we build:
  //
  // Time:            -1.0         0.0           1.0        2.0        3.0
  //                 
  //                                +A1N-- A0C -- a (CNNN)
  //                                |
  //                    +A2N A1G ---+ x (AGNN)
  //                    |           |
  //                    |           +----------- G1C ------- b (ACNN)
  // (AACA) A3N-- C2A --+ r (AAAN)
  //                    |
  //                    +A1N-------------------- A0G ------------------ c (GNAN)
  //
  // At some point, we tried enforcing the invariant that logical branches (single branches on the site-specific trees)
  // had at most 1 mutation.  This idea ended up being unworkable (the invariant cannot be maintained when a logical
  // branch is split over two different subruns).  This test preserved to make sure that fix_up_missations never merges 2+
  // mutations on a logical branch.

  tree.at(a).missations = {{Missation{1, rA}}, ref_sequence};  // - -> r -> x -> b becomes a single logical branch from t=-1.0 to t=2.0
  tree.at(a).mutations = {Mutation{rA, 0, rC, +0.5}};
  tree.at(x).mutations = {Mutation{rA, 1, rG, -0.5}};  // should not be merged into a single A1C mutation in the middle
  tree.at(b).mutations = {Mutation{rG, 1, rC, 1.0}};   // of the above logical branch (t = 0.5)

  fix_up_missations(tree);
  assert_phylo_tree_integrity(tree);

  EXPECT_THAT(tree.at(x).mutations, testing::ElementsAre(Mutation{rA, 1, rG, -0.5}));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rG, 1, rC, 1.0}));
}

TEST_F(Phylo_tree_complex_test, randomize_tree) {
  // This test isn't guaranteed to fail if there's something wrong, but we can hope...
  auto bitgen = absl::BitGen{};
  randomize_tree(tree, bitgen);
}

}  // namespace delphy
