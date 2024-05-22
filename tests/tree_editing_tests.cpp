#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "tree_editing.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Tree_editing_test_base : public testing::Test {
 protected:
  Real_sequence ref_sequence{rA, rC, rA, rA};
  
  Global_evo_model evo = make_global_evo_model({0, 1, 0, 1});  // 2 site partitions
  
  Tree_editing_test_base() {
    auto q_0_ab = Seq_matrix{0.0};
    q_0_ab[rA][rA] = -0.6-0.7-0.8; q_0_ab[rA][rC] = 0.6; q_0_ab[rA][rG] = 0.7; q_0_ab[rA][rT] = 0.8;
    q_0_ab[rC][rA] = 0.9; q_0_ab[rC][rC] = -0.9-1.0-1.1; q_0_ab[rC][rG] = 1.0; q_0_ab[rC][rT] = 1.1;
    q_0_ab[rG][rA] = 1.2; q_0_ab[rG][rC] = 1.3; q_0_ab[rG][rG] = -1.2-1.3-1.4; q_0_ab[rG][rT] = 1.4;
    q_0_ab[rT][rA] = 1.5; q_0_ab[rT][rC] = 1.6; q_0_ab[rT][rG] = 1.7; q_0_ab[rT][rT] = -1.5-1.6-1.7;
    EXPECT_TRUE(is_transition_rate_matrix(q_0_ab));
    
    auto q_1_ab = Seq_matrix{0.0};
    q_1_ab[rA][rA] = -2.6-2.7-2.8; q_1_ab[rA][rC] = 2.6; q_1_ab[rA][rG] = 2.7; q_1_ab[rA][rT] = 2.8;
    q_1_ab[rC][rA] = 2.9; q_1_ab[rC][rC] = -2.9-3.0-3.1; q_1_ab[rC][rG] = 3.0; q_1_ab[rC][rT] = 3.1;
    q_1_ab[rG][rA] = 3.2; q_1_ab[rG][rC] = 3.3; q_1_ab[rG][rG] = -3.2-3.3-3.4; q_1_ab[rG][rT] = 3.4;
    q_1_ab[rT][rA] = 3.5; q_1_ab[rT][rC] = 3.6; q_1_ab[rT][rG] = 3.7; q_1_ab[rT][rT] = -3.5-3.6-3.7;
    EXPECT_TRUE(is_transition_rate_matrix(q_1_ab));
    
    evo.nu_l = {0.2, 0.3, 0.4, 0.5};
    evo.partition_evo_model[0] = {
      .mu = 0.1,
      .pi_a = {0.05, 0.15, 0.25, 0.55},
      .q_ab = q_0_ab
    };
    evo.partition_evo_model[1] = {
      .mu = 1.1,
      .pi_a = {0.07, 0.17, 0.23, 0.53},
      .q_ab = q_1_ab
    };
  }

  auto mu_l(Site_index l) const -> double { return evo.mu_l(l); }
  auto nu_l(Site_index l) const -> double { return evo.nu_l[l]; }
  auto pi_l_a(Site_index l, Real_seq_letter a) const -> double { return evo.pi_l_a(l, a); }
  auto q_l_ab(Site_index l, Real_seq_letter a, Real_seq_letter b) const -> double { return evo.q_l_ab(l, a, b); }
  auto q_l_a(Site_index l, Real_seq_letter a) const -> double { return evo.q_l_a(l, a); }
};

class Tree_editing_simple_test : public Tree_editing_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                              +-- T0C -- a (CANA)
  //                              |
  //                 +A2N- A0T ---+ x (TANA)
  //                 |            |
  // (ACAA) -- C1A --+ r (AAAA)   +-------- A1G ------- b (TGNA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  // 

  Real_sequence ref_sequence{rA, rC, rA, rA};
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Phylo_tree old_tree;  // Initialized below

  Tree_editing_simple_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};

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

    old_tree = tree;
  }
};

TEST_F(Tree_editing_simple_test, slide_up_simple) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                              +-- T0C -- a (CANA)
  //                              |
  //                 +A2N- A0T ---+ x (TANA)
  //                 |            |
  // (ACAA) -- C1A --+ r (AAAA)   +-------- A1G ------- b (TGNA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  //
  // New
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                     +-------A0C-------- a (CANA)
  //                     |
  //                 +A2N+ x (AANA)
  //                 |   |
  //                 |   +-A0T------------- A1G ------- b (TGNA)
  // (ACAA) -- C1A --+ r (AAAA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  //

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(-1.0);
  edit.end();

  assert_phylo_tree_integrity(tree);
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(x).t, testing::DoubleNear(-1.0, 1e-6));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rA, 0, rC,  0.0}));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rT, -0.5}, Mutation{rA, 1, rG, 1.0}));
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(tree.at(a).missations, testing::Eq(old_tree.at(a).missations));
  EXPECT_THAT(tree.at(b).missations, testing::Eq(old_tree.at(b).missations));
  EXPECT_THAT(tree.at(c).missations, testing::Eq(old_tree.at(c).missations));
  EXPECT_THAT(tree.at(x).missations, testing::Eq(old_tree.at(x).missations));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

TEST_F(Tree_editing_simple_test, slide_up_missation_kills_mutation) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                              +T0N------ a (NANA)
  //                              |
  //                 +A2N- A0T ---+ x (TANA)
  //                 |            |
  // (ACAA) -- C1A --+ r (AAAA)   +-------- A1G ------- b (TGNA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  //
  // New
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                     +A0N--------------- a (NANA)
  //                     |
  //                 +A2N+ x (AANA)
  //                 |   |
  //                 |   +-A0T------------- A1G ------- b (TGNA)
  // (ACAA) -- C1A --+ r (AAAA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  //

  tree.at(a).mutations = {};
  tree.at(a).missations = {{Missation{0, rT}}, ref_sequence};
  assert_phylo_tree_integrity(tree);
  old_tree = tree;
  
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(-1.0);
  edit.end();

  assert_phylo_tree_integrity(tree);
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(x).t, testing::DoubleNear(-1.0, 1e-6));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rT, -0.5}, Mutation{rA, 1, rG, 1.0}));
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}));
  EXPECT_THAT(tree.at(b).missations, testing::Eq(old_tree.at(b).missations));
  EXPECT_THAT(tree.at(c).missations, testing::Eq(old_tree.at(c).missations));
  EXPECT_THAT(tree.at(x).missations, testing::Eq(old_tree.at(x).missations));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

class Tree_editing_hop_flip_test : public Tree_editing_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                                +A0N---------- a (NCNA)
  //                                |
  //                            +G2N+ x (ACNA)
  //                            |   |
  //                            |   +C1N--------- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)

  Real_sequence ref_sequence{rA, rC, rA, rT};
  Phylo_tree tree{7};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index y = 2;
  static constexpr Node_index a = 3;
  static constexpr Node_index b = 4;
  static constexpr Node_index c = 5;
  static constexpr Node_index d = 6;

  Phylo_tree old_tree;  // Initialized below

  Tree_editing_hop_flip_test() {

    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {y, d};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;

    tree.at(y).parent = r;
    tree.at(y).children = {x, c};
    tree.at(y).name = "y";
    tree.at(y).t = 0.0;
    tree.at(y).mutations = {Mutation{rA, 2, rG, -0.75}, Mutation{rT, 3, rA, -0.25}};

    tree.at(x).parent = y;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t = 0.0;
    tree.at(x).missations = {{Missation{2, rG}}, ref_sequence};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 2.0;
    tree.at(a).missations = {{Missation{0, rA}}, ref_sequence};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 3.0;
    tree.at(b).mutations = {Mutation{rA, 3, rC, 2.0}};
    tree.at(b).missations = {{Missation{1, rC}}, ref_sequence};

    tree.at(c).parent = y;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 1.0;
    tree.at(c).missations = {{Missation{1, rC}}, ref_sequence};

    tree.at(d).parent = r;
    tree.at(d).children = {};
    tree.at(d).name = "d";
    tree.at(d).t = -0.5;

    assert_phylo_tree_integrity(tree);

    old_tree = tree;
  }
};

TEST_F(Tree_editing_hop_flip_test, hop_up) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                                +A0N---------- a (NCNA)  (<--)
  //                                |
  //                            +G2N+ x (ACNA)
  //                            |   |
  //                            |   +C1N--------- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)
  //
  // New
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                            +A0N,G2N---------- a (NCNA)  (<--)
  //                            |
  //          +--- A2G -- T3A --+ x (ACGA)
  //          |                 |
  //          |                 |   +G2N--------- A3C -- b (ANNC)
  // (ACAT) --+ r (ACAT)        |   |
  //          |                 +C1N+ y (ANGA)
  //          |                     |
  //          |                     +------- c (ANGA)
  //          |
  //          +------- d (ACAT)

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.hop_up();
  edit.end();
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(d).t, testing::Eq(old_tree.at(d).t));
  EXPECT_THAT(tree.at(x).t, testing::Eq(old_tree.at(x).t));
  EXPECT_THAT(tree.at(y).t, testing::Eq(old_tree.at(y).t));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::Eq(old_tree.at(a).mutations));
  EXPECT_THAT(tree.at(b).mutations, testing::Eq(old_tree.at(b).mutations));
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(d).mutations, testing::Eq(old_tree.at(d).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::Eq(old_tree.at(y).mutations));  // Note swap
  EXPECT_THAT(tree.at(y).mutations, testing::Eq(old_tree.at(x).mutations));  // Note swap
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}, Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{2, rG}));
  EXPECT_THAT(tree.at(c).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(d).missations, testing::Eq(old_tree.at(d).missations));
  EXPECT_THAT(tree.at(x).missations, testing::IsEmpty());
  EXPECT_THAT(estd::ranges::to_vec(tree.at(y).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));

  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(tree.at(a).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(b).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(c).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(d).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(x).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(y).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(r).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.root, testing::Eq(r));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

TEST_F(Tree_editing_hop_flip_test, flip) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                                +A0N---------- a (NCNA)
  //                                |
  //                            +G2N+ x (ACNA)
  //                            |   |
  //                            |   +C1N--------- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)
  //
  // New
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                            +A0N,G2N---------- a (NCNA)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 |   +G2N---- A3C ------- b (ANNC)
  // (ACAT) --+ r (ACAT)        |   |
  //          |                 +C1N+ x (ANGA)
  //          |                     |
  //          |                     +------- c (ANGA)
  //          |
  //          +------- d (ACAT)

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, b, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.flip();
  edit.end();
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(d).t, testing::Eq(old_tree.at(d).t));
  EXPECT_THAT(tree.at(x).t, testing::Eq(old_tree.at(x).t));
  EXPECT_THAT(tree.at(y).t, testing::Eq(old_tree.at(y).t));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::Eq(old_tree.at(a).mutations));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 3, rC, 1.5}));
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(d).mutations, testing::Eq(old_tree.at(d).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::Eq(old_tree.at(x).mutations));  // No swap
  EXPECT_THAT(tree.at(y).mutations, testing::Eq(old_tree.at(y).mutations));  // No swap
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}, Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{2, rG}));
  EXPECT_THAT(tree.at(c).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(d).missations, testing::Eq(old_tree.at(d).missations));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(x).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}));
  EXPECT_THAT(tree.at(y).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));

  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(tree.at(a).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(b).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(c).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(d).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(x).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(y).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(r).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.root, testing::Eq(r));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

TEST_F(Tree_editing_hop_flip_test, hop_down) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                                +A0N---------- a (NCNA)
  //                                |
  //                            +G2N+ x (ACNA)
  //                            |   |
  //                            |   +C1N--------- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)
  //
  // New
  // ---
  //
  // Time:  -1.0               0.0 0.0   1.0      2.0   3.0
  //       
  //                            +A0N,G2N---------- a (NCNA)
  //                            |
  //          +--- A2G -- T3A --+ x (ACGA)
  //          |                 |
  //          |                 |   +G2N--------- A3C -- b (ANNC)
  // (ACAT) --+ r (ACAT)        |   |
  //          |                 +C1N+ y (ANGA)
  //          |                     |
  //          |                     +------- c (ANGA)
  //          |
  //          +------- d (ACAT)

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, c, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.hop_down(b);  // Note: hop_down(b) == do_hop_up(a) !!!
  edit.end();
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(d).t, testing::Eq(old_tree.at(d).t));
  EXPECT_THAT(tree.at(x).t, testing::Eq(old_tree.at(x).t));
  EXPECT_THAT(tree.at(y).t, testing::Eq(old_tree.at(y).t));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::Eq(old_tree.at(a).mutations));
  EXPECT_THAT(tree.at(b).mutations, testing::Eq(old_tree.at(b).mutations));
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(d).mutations, testing::Eq(old_tree.at(d).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::Eq(old_tree.at(y).mutations));  // Note swap
  EXPECT_THAT(tree.at(y).mutations, testing::Eq(old_tree.at(x).mutations));  // Note swap
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}, Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{2, rG}));
  EXPECT_THAT(tree.at(c).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(d).missations, testing::Eq(old_tree.at(d).missations));
  EXPECT_THAT(tree.at(x).missations, testing::IsEmpty());
  EXPECT_THAT(estd::ranges::to_vec(tree.at(y).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));

  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(tree.at(a).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(b).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(c).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(d).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(x).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(y).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(r).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.root, testing::Eq(r));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

class Tree_editing_slide_down_test : public Tree_editing_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:  0.0         1.0       2.0        3.0             4.0
  //
  //                     +T0N ------------------------- A1G - a (NGAA)
  //                     |
  //         +---- A0T --+ x (TAAA)
  //         |           |
  // (AAAA) -+ r (AAAA)  +-- A1C ----- T0G -- b (GCAA)
  //         |
  //         +---- c (AAAA)

  Real_sequence ref_sequence{rA, rA, rA, rA};
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Phylo_tree old_tree;  // Initialized below

  Tree_editing_slide_down_test() {

    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t = 0.0;

    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t = 1.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, 0.5}};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 4.0;
    tree.at(a).missations = {{Missation{0, rT}}, ref_sequence};
    tree.at(a).mutations = {Mutation{rA, 1, rG, 4.0}};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 3.0;
    tree.at(b).mutations = {Mutation{rA, 1, rC, 1.5}, Mutation{rT, 0, rG, 2.5}};

    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 0.5;

    assert_phylo_tree_integrity(tree);

    old_tree = tree;
  }
};

TEST_F(Tree_editing_slide_down_test, slide_down) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  // Time:  0.0         1.0       2.0        3.0             4.0
  //
  //                     +T0N ------------------------- A1G - a (NGAA)
  //                     |
  //         +---- A0T --+ x (TAAA)
  //         |           |
  // (AAAA) -+ r (AAAA)  +-- A1C ----- T0G -- b (GCAA)
  //         |
  //         +---- c (AAAA)
  //
  // New
  // ---
  // Time:  0.0         1.0       2.0        3.0             4.0
  //
  //                                          +G0N - C1G ---- a (NGAA)
  //                                          |
  //         +---- A0T ----- A1C ----- T0G ---+ x (GCAA)
  //         |                                |
  // (AAAA) -+ r (AAAA)                       +- b (GCAA)
  //         |
  //         +---- c (AAAA)
  //

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(3.0);
  edit.end();

  assert_phylo_tree_integrity(tree);
  
  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(x).t, testing::DoubleNear(3.0, 1e-6));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rC, 1, rG, 3.5}));
  EXPECT_THAT(tree.at(b).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(c).mutations, testing::Eq(old_tree.at(c).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::ElementsAre(
      Mutation{rA, 0, rT, 0.5},
      Mutation{rA, 1, rC, 1.5},
      Mutation{rT, 0, rG, 2.5}));
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rG}));
  EXPECT_THAT(tree.at(b).missations, testing::Eq(old_tree.at(b).missations));
  EXPECT_THAT(tree.at(c).missations, testing::Eq(old_tree.at(c).missations));
  EXPECT_THAT(tree.at(x).missations, testing::Eq(old_tree.at(x).missations));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

class Tree_editing_spr_test : public Tree_editing_test_base {
 protected:
  //
  // The tree that we build for a complex SPR move:
  //
  // Time:  -1.0               0.0       1.0      2.0   3.0
  //       
  //                                      +A0N---- a (NCNA)
  //                                      |
  //                            +G2N------+ x (ACNA)
  //                            |         |
  //                            |         +C1N--- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)

  Real_sequence ref_sequence{rA, rC, rA, rT};
  Phylo_tree tree{7};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index y = 2;
  static constexpr Node_index a = 3;
  static constexpr Node_index b = 4;
  static constexpr Node_index c = 5;
  static constexpr Node_index d = 6;

  Phylo_tree old_tree;  // Initialized below

  Tree_editing_spr_test() {

    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {y, d};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;

    tree.at(y).parent = r;
    tree.at(y).children = {x, c};
    tree.at(y).name = "y";
    tree.at(y).t = 0.0;
    tree.at(y).mutations = {Mutation{rA, 2, rG, -0.75}, Mutation{rT, 3, rA, -0.25}};

    tree.at(x).parent = y;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t = 1.0;
    tree.at(x).missations = {{Missation{2, rG}}, ref_sequence};

    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 2.0;
    tree.at(a).missations = {{Missation{0, rA}}, ref_sequence};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 3.0;
    tree.at(b).mutations = {Mutation{rA, 3, rC, 2.0}};
    tree.at(b).missations = {{Missation{1, rC}}, ref_sequence};

    tree.at(c).parent = y;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 1.0;
    tree.at(c).missations = {{Missation{1, rC}}, ref_sequence};

    tree.at(d).parent = r;
    tree.at(d).children = {};
    tree.at(d).name = "d";
    tree.at(d).t = -0.5;

    assert_phylo_tree_integrity(tree);

    old_tree = tree;
  }
};

TEST_F(Tree_editing_spr_test, spr) {
  auto scratch = Scratch_space{};
  
  // Old
  // ---
  //
  //
  // Time:  -1.0               0.0       1.0      2.0   3.0
  //       
  //                                      +A0N---- a (NCNA)
  //                                      |
  //                            +G2N------+ x (ACNA)
  //                            |         |
  //                            |         +C1N--- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 +C1N---- c (ANGA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------- d (ACAT)
  //
  // New
  // ---
  //
  // Time:  -1.0               0.0       1.0      2.0   3.0
  //       
  //                            +C1N,G2N--------- A3C -- b (ANNC)
  //                            |
  //          +--- A2G -- T3A --+ y (ACGA)
  //          |                 |
  //          |                 |   +C1N c (ANGA)
  //          |                 |   |
  //          |                 +---+ x (ACGA)
  //          |                     |
  //          |                     +A0N,G2N------ a (NCNA)
  // (ACAT) --+ r (ACAT)
  //          |
  //          +------ d (ACAT)

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  // Slide x all the way up to just below r
  edit.slide_P_along_branch(tree.at(y).t);
  edit.hop_up();
  edit.slide_P_along_branch(tree.at(r).t);
  
  // Slide down to start of y-b
  edit.slide_P_along_branch(tree.at(y).t);
  edit.hop_down(b);

  // Switch to y-c
  edit.flip();

  // Slide down to middle of the y-c branch
  edit.slide_P_along_branch(0.5);
  edit.end();
  
  assert_phylo_tree_integrity(tree);

  EXPECT_THAT(tree.at(a).t, testing::Eq(old_tree.at(a).t));
  EXPECT_THAT(tree.at(b).t, testing::Eq(old_tree.at(b).t));
  EXPECT_THAT(tree.at(c).t, testing::Eq(old_tree.at(c).t));
  EXPECT_THAT(tree.at(d).t, testing::Eq(old_tree.at(d).t));
  EXPECT_THAT(tree.at(x).t, testing::DoubleNear(0.5, 1e-6));
  EXPECT_THAT(tree.at(y).t, testing::Eq(old_tree.at(y).t));
  EXPECT_THAT(tree.at(r).t, testing::Eq(old_tree.at(r).t));
  
  EXPECT_THAT(tree.at(a).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 3, rC, 2.0}));
  EXPECT_THAT(tree.at(c).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(d).mutations, testing::Eq(old_tree.at(d).mutations));
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(y).mutations, testing::Eq(old_tree.at(y).mutations));
  EXPECT_THAT(tree.at(r).mutations, testing::Eq(old_tree.at(r).mutations));
  
  EXPECT_THAT(estd::ranges::to_vec(tree.at(a).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{0, rA}, Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(b).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}, Missation{2, rG}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(c).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}));
  EXPECT_THAT(tree.at(d).missations, testing::Eq(old_tree.at(d).missations));
  EXPECT_THAT(tree.at(x).missations, testing::IsEmpty());
  EXPECT_THAT(tree.at(y).missations, testing::Eq(old_tree.at(y).missations));
  EXPECT_THAT(tree.at(r).missations, testing::Eq(old_tree.at(r).missations));

  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(tree.at(a).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(b).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(c).parent, testing::Eq(x));
  EXPECT_THAT(tree.at(d).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(x).parent, testing::Eq(y));
  EXPECT_THAT(tree.at(y).parent, testing::Eq(r));
  EXPECT_THAT(tree.at(r).parent, testing::Eq(k_no_node));
  EXPECT_THAT(tree.root, testing::Eq(r));
  
  // Check lambda_i and num_sites_missing_at_every_node are updated
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

TEST_F(Tree_editing_simple_test, slide_up_against_edge_mutation) {
  auto scratch = Scratch_space{};

  // If we slide x right up against the beginning of its parent branch, all the mutations
  // before it should migrate to the child branch (edge case for zero-length branches).
  
  // Old
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                              +-- T0C -- a (CANA)
  //                              |
  //                 +A2N,A0T ----+ x (TANA)
  //                 |            |
  // (ACAA) -- C1A --+ r (AAAA)   +-------- A1G ------- b (TGNA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)
  //
  // New
  // ---
  // Time:         -1.0          0.0        1.0        2.0        3.0
  //          
  //                     +-------A0C-------- a (CANA)
  //                     |
  //                 +A2N+ x (AANA)
  //                 |   |
  //                 |   +A0T-------------- A1G ------- b (TGNA)
  // (ACAA) -- C1A --+ r (AAAA)
  //                 |
  //                 +A1N------------------ A0G ------------------ c (GNAA)

  tree.at(x).mutations = {Mutation{rA, 0, rT, tree.at(r).t}};

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(tree.at(r).t);
  edit.end();

  assert_phylo_tree_integrity(tree);

  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rA, 0, rC,  0.0}));
  EXPECT_THAT(tree.at(b).mutations, testing::ElementsAre(Mutation{rA, 0, rT, -1.0}, Mutation{rA, 1, rG, 1.0}));
  EXPECT_THAT(tree.at(x).mutations, testing::IsEmpty());
}

TEST_F(Tree_editing_slide_down_test, slide_down_against_edge_mutation) {
  auto scratch = Scratch_space{};

  // If we slide x right down against the end of its parent branch, all the mutations
  // after it should migrate to the parent branch (edge case for zero-length branches).
  
  // Old
  // ---
  // Time:  0.0         1.0       2.0        3.0             4.0
  //
  //                     +T0N ------------------------- A1G - a (NGAA)
  //                     |
  //         +---- A0T --+ x (TAAA)
  //         |           |
  // (AAAA) -+ r (AAAA)  +-- A1C -------- T0G b (GCAA)
  //         |
  //         +---- c (AAAA)
  //
  // New
  // ---
  // Time:  0.0         1.0       2.0        3.0             4.0
  //
  //                                          +G0N - C1G ---- a (NGAA)
  //                                          |
  //         +---- A0T ----- A1C -------- T0G + x (GCAA)
  //         |                                |
  // (AAAA) -+ r (AAAA)                       +- b (GCAA)
  //         |
  //         +---- c (AAAA)
  //

  tree.at(b).mutations = {Mutation{rA, 1, rC, 1.5}, Mutation{rT, 0, rG, tree.at(b).t}};

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, a, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(tree.at(b).t);
  edit.end();

  assert_phylo_tree_integrity(tree);

  EXPECT_THAT(tree.at(a).mutations, testing::ElementsAre(Mutation{rC, 1, rG,  3.5}));
  EXPECT_THAT(tree.at(b).mutations, testing::IsEmpty());
  EXPECT_THAT(tree.at(x).mutations, testing::ElementsAre(
      Mutation{rA, 0, rT, 0.5},
      Mutation{rA, 1, rC, 1.5},
      Mutation{rT, 0, rG, tree.at(b).t}));
}

TEST_F(Tree_editing_slide_down_test, slide_down_root) {
  auto scratch = Scratch_space{};

  // As always, handling the root is very tricky, because mutations on the sibling branch
  // don't simply migrate intact to above the root...
  
  // Old
  // ---
  // Time:         0.0         1.0       2.0        3.0             4.0
  //
  //                            +T0N ------------------------- A1G - a (NGAA)
  //                            |
  //                +--- C1A ---+ x (TAAA)
  //                |           |
  // (AAAA) -A0T,A1C+ r (TCAA)  +-- A1C -------- T0G b (GCAA)
  //                |
  //                +C1A,A2G-- c (TAGA)
  //
  // New
  // ---
  // Time:         0.0         1.0       2.0        3.0             4.0
  //
  //                            +T0N ------------------------- A1G - a (NGAA)
  //                            |
  //                  +-- G2A --+ x (TAAA)
  //                  |         |
  // (AAAA)   -A0T,A2G+ r(TAGA) +-- A1C -------- T0G b (GCAA)
  //                  |
  //                  +------- c (TAGA)
  //

  tree.at(c).t = 0.9;
  tree.at(x).mutations = {Mutation{rC, 1, rA, 0.5}};
  tree.at(r).mutations = {
    Mutation{rA, 0, rT, -std::numeric_limits<double>::max()},
    Mutation{rA, 1, rC, -std::numeric_limits<double>::max()}};
  tree.at(c).mutations = {Mutation{rC, 1, rA, 0.02}, Mutation{rA, 2, rG, 0.05}};

  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  
  auto edit = Tree_editing_session{tree, x, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  edit.slide_P_along_branch(0.5);
  edit.end();

  assert_phylo_tree_integrity(tree);

  EXPECT_THAT(tree.at(x).mutations, testing::ElementsAre(Mutation{rG, 2, rA, 0.75}));
  EXPECT_THAT(tree.at(r).mutations, testing::UnorderedElementsAre(
    Mutation{rA, 0, rT, -std::numeric_limits<double>::max()},
    Mutation{rA, 2, rG, -std::numeric_limits<double>::max()}));
  EXPECT_THAT(tree.at(c).mutations, testing::IsEmpty());
}

}  // namespace delphy
