#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "phylo_tree_calc.h"
#include "seq_vector_matchers.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Phylo_tree_calc_complex_test : public testing::Test {
 protected:
  //
  // The tree that we build:
  //
  // Time:             -1.0          0.0        1.0        2.0        3.0
  //                   
  //                                  +-- T0C -- a (CANN)
  //                                  |
  //                     +A2N- A0T ---+ x (TANN)
  //                     |            |
  // (AACA) A3N- C2A --+ r (AAAN)     +-------- A1G ------- b (TGNN)
  //                     |
  //                     +A1N--------A0T------- T0G ------------------ c (GNAN)
  // 

  Real_sequence ref_sequence{rA, rA, rC, rA};
  Phylo_tree tree{5};

  Global_evo_model evo = make_global_evo_model({0, 1, 0, 1});  // 2 site partitions

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  auto mu_l(Site_index l) const -> double { return evo.mu_l(l); }
  auto nu_l(Site_index l) const -> double { return evo.nu_l[l]; }
  auto pi_l_a(Site_index l, Real_seq_letter a) const -> double { return evo.pi_l_a(l, a); }
  auto q_l_ab(Site_index l, Real_seq_letter a, Real_seq_letter b) const -> double { return evo.q_l_ab(l, a, b); }
  auto q_l_a(Site_index l, Real_seq_letter a) const -> double { return evo.q_l_a(l, a); }
  
  Phylo_tree_calc_complex_test() {
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
    tree.at(c).mutations = {Mutation{rA, 0, rT, 0.0}, Mutation{rT, 0, rG, 1.0}};
    tree.at(c).missations = {{Missation{1, rA}}, ref_sequence};

    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Phylo_tree_calc_complex_test, view_of_sequence_at_loc) {
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{r, -1.5}), testing::ElementsAre(rA, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{x, -0.8}), testing::ElementsAre(rA, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{x,  0.0}), testing::ElementsAre(rT, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{a,  0.2}), testing::ElementsAre(rT, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{a,  1.0}), testing::ElementsAre(rC, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{b,  0.5}), testing::ElementsAre(rT, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{b,  2.0}), testing::ElementsAre(rT, rG, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{c, -0.5}), testing::ElementsAre(rA, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{c,  0.5}), testing::ElementsAre(rT, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, Phylo_tree_loc{c,  1.5}), testing::ElementsAre(rG, rA, rA, rA));
}

TEST_F(Phylo_tree_calc_complex_test, view_of_sequence_at_node) {
  EXPECT_THAT(view_of_sequence_at(tree, r), testing::ElementsAre(rA, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, x), testing::ElementsAre(rT, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, a), testing::ElementsAre(rC, rA, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, b), testing::ElementsAre(rT, rG, rA, rA));
  EXPECT_THAT(view_of_sequence_at(tree, c), testing::ElementsAre(rG, rA, rA, rA));
}

TEST_F(Phylo_tree_calc_complex_test, view_of_sequence_at_with_reversion) {
  // This used to crash!
  tree.at(a).mutations = {Mutation{rT, 0, rA, 0.5}};  // r -- A0T -- x -- T0A -- a
  EXPECT_THAT(view_of_sequence_at(tree, a), testing::ElementsAre(rA, rA, rA, rA));
}

TEST_F(Phylo_tree_calc_complex_test, reconstruct_missing_sites_at) {
  EXPECT_THAT(estd::ranges::to_vec(reconstruct_missing_sites_at(tree, r).slow_elements()), testing::ElementsAre(3));
  EXPECT_THAT(estd::ranges::to_vec(reconstruct_missing_sites_at(tree, x).slow_elements()), testing::ElementsAre(2, 3));
  EXPECT_THAT(estd::ranges::to_vec(reconstruct_missing_sites_at(tree, a).slow_elements()), testing::ElementsAre(2, 3));
  EXPECT_THAT(estd::ranges::to_vec(reconstruct_missing_sites_at(tree, b).slow_elements()), testing::ElementsAre(2, 3));
  EXPECT_THAT(estd::ranges::to_vec(reconstruct_missing_sites_at(tree, c).slow_elements()), testing::ElementsAre(1, 3));
}

TEST_F(Phylo_tree_calc_complex_test, is_site_missing_at) {
  EXPECT_THAT(is_site_missing_at(tree, r, 0), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, r, 1), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, r, 2), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, r, 3), testing::Eq(true));
  
  EXPECT_THAT(is_site_missing_at(tree, x, 0), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, x, 1), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, x, 2), testing::Eq(true));
  EXPECT_THAT(is_site_missing_at(tree, x, 3), testing::Eq(true));
  
  EXPECT_THAT(is_site_missing_at(tree, a, 0), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, a, 1), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, a, 2), testing::Eq(true));
  EXPECT_THAT(is_site_missing_at(tree, a, 3), testing::Eq(true));
  
  EXPECT_THAT(is_site_missing_at(tree, b, 0), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, b, 1), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, b, 2), testing::Eq(true));
  EXPECT_THAT(is_site_missing_at(tree, b, 3), testing::Eq(true));
  
  EXPECT_THAT(is_site_missing_at(tree, c, 0), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, c, 1), testing::Eq(true));
  EXPECT_THAT(is_site_missing_at(tree, c, 2), testing::Eq(false));
  EXPECT_THAT(is_site_missing_at(tree, c, 3), testing::Eq(true));
}

TEST(Phylo_tree_calc_test, calc_state_frequencies_per_partition_of) {
  auto seq = Real_sequence         {rA, rC, rC, rG, rG, rG, rT, rT, rT, rT};
  auto evo = make_global_evo_model({ 0,  1,  0,  1,  2,  0,  1,  2,  3,  0});
  EXPECT_THAT(calc_state_frequencies_per_partition_of(seq, evo), testing::Eq(
      Partition_vector<Seq_vector<int>>{
        Seq_vector{1, 1, 1, 1},  // Partition 0: A-C--G---T
        Seq_vector{0, 1, 1, 1},  // Partition 1: -C-G--T---
        Seq_vector{0, 0, 1, 1},  // Partition 2: ----G--T--
        Seq_vector{0, 0, 0, 1},  // Partition 3: --------T-
      }));
}

TEST_F(Phylo_tree_calc_complex_test, calc_site_state_at) {
  // First at nodes
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 0), testing::Eq(rC));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 1), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 2), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 3), testing::Eq(rA));

  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(b), 0), testing::Eq(rT));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(b), 1), testing::Eq(rG));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(b), 2), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(b), 3), testing::Eq(rA));

  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(c), 0), testing::Eq(rG));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(c), 1), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(c), 2), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(c), 3), testing::Eq(rA));

  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(x), 0), testing::Eq(rT));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 1), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 2), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(a), 3), testing::Eq(rA));

  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(r), 0), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(r), 1), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(r), 2), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, tree.node_loc(r), 3), testing::Eq(rA));

  // Then around mutations
  EXPECT_THAT(calc_site_state_at(tree, {a,  0.4}, 0), testing::Eq(rT));
  EXPECT_THAT(calc_site_state_at(tree, {a,  0.6}, 0), testing::Eq(rC));
  
  EXPECT_THAT(calc_site_state_at(tree, {b,  0.9}, 1), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, {b,  1.1}, 1), testing::Eq(rG));
  
  EXPECT_THAT(calc_site_state_at(tree, {c, -0.1}, 0), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, {c,  0.1}, 0), testing::Eq(rT));
  
  EXPECT_THAT(calc_site_state_at(tree, {c,  0.9}, 0), testing::Eq(rT));
  EXPECT_THAT(calc_site_state_at(tree, {c,  1.0}, 0), testing::Eq(rG));

  EXPECT_THAT(calc_site_state_at(tree, {x, -0.6}, 0), testing::Eq(rA));
  EXPECT_THAT(calc_site_state_at(tree, {x, -0.4}, 0), testing::Eq(rT));
}

TEST_F(Phylo_tree_calc_complex_test, calc_T) {
  auto result = calc_T(tree);

  auto expected = 0.0;
  expected += 1.0;   // r->x branch
  expected += 1.0;   // x->a branch
  expected += 2.0;   // x->b branch
  expected += 4.0;   // r->c branch

  EXPECT_THAT(result, testing::DoubleNear(expected, 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_Ttwiddle_beta_a) {
  auto results = calc_Ttwiddle_beta_a(tree, evo);

  auto expected = Partition_vector<Seq_vector<double>>(evo.num_partitions(), Seq_vector<double>{});

  // Site 0 - partition 0
  expected[0][rA] += 0.5*nu_l(0);  // r->x branch up to A0T
  expected[0][rT] += 0.5*nu_l(0);  // r->x branch from A0T
  expected[0][rT] += 0.5*nu_l(0);  // x->a branch to T0C
  expected[0][rC] += 0.5*nu_l(0);  // x->a branch from T0C
  expected[0][rT] += 2.0*nu_l(0);  // x->b branch
  expected[0][rA] += 1.0*nu_l(0);  // r->c branch to A0T
  expected[0][rT] += 1.0*nu_l(0);  // r->c branch from A0T to T0G
  expected[0][rG] += 2.0*nu_l(0);  // r->c branch from T0G

  // Site 1 - partition 1
  expected[1][rA] += 1.0*nu_l(1);  // r->x branch
  expected[1][rA] += 1.0*nu_l(1);  // x->a branch
  expected[1][rA] += 1.0*nu_l(1);  // x->b branch to A1G
  expected[1][rG] += 1.0*nu_l(1);  // x->b branch from A1G
  //expected[1][rA] += 4.0*nu_l(1);  // r->c branch: missing!

  // Site 2 - partition 0
  //expected[0][rA] += 1.0*nu_l(2);  // r->x branch: missing!
  //expected[0][rA] += 1.0*nu_l(2);  // x->a branch: missing!
  //expected[0][rA] += 2.0*nu_l(2);  // x->b branch: missing!
  expected[0][rA] += 4.0*nu_l(2);  // r->c branch

  // Site 3 - partition 1
  //expected[1][rA] += 1.0*nu_l(3);  // r->x branch: missing!
  //expected[1][rA] += 1.0*nu_l(3);  // x->a branch: missing!
  //expected[1][rA] += 2.0*nu_l(3);  // x->b branch: missing!
  //expected[1][rA] += 4.0*nu_l(3);  // r->c branch: missing!
  
  EXPECT_THAT(results, testing::Pointwise(pointwise_seq_vector_double_near(1e-6), expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_T_l_a) {
  auto expected = std::vector(4, Seq_vector{0.0, 0.0, 0.0, 0.0});

  // Site 0
  expected[0][rA] += 0.5;  // r->x branch up to A0T mutation
  expected[0][rT] += 0.5;  // r->x branch from A0T mutation
  expected[0][rT] += 0.5;  // x->a branch up to T0C mutation
  expected[0][rC] += 0.5;  // x->a branch from T0C mutation
  expected[0][rT] += 2.0;  // x->b branch
  expected[0][rA] += 1.0;  // r->c branch up to A0T mutation
  expected[0][rT] += 1.0;  // r->c branch from A0T mutation to T0G mutation
  expected[0][rG] += 2.0;  // r->c branch from A0G mutation

  // Site 1
  expected[1][rA] += 1.0;  // r->x branch
  expected[1][rA] += 1.0;  // x->a branch
  expected[1][rA] += 1.0;  // x->b branch up to A1G mutation
  expected[1][rG] += 1.0;  // x->b branch from A1G mutation
  // Site 1 missing below r->c attachment point

  // Site 2
  //Site 2 missing below r->x attachment point
  expected[2][rA] += 4.0;  // r->c branch

  // Site 3
  //Site 3 missing everywhere!

  EXPECT_THAT(calc_T_l_a(tree), testing::Pointwise(pointwise_seq_vector_double_near(1e-6), expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_Ttwiddle_l) {
  auto T_l_a = calc_T_l_a(tree);
  auto expected = std::vector(tree.num_sites(), 0.0);
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    expected[l] =
        q_l_a(l, rA) * T_l_a[l][rA] +
        q_l_a(l, rC) * T_l_a[l][rC] +
        q_l_a(l, rG) * T_l_a[l][rG] +
        q_l_a(l, rT) * T_l_a[l][rT];
  }

  EXPECT_THAT(calc_Ttwiddle_l(tree, evo), testing::Pointwise(pointwise_double_near(1e-6), expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_Ttwiddle_a) {
  auto nu_l = std::array{0.1, 0.2, 0.3, 0.4};
  auto T_l_a = calc_T_l_a(tree);

  auto expected = Seq_vector{0.0, 0.0, 0.0, 0.0};
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    for (auto a : k_all_real_seq_letters) {
      expected[a] += nu_l[l] * T_l_a[l][a];
    }
  }

  EXPECT_THAT(calc_Ttwiddle_a(tree, nu_l), seq_vector_double_near(expected, 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_dTtwiddle_a_dt_for_sequence) {
  auto nu_l = std::array{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  EXPECT_THAT(calc_dTtwiddle_a_dt_for_sequence(Real_sequence{rA, rC, rC, rG, rG, rG, rT, rT, rT, rT}, nu_l),
              seq_vector_double_near(Seq_vector{
                  nu_l[0],                                // 1 A
                  nu_l[1] + nu_l[2],                      // 2 C's
                  nu_l[3] + nu_l[4] + nu_l[5],            // 3 G's
                  nu_l[6] + nu_l[7] + nu_l[8] + nu_l[9]   // 4 T's
                }, 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_log_root_prior_zero_p_is_ok) {
  // If pi_a[a] == 0.0, then either n_a[a] == 0 too, or the root sequence is impossible
  evo.partition_evo_model[0].pi_a = Seq_vector{0.3, 0.7, 0.0, 0.0};  // A-A-
  evo.partition_evo_model[1].pi_a = Seq_vector{0.3, 0.0, 0.7, 0.0};  // -A-N
  
  auto result = calc_log_root_prior(tree, evo);

  auto expected =
      std::log(pi_l_a(0, rA)) +
      std::log(pi_l_a(1, rA)) +
      std::log(pi_l_a(2, rA)) +
      0.0;  // Note std::log(pi_l_a(3, rA)), site 3 is missing

  EXPECT_THAT(result, testing::DoubleNear(expected, 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_log_root_prior_zero_p_is_impossible) {
  // If pi_a[a] == 0.0, then either n_a[a] == 0 too, or the root sequence is impossible
  evo.partition_evo_model[0].pi_a = Seq_vector{0.0, 0.3, 0.7, 0.0};  // A-A-
  evo.partition_evo_model[1].pi_a = Seq_vector{0.3, 0.0, 0.7, 0.0};  // -A-N

  auto result = calc_log_root_prior(tree, evo);

  EXPECT_THAT(result, testing::Eq(-std::numeric_limits<double>::infinity()));
}

TEST_F(Phylo_tree_calc_complex_test, calc_log_G_below_root) {
  auto result = calc_log_G_below_root(tree, evo);

  auto expected = 0.0;

  // Site 0
  expected +=
      // r-x
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      
      // x-a
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rC)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rC)*0.5 +
      
      // x-b
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0 +
      
      // r-c
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rG)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rG)*2.0;

  // Site 1
  expected +=
      // r-x
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      
      // x-a
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      
      // x-b
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rA, rG)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rG)*1.0 +

      // r-c (missing)
      0.0;

  // Site 2
  expected +=
      // r-x (missing)
      0.0 +
      
      // x-a (missing)
      0.0 +
      
      // x-b (missing)
      0.0 +

      // r-c
      -mu_l(2)*nu_l(2)*q_l_a(2, rA)*4.0;

  ASSERT_THAT(result, testing::DoubleNear(expected, 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_muts) {
  EXPECT_THAT(calc_num_muts(tree), testing::Eq(5));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_muts_ab) {
  auto result = calc_num_muts_ab(tree);

  auto expected = Seq_matrix<int>{0};
  ++expected[rA][rT];  // A0T on r->x branch
  ++expected[rT][rC];  // T0C on x->a branch
  ++expected[rA][rG];  // A1G on x->b branch
  ++expected[rA][rT];  // A0T on r->c branch
  ++expected[rT][rG];  // T0G on r->c branch

  EXPECT_THAT(result, testing::Eq(expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_muts_beta_ab) {
  auto result = calc_num_muts_beta_ab(tree, evo);

  auto expected = Partition_vector<Seq_matrix<int>>(evo.num_partitions(), Seq_matrix<int>{});
  ++expected[0][rA][rT];  // A0T on r->x branch (site 0 = partition 0)
  ++expected[0][rT][rC];  // T0C on x->a branch (site 0 = partition 0)
  ++expected[1][rA][rG];  // A1G on x->b branch (site 1 = partition 1)
  ++expected[0][rA][rT];  // A0T on r->c branch (site 0 = partition 0)
  ++expected[0][rT][rG];  // T0G on r->c branch (site 0 = partition 0)

  EXPECT_THAT(result, testing::Eq(expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_muts_l) {
  auto result = calc_num_muts_l(tree);

  auto expected = Node_vector<int>{
    4,  // A0T on r->x branch, T0C on x->a branch, A0T and T0G on r->c branch
    1,  // A1G on x->b branch
    0,  // (C2A above root is not a real mutation)
    0
  };

  EXPECT_THAT(result, testing::Eq(expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_muts_l_ab) {
  auto result = calc_num_muts_l_ab(tree);

  auto expected = std::vector{Seq_matrix<int>{}, Seq_matrix<int>{}, Seq_matrix<int>{}, Seq_matrix<int>{}};
  ++expected[0][rA][rT];  // A0T on r->x branch
  ++expected[0][rT][rC];  // T0C on x->a branch
  ++expected[1][rA][rG];  // A1G on x->b branch
  ++expected[0][rA][rT];  // A0T on r->c branch
  ++expected[0][rT][rG];  // T0G on r->c branch

  EXPECT_THAT(result, testing::Eq(expected));
}

TEST_F(Phylo_tree_calc_complex_test, calc_num_sites_missing_at_every_node) {
  auto result = calc_num_sites_missing_at_every_node(tree);

  EXPECT_THAT(result[r], testing::Eq(1));     // r: AAAN
  EXPECT_THAT(result[x], testing::Eq(2));     // x: TANN
  EXPECT_THAT(result[a], testing::Eq(2));     // a: CANN
  EXPECT_THAT(result[b], testing::Eq(2));     // b: TGNN
  EXPECT_THAT(result[c], testing::Eq(2));     // c: GNAN
}

TEST_F(Phylo_tree_calc_complex_test, recalc_num_sites_missing_upstream) {
  auto truth = calc_num_sites_missing_at_every_node(tree);

  for (const auto& ancestor : {k_no_node, r, x, a, b, c}) {
    SCOPED_TRACE(absl::StrFormat("ancestor = %s", ancestor == k_no_node ? "-" : tree.at(ancestor).name));
    
    for (const auto& node : {r, x, a, b, c}) {
      SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
      
      if (not descends_from(tree, node, ancestor) || node == ancestor) { continue; }

      auto result = Node_vector<int>(std::ssize(tree), -100);
      result[node] = reconstruct_missing_sites_at(tree, node).num_sites();

      recalc_num_sites_missing_upstream(tree, node, ancestor, result);

      for (auto cur = node; cur != ancestor; cur = tree.at(cur).parent) {
        SCOPED_TRACE(absl::StrFormat("cur = %s", tree.at(cur).name));
        EXPECT_THAT(result[cur], testing::Eq(truth[cur]));
      }
    }
  }
}

TEST_F(Phylo_tree_calc_complex_test, calc_cum_Q_l_for_sequence) {
  auto result = calc_cum_Q_l_for_sequence({rA, rC, rG, rT}, evo);

  auto expected_cum_Q_l = std::vector<double>{
    0.0,                               // no site before Site 0
    mu_l(0) * nu_l(0) * q_l_a(0, rA),  // Site 0: A
    mu_l(1) * nu_l(1) * q_l_a(1, rC),  // Site 1: C
    mu_l(2) * nu_l(2) * q_l_a(2, rG),  // Site 2: G
    mu_l(3) * nu_l(3) * q_l_a(3, rT)   // Site 3: T
  };
  std::partial_sum(expected_cum_Q_l.begin(), expected_cum_Q_l.end(), expected_cum_Q_l.begin());
    
  EXPECT_THAT(result, testing::Pointwise(testing::DoubleNear(1e-6), expected_cum_Q_l));
}

TEST_F(Phylo_tree_calc_complex_test, calc_lambda_for_sequence) {
  auto result = calc_lambda_for_sequence({rA, rC, rG, rT}, evo);

  EXPECT_THAT(result, testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rA) +  // Site 0: A
      mu_l(1) * nu_l(1) * q_l_a(1, rC) +  // Site 1: C
      mu_l(2) * nu_l(2) * q_l_a(2, rG) +  // Site 2: G
      mu_l(3) * nu_l(3) * q_l_a(3, rT)    // Site 3: T
      , 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_lambda_i) {
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);

  // r
  EXPECT_THAT(lambda_i[r], testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rA) +  // Site 0: A
      mu_l(1) * nu_l(1) * q_l_a(1, rA) +  // Site 1: A
      mu_l(2) * nu_l(2) * q_l_a(2, rA) +  // Site 2: A
      0.0                                 // Site 3: N
      , 1e-6));

  // x
  EXPECT_THAT(lambda_i[x], testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rT) +  // Site 0: T
      mu_l(1) * nu_l(1) * q_l_a(1, rA) +  // Site 1: A
      0.0 +                               // Site 2: N
      0.0                                 // Site 3: N
      , 1e-6));

  // a
  EXPECT_THAT(lambda_i[a], testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rC) +  // Site 0: C
      mu_l(1) * nu_l(1) * q_l_a(1, rA) +  // Site 1: A
      0.0 +                               // Site 2: N
      0.0                                 // Site 3: N
      , 1e-6));

  // b
  EXPECT_THAT(lambda_i[b], testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rT) +  // Site 0: T
      mu_l(1) * nu_l(1) * q_l_a(1, rG) +  // Site 1: G
      0.0 +                               // Site 2: N
      0.0                                 // Site 3: N
      , 1e-6));

  // c
  EXPECT_THAT(lambda_i[c], testing::DoubleNear(
      mu_l(0) * nu_l(0) * q_l_a(0, rG) +  // Site 0: G
      0.0 +                               // Site 1: N
      mu_l(2) * nu_l(2) * q_l_a(2, rA) +  // Site 2: A
      0.0                                 // Site 3: N
      , 1e-6));

  // All from scratch
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
    EXPECT_THAT(lambda_i[node], testing::DoubleNear(
        calc_lambda_at_node(tree, node, evo), 1e-6));
  }
}

TEST_F(Phylo_tree_calc_complex_test, recalc_lambda_i_upstream) {
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  for (const auto& ancestor : {k_no_node, r, x, a, b, c}) {
    SCOPED_TRACE(absl::StrFormat("ancestor = %s", ancestor == k_no_node ? "-" : tree.at(ancestor).name));
    
    for (const auto& node : {r, x, a, b, c}) {
      SCOPED_TRACE(absl::StrFormat("node = %s", tree.at(node).name));
      
      if (not descends_from(tree, node, ancestor) || node == ancestor) { continue; }

      auto result = Node_vector<double>(std::ssize(tree), 1.23456);
      result[node] = calc_lambda_at_node(tree, node, evo);

      recalc_lambda_i_upstream(tree, node, ancestor, evo, result, ref_cum_Q_l);

      for (auto cur = node; cur != ancestor; cur = tree.at(cur).parent) {
        SCOPED_TRACE(absl::StrFormat("cur = %s", tree.at(cur).name));
        EXPECT_THAT(result[cur], testing::DoubleNear(
            calc_lambda_at_node(tree, cur, evo), 1e-6));
      }
    }
  }
}

TEST_F(Phylo_tree_calc_complex_test, calc_path_log_G) {
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto state_frequencies_of_ref_sequence_per_partition =
      calc_state_frequencies_per_partition_of(ref_sequence, evo);

  // Y -> Y for all Y
  for (const auto& Y : index_order_traversal(tree)) {
    SCOPED_TRACE(absl::StrFormat("Y = %s", tree.at(Y).name));
    EXPECT_THAT(calc_path_log_G(tree, Y, Y, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
                testing::Eq(0.0));
  }
  
  // - -> r
  EXPECT_THAT(calc_path_log_G(tree, k_no_node, r, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // root prior for AAAN
                  std::log(pi_l_a(0, rA)) +
                  std::log(pi_l_a(1, rA)) +
                  std::log(pi_l_a(2, rA)) +
                  0.0
                  , 1e-6));

  // - -> x
  EXPECT_THAT(calc_path_log_G(tree, k_no_node, x, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // root prior for AAAN
                  std::log(pi_l_a(0, rA)) +
                  std::log(pi_l_a(1, rA)) +
                  std::log(pi_l_a(2, rA)) +
                  0.0
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  , 1e-6));

  // - -> a
  EXPECT_THAT(calc_path_log_G(tree, k_no_node, a, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // root prior for AAAN
                  std::log(pi_l_a(0, rA)) +
                  std::log(pi_l_a(1, rA)) +
                  std::log(pi_l_a(2, rA)) +
                  0.0
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  // x-a
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rT,rC)) - mu_l(0)*nu_l(0)*q_l_a(0, rC)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  , 1e-6));

  // - -> b
  EXPECT_THAT(calc_path_log_G(tree, k_no_node, b, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // root prior for AAAN
                  std::log(pi_l_a(0, rA)) +
                  std::log(pi_l_a(1, rA)) +
                  std::log(pi_l_a(2, rA)) +
                  0.0
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  // x-b
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 + std::log(mu_l(1)*nu_l(1)*q_l_ab(1,rA,rG)) - mu_l(1)*nu_l(1)*q_l_a(1, rG)*1.0
                  , 1e-6));

  // - -> c
  EXPECT_THAT(calc_path_log_G(tree, k_no_node, c, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // root prior for AAAN
                  std::log(pi_l_a(0, rA)) +
                  std::log(pi_l_a(1, rA)) +
                  std::log(pi_l_a(2, rA)) +
                  0.0
                  // r-c
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT))
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*1.0 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rT,rG))
                  - mu_l(0)*nu_l(0)*q_l_a(0, rG)*2.0
                  - mu_l(2)*nu_l(2)*q_l_a(2, rA)*4.0
                  , 1e-6));

  // r -> x
  EXPECT_THAT(calc_path_log_G(tree, r, x, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  , 1e-6));

  // r -> a
  EXPECT_THAT(calc_path_log_G(tree, r, a, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  // x-a
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rT,rC)) - mu_l(0)*nu_l(0)*q_l_a(0, rC)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  , 1e-6));

  // r -> b
  EXPECT_THAT(calc_path_log_G(tree, r, b, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // r-x
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT)) - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  // x-b
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 + std::log(mu_l(1)*nu_l(1)*q_l_ab(1,rA,rG)) - mu_l(1)*nu_l(1)*q_l_a(1, rG)*1.0
                  , 1e-6));

  // r -> c
  EXPECT_THAT(calc_path_log_G(tree, r, c, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // r-c
                  - mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rA,rT))
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*1.0 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rT,rG))
                  - mu_l(0)*nu_l(0)*q_l_a(0, rG)*2.0
                  - mu_l(2)*nu_l(2)*q_l_a(2, rA)*4.0
                  , 1e-6));

  // x -> a
  EXPECT_THAT(calc_path_log_G(tree, x, a, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // x-a
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 + std::log(mu_l(0)*nu_l(0)*q_l_ab(0,rT,rC)) - mu_l(0)*nu_l(0)*q_l_a(0, rC)*0.5
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0
                  , 1e-6));

  // x -> b
  EXPECT_THAT(calc_path_log_G(tree, x, b, evo, lambda_i, state_frequencies_of_ref_sequence_per_partition),
              testing::DoubleNear(
                  // x-b
                  - mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0
                  - mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 + std::log(mu_l(1)*nu_l(1)*q_l_ab(1,rA,rG)) - mu_l(1)*nu_l(1)*q_l_a(1, rG)*1.0
                  , 1e-6));
}

TEST_F(Phylo_tree_calc_complex_test, calc_max_tip_time) {
  EXPECT_THAT(calc_max_tip_time(tree), testing::Eq(tree.at(c).t));
}

}  // namespace delphy
