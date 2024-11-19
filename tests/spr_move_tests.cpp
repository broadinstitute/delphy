#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "spr_move.h"

#include "evo_model.h"
#include "phylo_tree_calc.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Spr_move_test_base : public testing::Test {
 protected:
  Real_sequence ref_sequence{rA, rC, rA, rA};
  
  Global_evo_model evo = make_global_evo_model({0, 1, 0, 1});  // 2 site partitions
  double mu_JC = 0.125;  // For proposing mutational histories
  
  Spr_move_test_base() {
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

  auto P_JC(Real_seq_letter a, Real_seq_letter b, double t) -> double {
    if (a == b) { return 1.0 + 3./4*std::expm1(-4./3.*mu_JC*t); }
    else { return -1./4.*std::expm1(-4./3.*mu_JC*t); }
  }
  
  auto run_propose_new_graft_test(Phylo_tree& tree) -> void;
  
  auto run_full_spr_move_test(
      Phylo_tree& tree,
      Node_index X,
      Node_index SS,
      double new_t_P,
      absl::BitGenRef bitgen)
      -> void;
};

class Spr_move_simple_test : public Spr_move_test_base {
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
  // (ACAA) A3N- C1A --+ r (AAAN)     +-------- A1G ------- b (TGNN)
  //                     |
  //                     +A1N--------A0T------- T0G ------------------ c (GNAN)
  // 

  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Spr_move_simple_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;

    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
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

TEST_F(Spr_move_simple_test, analyze_graft_a) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(a);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(1));

  const auto& b0 = analysis.branch_infos[0];
  EXPECT_THAT(b0.A, testing::Eq(x));
  EXPECT_THAT(b0.B, testing::Eq(a));
  EXPECT_THAT(b0.is_open, testing::IsFalse());
  EXPECT_THAT(b0.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b0.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b0.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rC) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(estd::ranges::to_vec(b0.warm_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(b0.hot_muts_to_X, testing::ElementsAre(Mutation{rT, 0, rC, 0.5}));
  EXPECT_THAT(b0.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rC)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -std::log(P_JC(rT, rC, 1.0)) +
      // Site 1
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rA, 1.0)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rC)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rC)*0.5 +
      // Site 1
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_simple_test, analyze_graft_b) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(b);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(1));

  const auto& b0 = analysis.branch_infos[0];
  EXPECT_THAT(b0.A, testing::Eq(x));
  EXPECT_THAT(b0.B, testing::Eq(b));
  EXPECT_THAT(b0.is_open, testing::IsFalse());
  EXPECT_THAT(b0.T_to_X, testing::Eq(2.0));
  EXPECT_THAT(b0.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b0.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rG),
      1e-6));
  EXPECT_THAT(estd::ranges::to_vec(b0.warm_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(b0.hot_muts_to_X, testing::ElementsAre(Mutation{rA, 1, rG, 1.0}));
  EXPECT_THAT(b0.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(1, rA, rG)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0
      -mu_JC*2.0 +
      -std::log(P_JC(rT, rT, 2.0)) +
      // Site 1
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rG, 2.0)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0 +
      // Site 1
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rA, rG)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rG)*1.0 +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_simple_test, analyze_graft_c) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(c);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // Site 2 is only present in P->X
  const auto& b_P_X = analysis.branch_infos[Spr_graft::k_branch_info_P_X];
  EXPECT_THAT(b_P_X.A, testing::Eq(r));
  EXPECT_THAT(b_P_X.B, testing::Eq(c));
  EXPECT_THAT(b_P_X.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_X.T_to_X, testing::Eq(4.0));
  EXPECT_THAT(b_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_X.warm_sites, testing::ElementsAre(Site_interval{2, 3}));
  EXPECT_THAT(b_P_X.hot_sites, testing::ElementsAre(Site_interval{2, 3}));
  EXPECT_THAT(b_P_X.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_X.hot_deltas_to_X, testing::IsEmpty());

  // Site 1 is only present in P->S
  const auto& b_P_S = analysis.branch_infos[Spr_graft::k_branch_info_P_S];
  EXPECT_THAT(b_P_S.A, testing::Eq(r));
  EXPECT_THAT(b_P_S.B, testing::Eq(x));
  EXPECT_THAT(b_P_S.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_S.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b_P_S.partial_lambda_at_A, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_P_S.partial_lambda_at_X, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_P_S.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b_P_S.hot_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b_P_S.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_S.hot_deltas_to_X, testing::IsEmpty());

  // Site 0 is present all the way along S->P->X
  const auto& b_S_P_X = analysis.branch_infos[Spr_graft::k_branch_info_S_P_X];
  EXPECT_THAT(b_S_P_X.A, testing::Eq(x));
  EXPECT_THAT(b_S_P_X.B, testing::Eq(r));
  EXPECT_THAT(b_S_P_X.is_open, testing::IsFalse());
  EXPECT_THAT(b_S_P_X.T_to_X, testing::Eq(5.0));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT),
      1e-6));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rG),
      1e-6));
  EXPECT_THAT(b_S_P_X.warm_sites, testing::ElementsAre(Site_interval{0, 1}));
  EXPECT_THAT(b_S_P_X.hot_sites, testing::ElementsAre(Site_interval{0, 1}));
  EXPECT_THAT(b_S_P_X.hot_muts_to_X, testing::ElementsAre(
      Mutation{rT, 0, rA, -1.5},
      Mutation{rA, 0, rT,  0.0},
      Mutation{rT, 0, rG,  1.0}));
  EXPECT_THAT(b_S_P_X.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rG)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0 (S->P->X)
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*2.0 +
      -std::log(P_JC(rT, rG, 5.0)) +
      // Site 1 (P->S)
      -mu_JC*1.0 +
      // Site 2 (P->X)
      -mu_JC*4.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0 (P->S part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      // Site 0 (P->X part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rG)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rG)*2.0 +
      // Site 0 (root sequence change)
      std::log(pi_l_a(0, rA) / pi_l_a(0, rT)) +
      
      // Site 1 (P->S)
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      // Site 2 (P->X)
      -mu_l(2)*nu_l(2)*q_l_a(2, rA)*4.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_simple_test, analyze_graft_x) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(x);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // Site 1 is only present in P->X
  const auto& b_P_X = analysis.branch_infos[Spr_graft::k_branch_info_P_X];
  EXPECT_THAT(b_P_X.A, testing::Eq(r));
  EXPECT_THAT(b_P_X.B, testing::Eq(x));
  EXPECT_THAT(b_P_X.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_X.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_P_X.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b_P_X.hot_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b_P_X.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_X.hot_deltas_to_X, testing::IsEmpty());

  // Site 2 is only present in P->S
  const auto& b_P_S = analysis.branch_infos[Spr_graft::k_branch_info_P_S];
  EXPECT_THAT(b_P_S.A, testing::Eq(r));
  EXPECT_THAT(b_P_S.B, testing::Eq(c));
  EXPECT_THAT(b_P_S.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_S.T_to_X, testing::Eq(4.0));
  EXPECT_THAT(b_P_S.partial_lambda_at_A, testing::DoubleNear(
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_S.partial_lambda_at_X, testing::DoubleNear(
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_S.warm_sites, testing::ElementsAre(Site_interval{2, 3}));
  EXPECT_THAT(b_P_S.hot_sites, testing::ElementsAre(Site_interval{2, 3}));
  EXPECT_THAT(b_P_S.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_S.hot_deltas_to_X, testing::IsEmpty());

  // Site 0 is present all the way along S->P->X
  const auto& b_S_P_X = analysis.branch_infos[Spr_graft::k_branch_info_S_P_X];
  EXPECT_THAT(b_S_P_X.A, testing::Eq(c));
  EXPECT_THAT(b_S_P_X.B, testing::Eq(r));
  EXPECT_THAT(b_S_P_X.is_open, testing::IsFalse());
  EXPECT_THAT(b_S_P_X.T_to_X, testing::Eq(5.0));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rG),
      1e-6));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT),
      1e-6));
  EXPECT_THAT(b_S_P_X.warm_sites, testing::ElementsAre(Site_interval{0, 1}));
  EXPECT_THAT(b_S_P_X.hot_sites, testing::ElementsAre(Site_interval{0, 1}));
  EXPECT_THAT(b_S_P_X.hot_muts_to_X, testing::ElementsAre(
      Mutation{rG, 0, rT, -3.0},
      Mutation{rT, 0, rA, -2.0},
      Mutation{rA, 0, rT, -0.5}));
  EXPECT_THAT(b_S_P_X.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(0, rG, rT)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0 (S->P->X)
      -mu_JC*2.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -std::log(P_JC(rG, rT, 5.0)) +
      // Site 1 (P->X)
      -mu_JC*1.0 +
      // Site 2 (P->S)
      -mu_JC*4.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0 (P->S part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*1.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rG)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rG)*2.0 +
      // Site 0 (P->X part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      // Site 0 (root sequence change)
      std::log(pi_l_a(0, rA) / pi_l_a(0, rG)) +
      
      // Site 1 (P->X)
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      // Site 2 (P->S)
      -mu_l(2)*nu_l(2)*q_l_a(2, rA)*4.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

auto Spr_move_test_base::run_propose_new_graft_test(Phylo_tree& tree) -> void {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto old_tree = Phylo_tree{tree};  // Copy
  
  for (auto seed = 0; seed != 1000; ++seed) {
    auto bitgen = std::mt19937_64(seed + 12345);
    SCOPED_TRACE(absl::StrFormat("Seed %d", seed));

    // Restore old tree
    tree = old_tree;  // Copy
    lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
    num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);

    // Pick a random non-root node
    auto X = tree.root;
    while (X == tree.root) {
      X = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0, std::ssize(tree));
    }
    SCOPED_TRACE(absl::StrFormat("X = %d", X));

    // Analyze the old situation
    auto old_graft = spr.analyze_graft(X);
    auto old_log_G = calc_log_root_prior(tree, evo) + calc_log_G_below_root(tree, evo);

    // Without moving X, propose and apply new changes
    spr.peel_graft(old_graft);
    auto new_graft = spr.propose_new_graft(X, bitgen);
    spr.apply_graft(new_graft);

    // Check that everything worked out
    assert_phylo_tree_integrity(tree);
    
    // Analyze the new situation
    auto new_graft_redux = spr.analyze_graft(X);
    auto new_log_G = calc_log_root_prior(tree, evo) + calc_log_G_below_root(tree, evo);

    EXPECT_THAT(new_graft_redux.X, testing::Eq(new_graft.X));
    EXPECT_THAT(new_graft_redux.S, testing::Eq(new_graft.S));
    EXPECT_THAT(new_graft_redux.t_P, testing::DoubleNear(new_graft.t_P, 1e-6));
    EXPECT_THAT(new_graft_redux.branch_infos, testing::SizeIs(std::ssize(new_graft.branch_infos)));
    for (auto i = 0; i != std::ssize(new_graft.branch_infos); ++i) {
      SCOPED_TRACE(absl::StrFormat("i = %d", i));
      const auto& bi = new_graft.branch_infos[i];
      const auto& bi_redux = new_graft_redux.branch_infos[i];

      EXPECT_THAT(bi_redux.A, testing::Eq(bi.A));
      EXPECT_THAT(bi_redux.B, testing::Eq(bi.B));
      EXPECT_THAT(bi_redux.is_open, testing::Eq(bi.is_open));
      EXPECT_THAT(bi_redux.T_to_X, testing::DoubleNear(bi.T_to_X, 1e-6));
      EXPECT_THAT(bi_redux.partial_lambda_at_A, testing::DoubleNear(bi.partial_lambda_at_A, 1e-6));
      EXPECT_THAT(bi_redux.partial_lambda_at_X, testing::DoubleNear(bi.partial_lambda_at_X, 1e-6));
      EXPECT_THAT(bi_redux.warm_sites, testing::Eq(bi_redux.warm_sites));
      EXPECT_THAT(bi_redux.hot_sites, testing::Eq(bi_redux.hot_sites));
      EXPECT_THAT(bi_redux.hot_muts_to_X, testing::SizeIs(std::ssize(bi.hot_muts_to_X)));
      for (auto j = 0; j != std::ssize(bi.hot_muts_to_X); ++j) {
        const auto& m = bi.hot_muts_to_X[j];
        const auto& m_redux = bi_redux.hot_muts_to_X[j];
        SCOPED_TRACE(absl::StrFormat("j = %d, m = %s, m_redux = %s",
                                     j, absl::FormatStreamed(m), absl::FormatStreamed(m_redux)));
        EXPECT_THAT(m.from, testing::Eq(m_redux.from));
        EXPECT_THAT(m.to, testing::Eq(m_redux.to));
        EXPECT_THAT(m.site, testing::Eq(m_redux.site));
        EXPECT_THAT(m.t, testing::DoubleNear(m_redux.t, 1e-6));
      }
      EXPECT_THAT(bi_redux.hot_deltas_to_X, testing::Eq(bi.hot_deltas_to_X));
    }
    EXPECT_THAT(new_graft_redux.delta_log_G, testing::DoubleNear(new_graft.delta_log_G, 1e-6));
    EXPECT_THAT(new_graft_redux.log_alpha_mut, testing::DoubleNear(new_graft.log_alpha_mut, 1e-6));
    
    EXPECT_THAT(new_log_G, testing::DoubleNear(old_log_G - old_graft.delta_log_G + new_graft.delta_log_G, 1e-6));

    for (const auto& node : index_order_traversal(tree)) {
      if (tree.at(node).is_tip()) {
        for (auto l = 0; l != tree.num_sites(); ++l) {
          auto missing_now = is_site_missing_at(tree, node, l);
          EXPECT_THAT(missing_now, testing::Eq(is_site_missing_at(old_tree, node, l)));
          if (not missing_now) {
            EXPECT_THAT(view_of_sequence_at(tree, node)[l],
                        testing::Eq(view_of_sequence_at(old_tree, node)[l]));
          }
        }
      }
    }
    EXPECT_THAT(lambda_i, testing::Pointwise(testing::DoubleNear(1e-6), calc_lambda_i(tree, evo, ref_cum_Q_l)));
    EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
  }
}

TEST_F(Spr_move_simple_test, propose_new_graft) {
  run_propose_new_graft_test(tree);
}

class Spr_move_superfluous_mutation_at_root_test : public Spr_move_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:                 -1.0        0.0        1.0        2.0        3.0
  //                       
  //                         +-- A0T ---+ x (TANN)
  //                         |          
  //  (ACAA) A2N,A3N-- C1A --+ r (AANN)
  //                         |
  //                         +------------------- A0T ------------------ s (TANN)
  // 

  Phylo_tree tree{3};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index s = 2;

  Spr_move_superfluous_mutation_at_root_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;
  
    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, s};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).missations = {{Missation{2, rA}, Missation{3, rA}}, ref_sequence};
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
  
    tree.at(x).parent = r;
    tree.at(x).children = {};
    tree.at(x).name = "x";
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, -0.5}};

    tree.at(s).parent = r;
    tree.at(s).children = {};
    tree.at(s).name = "s";
    tree.at(s).t = 3.0;
    tree.at(s).mutations = {Mutation{rA, 0, rT, 1.0}};
  
    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Spr_move_superfluous_mutation_at_root_test, analyze_graft_x) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(x);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // No site is present only in P->X
  const auto& b_P_X = analysis.branch_infos[Spr_graft::k_branch_info_P_X];
  EXPECT_THAT(b_P_X.A, testing::Eq(r));
  EXPECT_THAT(b_P_X.B, testing::Eq(x));
  EXPECT_THAT(b_P_X.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_X.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b_P_X.partial_lambda_at_A, testing::Eq(0.0));
  EXPECT_THAT(b_P_X.partial_lambda_at_X, testing::Eq(0.0));
  EXPECT_THAT(b_P_X.warm_sites, testing::IsEmpty());
  EXPECT_THAT(b_P_X.hot_sites, testing::IsEmpty());
  EXPECT_THAT(b_P_X.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_X.hot_deltas_to_X, testing::IsEmpty());

  // No site is present only in P->S
  const auto& b_P_S = analysis.branch_infos[Spr_graft::k_branch_info_P_S];
  EXPECT_THAT(b_P_S.A, testing::Eq(r));
  EXPECT_THAT(b_P_S.B, testing::Eq(s));
  EXPECT_THAT(b_P_S.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_S.T_to_X, testing::Eq(4.0));
  EXPECT_THAT(b_P_S.partial_lambda_at_A, testing::Eq(0.0));
  EXPECT_THAT(b_P_S.partial_lambda_at_X, testing::Eq(0.0));
  EXPECT_THAT(b_P_S.warm_sites, testing::IsEmpty());
  EXPECT_THAT(b_P_S.hot_sites, testing::IsEmpty());
  EXPECT_THAT(b_P_S.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b_P_S.hot_deltas_to_X, testing::IsEmpty());

  // Sites 0 and 1 are present all the way along S->P->X
  const auto& b_S_P_X = analysis.branch_infos[Spr_graft::k_branch_info_S_P_X];
  EXPECT_THAT(b_S_P_X.A, testing::Eq(s));
  EXPECT_THAT(b_S_P_X.B, testing::Eq(r));
  EXPECT_THAT(b_S_P_X.is_open, testing::IsFalse());
  EXPECT_THAT(b_S_P_X.T_to_X, testing::Eq(5.0));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b_S_P_X.warm_sites, testing::ElementsAre(Site_interval{0, 2}));
  EXPECT_THAT(b_S_P_X.hot_sites, testing::ElementsAre(Site_interval{0, 2}));
  EXPECT_THAT(b_S_P_X.hot_muts_to_X, testing::ElementsAre(
      Mutation{rT, 0, rA, -3.0},
      Mutation{rA, 0, rT, -0.5}));
  EXPECT_THAT(b_S_P_X.hot_deltas_to_X, testing::IsEmpty());

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0 (S->P->X)
      -mu_JC*2.0 +
      std::log(mu_JC/3) +
      -mu_JC*2.0 +
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -std::log(P_JC(rT, rT, 5.0)) +
      // Site 1 (S->P->X)
      -mu_JC*5.0 +
      -std::log(P_JC(rA, rA, 5.0)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0 (P->S part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*2.0 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*2.0 +
      // Site 0 (P->X part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rA, rT)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      // Site 0 (root sequence change)
      std::log(pi_l_a(0, rA) / pi_l_a(0, rT)) +
      
      // Site 1 (P->S part)
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*4.0 +
      // Site 1 (P->X part)
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_superfluous_mutation_at_root_test, propose_new_graft) {
  run_propose_new_graft_test(tree);
}

class Spr_move_superfluous_mutation_not_at_root_test : public Spr_move_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:                 -1.0        0.0        1.0        2.0        3.0
  //          
  //                                    +-- T0A -- a (AANN)
  //                                    |
  //                         +-- A0T ---+ x (TANN)
  //                         |          |
  //  (ACAA) A2N,A3N-- C1A --+ r (AANN) +-------- T0A ------- b (AANN)
  //                         |
  //                         +------------------------------------------ c (AANN)
  // 
  
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Spr_move_superfluous_mutation_not_at_root_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;
  
    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).missations = {{Missation{2, rA}, Missation{3, rA}}, ref_sequence};
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
  
    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, -0.5}};
  
    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 1.0;
    tree.at(a).mutations = {Mutation{rT, 0, rA, 0.5}};
  
    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 2.0;
    tree.at(b).mutations = {Mutation{rT, 0, rA, 1.0}};
  
    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 3.0;
  
    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Spr_move_superfluous_mutation_not_at_root_test, analyze_graft_a) {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(a);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(1));

  const auto& b0 = analysis.branch_infos[0];
  EXPECT_THAT(b0.A, testing::Eq(x));
  EXPECT_THAT(b0.B, testing::Eq(a));
  EXPECT_THAT(b0.is_open, testing::IsFalse());
  EXPECT_THAT(b0.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b0.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rT) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b0.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA) +
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(estd::ranges::to_vec(b0.warm_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(b0.hot_muts_to_X, testing::ElementsAre(Mutation{rT, 0, rA, 0.5}));
  EXPECT_THAT(b0.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -std::log(P_JC(rT, rA, 1.0)) +
      // Site 1
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rA, 1.0)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0
      -mu_l(0)*nu_l(0)*q_l_a(0, rT)*0.5 +
      std::log(mu_l(0)*nu_l(0)*q_l_ab(0, rT, rA)) +
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*0.5 +
      // Site 1
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_superfluous_mutation_not_at_root_test, propose_new_graft) {
  run_propose_new_graft_test(tree);
}

class Spr_move_precarious_path_without_root_test : public Spr_move_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:                 -1.0        0.0        1.0        2.0        3.0
  //
  //                                               +-------- C1T ------- a (ATNN)  <--- Precarious path r--->a on site 1
  //                                               |
  //                                    +----------+ y (ACNN)
  //                                    |          |
  //                         +-- A1C ---+ x (ACNN) +C1N----------------- b (ANNN)
  //                         |          |
  //  (ACAA) A2N,A3N-- C1A --+ r (AANN) +C1N---------------------------- c (ANNN)
  //                         |
  //                         +------------------------------------------ d (AANN)
  // 
  
  Phylo_tree tree{7};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index y = 2;
  static constexpr Node_index a = 3;
  static constexpr Node_index b = 4;
  static constexpr Node_index c = 5;
  static constexpr Node_index d = 6;

  Spr_move_precarious_path_without_root_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;
  
    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, d};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).missations = {{Missation{2, rA}, Missation{3, rA}}, ref_sequence};
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
  
    tree.at(x).parent = r;
    tree.at(x).children = {y, c};
    tree.at(x).name = "x";
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 1, rC, -0.5}};
  
    tree.at(y).parent = x;
    tree.at(y).children = {a, b};
    tree.at(y).name = "y";
    tree.at(y).t = 1.0;
  
    tree.at(a).parent = y;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 3.0;
    tree.at(a).mutations = {Mutation{rC, 1, rT, 2.0}};
  
    tree.at(b).parent = y;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 3.0;
    tree.at(b).missations = {{Missation{1, rC}}, ref_sequence};
  
    tree.at(c).parent = x;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 3.0;
    tree.at(c).missations = {{Missation{1, rC}}, ref_sequence};
  
    tree.at(d).parent = r;
    tree.at(d).children = {};
    tree.at(d).name = "d";
    tree.at(d).t = 3.0;
  
    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Spr_move_precarious_path_without_root_test, analyze_graft_a) {
  auto can_change_root = false;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(a);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // y->a:  Sites 0 & 1 are warm, and site 0 is hot
  const auto& b0 = analysis.branch_infos[0];
  EXPECT_THAT(b0.A, testing::Eq(y));
  EXPECT_THAT(b0.B, testing::Eq(a));
  EXPECT_THAT(b0.is_open, testing::IsFalse());
  EXPECT_THAT(b0.T_to_X, testing::Eq(2.0));
  EXPECT_THAT(b0.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA),
      1e-6));
  EXPECT_THAT(b0.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA),
      1e-6));
  EXPECT_THAT(estd::ranges::to_vec(b0.warm_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::Contains(0));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::Not(testing::Contains(1)));
  EXPECT_THAT(b0.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b0.hot_deltas_to_X, testing::IsEmpty());

  // x->y: Site 1 is warm but not hot
  const auto& b1 = analysis.branch_infos[1];
  EXPECT_THAT(b1.A, testing::Eq(x));
  EXPECT_THAT(b1.B, testing::Eq(y));
  EXPECT_THAT(b1.is_open, testing::IsFalse());
  EXPECT_THAT(b1.T_to_X, testing::Eq(3.0));
  EXPECT_THAT(b1.partial_lambda_at_A, testing::DoubleNear(
      0.0,
      1e-6));
  EXPECT_THAT(b1.partial_lambda_at_X, testing::DoubleNear(
      0.0,
      1e-6));
  EXPECT_THAT(b1.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b1.hot_sites, testing::IsEmpty());
  EXPECT_THAT(b1.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b1.hot_deltas_to_X, testing::IsEmpty());

  // r->x: Site 1 is warm and hot
  const auto& b2 = analysis.branch_infos[2];
  EXPECT_THAT(b2.A, testing::Eq(r));
  EXPECT_THAT(b2.B, testing::Eq(x));
  EXPECT_THAT(b2.is_open, testing::IsFalse());
  EXPECT_THAT(b2.T_to_X, testing::Eq(4.0));
  EXPECT_THAT(b2.partial_lambda_at_A, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b2.partial_lambda_at_X, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rT),
      1e-6));
  EXPECT_THAT(b2.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b2.hot_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b2.hot_muts_to_X, testing::ElementsAre(
      Mutation{rA, 1, rC, -0.5},
      Mutation{rC, 1, rT,  2.0}));
  EXPECT_THAT(b2.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(1, rA, rT)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0
      -mu_JC*2.0 +
      -std::log(P_JC(rA, rA, 2.0)) +
      // Site 1
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -mu_JC*1.0 +
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rT, 4.0)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*2.0 +
      // Site 1
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*0.5 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rA, rC)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*0.5 +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*1.0 +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*1.0 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rC, rT)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rT)*1.0 +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_precarious_path_without_root_test, propose_new_graft) {
  run_propose_new_graft_test(tree);
}

class Spr_move_precarious_path_with_root_test : public Spr_move_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:                 -1.0        0.0        1.0        2.0        3.0
  //
  //                                    +----------+ a (ACNN)  <--- Precarious path - --->a on site 1
  //                                    |          
  //                         +-- A1C ---+ x (ACNN)
  //                         |          |
  //  (ACAA) A2N,A3N-- C1A --+ r (AANN) +C1N------- b (ANNN)
  //                         |
  //                         +A1N------------------ c (ANNN)
  // 
  
  Phylo_tree tree{5};

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Spr_move_precarious_path_with_root_test() {
    tree.root = r;
    tree.ref_sequence = ref_sequence;
  
    tree.at(r).parent = k_no_node;
    tree.at(r).children = {x, c};
    tree.at(r).name = "r";
    tree.at(r).t = -1.0;
    tree.at(r).missations = {{Missation{2, rA}, Missation{3, rA}}, ref_sequence};
    tree.at(r).mutations = {Mutation{rC, 1, rA, -std::numeric_limits<double>::max()}};
  
    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 1, rC, -0.5}};
  
    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = 1.0;
  
    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = 1.0;
    tree.at(b).missations = {{Missation{1, rC}}, ref_sequence};
  
    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = 1.0;
    tree.at(c).missations = {{Missation{1, rA}}, ref_sequence};
  
    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Spr_move_precarious_path_with_root_test, analyze_graft_a) {
  auto can_change_root = true;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(a);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // x->a:  Sites 0 & 1 are warm, and site 0 is hot
  const auto& b0 = analysis.branch_infos[0];
  EXPECT_THAT(b0.A, testing::Eq(x));
  EXPECT_THAT(b0.B, testing::Eq(a));
  EXPECT_THAT(b0.is_open, testing::IsFalse());
  EXPECT_THAT(b0.T_to_X, testing::Eq(1.0));
  EXPECT_THAT(b0.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA),
      1e-6));
  EXPECT_THAT(b0.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA),
      1e-6));
  EXPECT_THAT(estd::ranges::to_vec(b0.warm_sites.slow_elements()),
              testing::IsSupersetOf({0, 1}));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::Contains(0));
  EXPECT_THAT(estd::ranges::to_vec(b0.hot_sites.slow_elements()),
              testing::Not(testing::Contains(1)));
  EXPECT_THAT(b0.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b0.hot_deltas_to_X, testing::IsEmpty());

  // r->x: Site 1 is warm but not hot
  const auto& b1 = analysis.branch_infos[1];
  EXPECT_THAT(b1.A, testing::Eq(r));
  EXPECT_THAT(b1.B, testing::Eq(x));
  EXPECT_THAT(b1.is_open, testing::IsFalse());
  EXPECT_THAT(b1.T_to_X, testing::Eq(2.0));
  EXPECT_THAT(b1.partial_lambda_at_A, testing::Eq(0.0));
  EXPECT_THAT(b1.partial_lambda_at_X, testing::Eq(0.0));
  EXPECT_THAT(b1.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b1.hot_sites, testing::IsEmpty());
  EXPECT_THAT(b1.hot_muts_to_X, testing::IsEmpty());
  EXPECT_THAT(b1.hot_deltas_to_X, testing::IsEmpty());

  // *->r: Site 1 is warm and hot
  const auto& b2 = analysis.branch_infos[2];
  EXPECT_THAT(b2.A, testing::Eq(k_no_node));
  EXPECT_THAT(b2.B, testing::Eq(r));
  EXPECT_THAT(b2.is_open, testing::IsTrue());
  EXPECT_THAT(b2.T_to_X, testing::Eq(2.0));
  EXPECT_THAT(b2.partial_lambda_at_A, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA),
      1e-6));
  EXPECT_THAT(b2.partial_lambda_at_X, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rC),
      1e-6));
  EXPECT_THAT(b2.warm_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b2.hot_sites, testing::ElementsAre(Site_interval{1, 2}));
  EXPECT_THAT(b2.hot_muts_to_X, testing::ElementsAre(
      Mutation{rA, 1, rC, -0.5}));
  EXPECT_THAT(b2.hot_deltas_to_X, testing::IsEmpty());  // open path => no hot deltas

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rA, 1.0)) +
      // Site 1
      -mu_JC*0.5 +
      std::log(mu_JC/3) +
      -mu_JC*0.5 +
      -mu_JC*1.0 +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*1.0 +
      // Site 1
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*0.5 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rA, rC)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*0.5 +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*1.0 +
      std::log(pi_l_a(1, rA) / pi_l_a(1, rC)) +
      // Site 2 (missing)
      0.0 +
      // Site 3 (missing)
      0.0,
      1e-6));
}

TEST_F(Spr_move_precarious_path_with_root_test, propose_new_graft) {
  run_propose_new_graft_test(tree);
}

class Spr_move_tricky_rooty_graft_test : public Spr_move_test_base {
 protected:
  //
  // The tree that we build:
  //
  // Time:                   0.0          1.0        2.0        3.0        4.0
  //
  //                          +C3N,C4N--- A1C ------ A5C --------+ X (ACA--CC-A)
  //                          |
  // (AAACCACTT) T7N-- T8A ---+ P (AAACCAC-A)
  //                          |
  //                          +A1N,A2N--- C3G ------ C6G -------------------+ S (A--GCAG-A)
  //
  // The graft at X has mutations on the P-X, P-S and S-P-X paths, and plenty of crazy root changes
  // 
  
  Phylo_tree tree{3};

  static constexpr Node_index P = 0;
  static constexpr Node_index X = 1;
  static constexpr Node_index S = 2;

  Spr_move_tricky_rooty_graft_test() {
    tree.ref_sequence = ref_sequence = {rA, rA, rA, rC, rC, rA, rC, rT, rT};  // ref_sequence differs from the test base
    evo.partition_for_site = {0, 1, 0, 1, 0, 1, 0, 1, 0};  // 2 site partitions as in the test base
    evo.nu_l = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  
    tree.root = P;
  
    tree.at(P).parent = k_no_node;
    tree.at(P).children = {X, S};
    tree.at(P).name = "P";
    tree.at(P).t = 0.0;
    tree.at(P).missations = {{Missation{7, rT}}, ref_sequence};
    tree.at(P).mutations = {Mutation{rT, 8, rA, -std::numeric_limits<double>::max()}};
  
    tree.at(X).parent = P;
    tree.at(X).children = {};
    tree.at(X).name = "X";
    tree.at(X).t = 3.0;
    tree.at(X).missations = {{Missation{3, rC}, Missation{4, rC}}, ref_sequence};
    tree.at(X).mutations = {Mutation{rA, 1, rC, 1.0}, Mutation{rA, 5, rC, 2.0}};
  
    tree.at(S).parent = P;
    tree.at(S).children = {};
    tree.at(S).name = "S";
    tree.at(S).t = 4.0;
    tree.at(S).missations = {{Missation{1, rA}, Missation{2, rA}}, ref_sequence};
    tree.at(S).mutations = {Mutation{rC, 3, rG, 1.0}, Mutation{rC, 6, rG, 2.0}};
  
    assert_phylo_tree_integrity(tree);
  }
};

TEST_F(Spr_move_tricky_rooty_graft_test, analyze_graft_X) {
  auto can_change_root = true;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(X);

  EXPECT_THAT(analysis.branch_infos, testing::SizeIs(3));

  // P->X: Sites 1,2
  const auto& b_P_X = analysis.branch_infos[Spr_graft::k_branch_info_P_X];
  EXPECT_THAT(b_P_X.A, testing::Eq(P));
  EXPECT_THAT(b_P_X.B, testing::Eq(X));
  EXPECT_THAT(b_P_X.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_X.T_to_X, testing::Eq(3.0));
  EXPECT_THAT(b_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rA) +
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(1)*nu_l(1)*q_l_a(1, rC) +
      mu_l(2)*nu_l(2)*q_l_a(2, rA),
      1e-6));
  EXPECT_THAT(b_P_X.warm_sites, testing::ElementsAre(Site_interval{1, 3}));
  EXPECT_THAT(b_P_X.hot_sites, testing::ElementsAre(Site_interval{1, 3}));
  EXPECT_THAT(b_P_X.hot_muts_to_X, testing::ElementsAre(Mutation{rA, 1, rC, 1.0}));
  EXPECT_THAT(b_P_X.hot_deltas_to_X, testing::IsEmpty());

  // P->S: Sites 3,4
  const auto& b_P_S = analysis.branch_infos[Spr_graft::k_branch_info_P_S];
  EXPECT_THAT(b_P_S.A, testing::Eq(P));
  EXPECT_THAT(b_P_S.B, testing::Eq(S));
  EXPECT_THAT(b_P_S.is_open, testing::IsTrue());
  EXPECT_THAT(b_P_S.T_to_X, testing::Eq(4.0));
  EXPECT_THAT(b_P_S.partial_lambda_at_A, testing::DoubleNear(
      mu_l(3)*nu_l(3)*q_l_a(3, rC) +
      mu_l(4)*nu_l(4)*q_l_a(4, rC),
      1e-6));
  EXPECT_THAT(b_P_S.partial_lambda_at_X, testing::DoubleNear(
      mu_l(3)*nu_l(3)*q_l_a(3, rG) +
      mu_l(4)*nu_l(4)*q_l_a(4, rC),
      1e-6));
  EXPECT_THAT(b_P_S.warm_sites, testing::ElementsAre(Site_interval{3, 5}));
  EXPECT_THAT(b_P_S.hot_sites, testing::ElementsAre(Site_interval{3, 5}));
  EXPECT_THAT(b_P_S.hot_muts_to_X, testing::ElementsAre(Mutation{rC, 3, rG, 1.0}));
  EXPECT_THAT(b_P_S.hot_deltas_to_X, testing::IsEmpty());

  // S->P->X: Sites 0,5,6,8
  const auto& b_S_P_X = analysis.branch_infos[Spr_graft::k_branch_info_S_P_X];
  EXPECT_THAT(b_S_P_X.A, testing::Eq(S));
  EXPECT_THAT(b_S_P_X.B, testing::Eq(P));
  EXPECT_THAT(b_S_P_X.is_open, testing::IsFalse());
  EXPECT_THAT(b_S_P_X.T_to_X, testing::Eq(7.0));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_A, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA) +
      mu_l(5)*nu_l(5)*q_l_a(5, rA) +
      mu_l(6)*nu_l(6)*q_l_a(6, rG) +
      mu_l(8)*nu_l(8)*q_l_a(8, rA),
      1e-6));
  EXPECT_THAT(b_S_P_X.partial_lambda_at_X, testing::DoubleNear(
      mu_l(0)*nu_l(0)*q_l_a(0, rA) +
      mu_l(5)*nu_l(5)*q_l_a(5, rC) +
      mu_l(6)*nu_l(6)*q_l_a(6, rC) +
      mu_l(8)*nu_l(8)*q_l_a(8, rA),
      1e-6));
  EXPECT_THAT(b_S_P_X.warm_sites, testing::ElementsAre(Site_interval{0, 1}, Site_interval{5, 7}, Site_interval{8, 9}));
  EXPECT_THAT(b_S_P_X.hot_sites, testing::ElementsAre(Site_interval{0, 1}, Site_interval{5, 7}, Site_interval{8, 9}));
  EXPECT_THAT(b_S_P_X.hot_muts_to_X, testing::ElementsAre(
      Mutation{rG, 6, rC, -2.0},
      Mutation{rA, 5, rC,  2.0}));
  EXPECT_THAT(b_S_P_X.hot_deltas_to_X, testing::UnorderedElementsAre(
      site_deltas_entry(5, rA, rC),
      site_deltas_entry(6, rG, rC)));

  EXPECT_THAT(analysis.log_alpha_mut, testing::DoubleNear(
      // Site 0 (S->P->X)
      -mu_JC*7.0 +
      -std::log(P_JC(rA, rA, 7.0)) +
      // Site 1 (P->X)
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*2.0 +
      // Site 2 (P->X)
      -mu_JC*3.0 +
      // Site 3 (P->S)
      -mu_JC*1.0 +
      std::log(mu_JC/3) +
      -mu_JC*3.0 +
      // Site 4 (P->X)
      -mu_JC*4.0 +
      // Site 5 (S->P->X)
      -mu_JC*4.0 +
      -mu_JC*2.0 +
      std::log(mu_JC/3) +
      -mu_JC*1.0 +
      -std::log(P_JC(rA, rC, 7.0)) +
      // Site 6 (S->P->X)
      -mu_JC*2.0 +
      std::log(mu_JC/3) +
      -mu_JC*2.0 +
      -mu_JC*3.0 +
      -std::log(P_JC(rG, rC, 7.0)) +
      // Site 7 (missing)
      0.0 +
      // Site 8 (S->P->X)
      -mu_JC*4.0 +
      -mu_JC*3.0 +
      -std::log(P_JC(rA, rA, 7.0)),
      1e-6));
  EXPECT_THAT(analysis.delta_log_G, testing::DoubleNear(
      // Site 0 (P->S part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*4.0 +
      // Site 0 (P->X part)
      -mu_l(0)*nu_l(0)*q_l_a(0, rA)*3.0 +
      // Site 0 (root sequence change)
      0.0 +

      // Site 1 (P->X part)
      -mu_l(1)*nu_l(1)*q_l_a(1, rA)*1.0 +
      std::log(mu_l(1)*nu_l(1)*q_l_ab(1, rA, rC)) +
      -mu_l(1)*nu_l(1)*q_l_a(1, rC)*2.0 +
      // Site 1 (root sequence change)
      std::log(pi_l_a(1, rA) / pi_l_a(1, rC)) +
      
      // Site 2 (P->X part)
      -mu_l(2)*nu_l(2)*q_l_a(2, rA)*3.0 +
      // Site 2 (root sequence change)
      0.0 +
      
      // Site 3 (P->S part)
      -mu_l(3)*nu_l(3)*q_l_a(3, rC)*1.0 +
      std::log(mu_l(3)*nu_l(3)*q_l_ab(3, rC, rG)) +
      -mu_l(3)*nu_l(3)*q_l_a(3, rG)*3.0 +
      // Site 3 (root sequence change)
      std::log(pi_l_a(3, rC) / pi_l_a(3, rG)) +
      
      // Site 4 (P->S part)
      -mu_l(4)*nu_l(4)*q_l_a(4, rC)*4.0 +
      // Site 4 (root sequence change)
      0.0 +
      
      // Site 5 (P->S part)
      -mu_l(5)*nu_l(5)*q_l_a(5, rA)*4.0 +
      // Site 5 (P->X part)
      -mu_l(5)*nu_l(5)*q_l_a(5, rA)*2.0 +
      std::log(mu_l(5)*nu_l(5)*q_l_ab(5, rA, rC)) +
      -mu_l(5)*nu_l(5)*q_l_a(5, rC)*1.0 +
      // Site 5 (root sequence change)
      0.0 +

      // Site 6 (P->S part)
      -mu_l(6)*nu_l(6)*q_l_a(6, rC)*2.0 +
      std::log(mu_l(6)*nu_l(6)*q_l_ab(6, rC, rG)) +
      -mu_l(6)*nu_l(6)*q_l_a(6, rG)*2.0 +
      // Site 6 (P->X part)
      -mu_l(6)*nu_l(6)*q_l_a(6, rC)*3.0 +
      // Site 6 (root sequence change)
      std::log(pi_l_a(6, rC) / pi_l_a(6, rG)) +
      
      // Site 7 (missing)
      0.0 +
      
      // Site 8 (P->S part)
      -mu_l(8)*nu_l(8)*q_l_a(8, rA)*4.0 +
      // Site 8 (P->X part)
      -mu_l(8)*nu_l(8)*q_l_a(8, rA)*3.0 +
      // Site 8 (root sequence change)
      0.0,
      1e-6));
}

TEST_F(Spr_move_tricky_rooty_graft_test, peel_graft_X) {
  auto can_change_root = true;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(X);
  spr.peel_graft(analysis);

  EXPECT_THAT(tree.ref_sequence, testing::Eq(Real_sequence{rA, rA, rA, rC, rC, rA, rC, rT, rT}));
  EXPECT_THAT(tree.ref_sequence, testing::Eq(ref_sequence));
  
  EXPECT_THAT(tree.at(P).mutations, testing::UnorderedElementsAre(
      Mutation{rA, 1, rC, -std::numeric_limits<double>::max()},
      Mutation{rC, 3, rG, -std::numeric_limits<double>::max()},
      Mutation{rC, 6, rG, -std::numeric_limits<double>::max()},
      Mutation{rT, 8, rA, -std::numeric_limits<double>::max()}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(P).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{7, rT}));
  
  EXPECT_THAT(tree.at(X).mutations, testing::UnorderedElementsAre(
      Mutation{rA, 5, rC, 1.5},
      Mutation{rG, 6, rC, 1.5}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(X).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{3, rG}, Missation{4, rC}));
  
  EXPECT_THAT(tree.at(S).mutations, testing::IsEmpty());
  EXPECT_THAT(estd::ranges::to_vec(tree.at(S).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rC}, Missation{2, rA}));
}

TEST_F(Spr_move_tricky_rooty_graft_test, closed_mutations_graft_X) {
  auto can_change_root = true;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(X);
  spr.peel_graft(analysis);

  EXPECT_THAT(spr.count_closed_mutations(analysis), testing::Eq(2));
  EXPECT_THAT(spr.summarize_closed_mutations(analysis), testing::UnorderedElementsAre(
      site_deltas_entry(5, rA, rC),
      site_deltas_entry(6, rG, rC)));
}

TEST_F(Spr_move_tricky_rooty_graft_test, peel_reapply_graft_X) {
  auto can_change_root = true;
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, can_change_root, evo, lambda_i,
    ref_cum_Q_l, num_sites_missing_at_every_node, scratch};
  
  auto analysis = spr.analyze_graft(X);
  spr.peel_graft(analysis);
  spr.apply_graft(analysis);

  EXPECT_THAT(tree.ref_sequence, testing::Eq(Real_sequence{rA, rA, rA, rC, rC, rA, rC, rT, rT}));
  EXPECT_THAT(tree.ref_sequence, testing::Eq(ref_sequence));
  
  EXPECT_THAT(tree.at(P).mutations, testing::UnorderedElementsAre(
      Mutation{rT, 8, rA, -std::numeric_limits<double>::max()}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(P).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{7, rT}));
  
  EXPECT_THAT(tree.at(X).mutations, testing::UnorderedElementsAre(
      Mutation{rA, 1, rC, 1.0},
      Mutation{rA, 5, rC, 2.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(X).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{3, rC}, Missation{4, rC}));
  
  EXPECT_THAT(tree.at(S).mutations, testing::UnorderedElementsAre(
      Mutation{rC, 3, rG, 1.0},
      Mutation{rC, 6, rG, 2.0}));
  EXPECT_THAT(estd::ranges::to_vec(tree.at(S).missations.slow_elements(ref_sequence)),
              testing::ElementsAre(Missation{1, rA}, Missation{2, rA}));
}

auto Spr_move_test_base::run_full_spr_move_test(
    Phylo_tree& tree,
    Node_index X,
    Node_index SS,
    double new_t_P,
    absl::BitGenRef bitgen)
    -> void {
  auto scratch = Scratch_space{};
  auto ref_cum_Q_l = calc_cum_Q_l_for_sequence(tree.ref_sequence, evo);
  auto lambda_i = calc_lambda_i(tree, evo, ref_cum_Q_l);
  auto num_sites_missing_at_every_node = calc_num_sites_missing_at_every_node(tree);
  auto spr = Spr_move{tree, mu_JC, true, evo, lambda_i, ref_cum_Q_l, num_sites_missing_at_every_node, scratch};

  auto old_tree = Phylo_tree{tree};  // Copy

  ASSERT_THAT(X, testing::Ne(tree.root));
  auto P = tree.at(X).parent;
  auto old_t_P = tree.at(P).t;
  auto S = tree.at(P).sibling_of(X);
  auto G = tree.at(P).parent;  // Could be empty
  auto GG = tree.at(SS).parent;  // Could be empty
  if (GG == P) { GG = G; }
  
  // Before move
  assert_phylo_tree_integrity(tree);
  auto old_graft = spr.analyze_graft(X);
  auto old_log_G = calc_log_root_prior(tree, evo) + calc_log_G_below_root(tree, evo);
  
  // Move
  spr.peel_graft(old_graft);
  //std::cerr << "Before: " << tree << "\n";
  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  
  spr.move(X, SS, new_t_P);
  
  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  //std::cerr << "After move: " << tree << "\n";
  
  auto new_graft = spr.propose_new_graft(X, bitgen);
  spr.apply_graft(new_graft);
  //std::cerr << "After regraft: " << tree << "\n";

  // After move
  assert_phylo_tree_integrity(tree);
  auto new_log_G = calc_log_root_prior(tree, evo) + calc_log_G_below_root(tree, evo);
  EXPECT_THAT(new_log_G, testing::DoubleNear(old_log_G - old_graft.delta_log_G + new_graft.delta_log_G, 1e-6));
  EXPECT_THAT(tree.ref_sequence, testing::Eq(old_tree.ref_sequence));
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));

  // Check that tree topology changes were minimal
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(tree.at(node).name);
    if (node == P || node == X || node == S || node == SS) { continue; }  // only some parents can change
    EXPECT_THAT(tree.at(node).parent, testing::Eq(old_tree.at(node).parent));
  }
  
  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(tree.at(node).name);
    if (node == G || node == GG || node == P) { continue; }  // only some children can change
    EXPECT_THAT(tree.at(node).children, testing::UnorderedElementsAreArray(old_tree.at(node).children));
  }

  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(tree.at(node).name);
    
    // All tips and everything except nodes upstream of S or SS (exclusive) should be identical
    EXPECT_THAT(tree.at(node).is_tip(), testing::Eq(old_tree.at(node).is_tip()));
    auto upstream_of_S = node != S && descends_from(tree, S, node);
    auto upstream_of_SS = node != SS && descends_from(tree, SS, node);
    if (tree.at(node).is_tip() || not (upstream_of_S || upstream_of_SS)) {
      auto old_seq = view_of_sequence_at(old_tree, node);
      auto new_seq = view_of_sequence_at(tree, node);
      auto old_missing = reconstruct_missing_sites_at(old_tree, node);
      auto new_missing = reconstruct_missing_sites_at(tree, node);
      EXPECT_THAT(old_missing, testing::Eq(new_missing));
      for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
        SCOPED_TRACE(absl::StrFormat("Site %d", l));
        if (not old_missing.contains(l) && not new_missing.contains(l)) {
          EXPECT_THAT(old_seq[l], testing::Eq(new_seq[l]));
        } else {
          // Missing sites can change as precarious branches grow and shrink
        }
      }
    }
  }

  // Revert everything
  spr.peel_graft(new_graft);
  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  
  spr.move(X, S, old_t_P);
  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  
  //std::cerr << "After move back: " << tree << "\n";
  spr.apply_graft(old_graft);
  assert_phylo_tree_integrity(tree);
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  //std::cerr << "After regraft back: " << tree << "\n";

  for (const auto& node : index_order_traversal(tree)) {
    SCOPED_TRACE(tree.at(node).name);
    EXPECT_THAT(view_of_sequence_at(tree, node), testing::ElementsAreArray(view_of_sequence_at(old_tree, node)));
    EXPECT_THAT(reconstruct_missing_sites_at(tree, node), testing::Eq(reconstruct_missing_sites_at(old_tree, node)));
    EXPECT_THAT(tree.at(node).parent, testing::Eq(old_tree.at(node).parent));
    EXPECT_THAT(tree.at(node).children, testing::UnorderedElementsAreArray(old_tree.at(node).children));
  }
  
  auto old_log_G_redux = calc_log_root_prior(tree, evo) + calc_log_G_below_root(tree, evo);
  EXPECT_THAT(old_log_G_redux, testing::DoubleNear(old_log_G, 1e-6));
  EXPECT_THAT(tree.ref_sequence, testing::Eq(old_tree.ref_sequence));
  EXPECT_THAT(calc_lambda_i(tree, evo, ref_cum_Q_l),
              testing::Pointwise(testing::DoubleNear(1e-6), lambda_i));
  EXPECT_THAT(num_sites_missing_at_every_node, testing::Eq(calc_num_sites_missing_at_every_node(tree)));
}

TEST_F(Spr_move_simple_test, full_spr_move) {
  auto name = [&](Node_index node) { return node == k_no_node ? "-" : tree.at(node).name; };
  struct Case { Node_index X; Node_index SS; double t; };
  for (auto& [X, SS, t] : {
      Case{a, b,  0.5},  // slide
      Case{a, x, -0.5},  // slide
      Case{a, c,  0.5},
      Case{a, r, -1.5},  // slide
      Case{b, a,  0.5},  // slide
      Case{b, x, -0.5},  // slide
      Case{b, c,  0.5},
      Case{b, r, -1.5},  // slide
      Case{c, x, -0.5},  // slide
      Case{c, r, -1.5},  // slide
      Case{c, a,  0.5},  // slide
      Case{c, b,  1.0},  // slide
      Case{x, c, -0.4},  // slide
      Case{x, r, -1.6}   // slide
    }) {
    for (auto seed = 0; seed != 200; ++seed) {
      auto bitgen = std::mt19937_64(seed + 12345);
      SCOPED_TRACE(absl::StrFormat("Seed %d", seed));
      SCOPED_TRACE(absl::StrFormat("%s->%s to %s->%s",
                                   name(tree.at(X).parent), name(X),
                                   name(tree.at(SS).parent), name(SS)));
      run_full_spr_move_test(tree, X, SS, t, bitgen);
    }
  }
}

TEST_F(Spr_move_superfluous_mutation_at_root_test, full_spr_move) {
  auto name = [&](Node_index node) { return node == k_no_node ? "-" : tree.at(node).name; };
  struct Case { Node_index X; Node_index SS; double t; };
  for (auto& [X, SS, t] : {
      Case{x, s, -0.5},     // slide
      Case{x, r, -1.5},     // slide
      Case{s, x, -0.5},     // slide
      Case{s, r, -1.5}}) {  // slide
    for (auto seed = 0; seed != 200; ++seed) {
      auto bitgen = std::mt19937_64(seed + 12345);
      SCOPED_TRACE(absl::StrFormat("Seed %d", seed));
      SCOPED_TRACE(absl::StrFormat("%s->%s to %s->%s",
                                   name(tree.at(X).parent), name(X),
                                   name(tree.at(SS).parent), name(SS)));
      run_full_spr_move_test(tree, X, SS, t, bitgen);
    }
  }
}

TEST_F(Spr_move_superfluous_mutation_not_at_root_test, full_spr_move) {
  auto name = [&](Node_index node) { return node == k_no_node ? "-" : tree.at(node).name; };
  struct Case { Node_index X; Node_index SS; double t; };
  for (auto& [X, SS, t] : {
      Case{a, b,  0.6},  // slide
      Case{a, x, -0.4},  // slide
      Case{a, c,  0.5},
      Case{a, r, -1.5},  // slide
      Case{b, a,  0.6},  // slide
      Case{b, x, -0.4},  // slide
      Case{b, c,  0.5},
      Case{b, r, -1.5},  // slide
      Case{c, x, -0.4},  // slide
      Case{c, r, -1.3},  // slide
      Case{c, a,  0.5},  // slide
      Case{c, b,  1.0},  // slide
      Case{x, c, -0.4},  // slide
      Case{x, r, -1.6}   // slide
    }) {
    for (auto seed = 0; seed != 200; ++seed) {
      auto bitgen = std::mt19937_64(seed + 12345);
      SCOPED_TRACE(absl::StrFormat("Seed %d", seed));
      SCOPED_TRACE(absl::StrFormat("%s->%s to %s->%s",
                                   name(tree.at(X).parent), name(X),
                                   name(tree.at(SS).parent), name(SS)));
      run_full_spr_move_test(tree, X, SS, t, bitgen);
    }
  }
}

TEST_F(Spr_move_precarious_path_without_root_test, full_spr_move) {
  auto name = [&](Node_index node) { return node == k_no_node ? "-" : tree.at(node).name; };
  struct Case { Node_index X; Node_index SS; double t; };
  for (auto& [X, SS, t] : {
      Case{a, b,  1.4},  // slide
      Case{a, y,  0.6},  // slide
      Case{a, c,  1.5},
      Case{a, d,  1.0},
      Case{a, r, -1.5},  // slide
      Case{b, a,  1.4},  // slide
      Case{b, y,  0.6},  // slide
      Case{b, c,  1.5},
      Case{b, d,  1.0},
      Case{b, r, -1.5},  // slide
      Case{c, y,  0.6},  // slide
      Case{c, x, -0.4},  // slide
      Case{c, a,  2.0},  // slide
      Case{c, b,  2.0},  // slide
      Case{c, d,  1.0},
      Case{c, r, -1.5},  // slide
      Case{d, x, -0.4},  // slide
      Case{d, r, -1.6},  // slide
      Case{d, a,  2.0},  // slide
      Case{d, b,  2.0},  // slide
      Case{d, c,  1.5},  // slide
      Case{y, c,  0.4},  // slide
      Case{y, x, -0.6},  // slide
      Case{y, d,  0.3},
      Case{y, r, -1.2},  // slide
      Case{x, d, -0.4},  // slide
      Case{x, r, -1.6}   // slide
    }) {
    for (auto seed = 0; seed != 200; ++seed) {
      auto bitgen = std::mt19937_64(seed + 12345);
      SCOPED_TRACE(absl::StrFormat("Seed %d", seed));
      SCOPED_TRACE(absl::StrFormat("%s->%s to %s->%s",
                                   name(tree.at(X).parent), name(X),
                                   name(tree.at(SS).parent), name(SS)));
      run_full_spr_move_test(tree, X, SS, t, bitgen);
    }
  }
}

TEST_F(Spr_move_precarious_path_with_root_test, full_spr_move) {
  auto name = [&](Node_index node) { return node == k_no_node ? "-" : tree.at(node).name; };
  struct Case { Node_index X; Node_index SS; double t; };
  for (auto& [X, SS, t] : {
      Case{a, b,  0.6},  // slide
      Case{a, x, -0.4},  // slide
      Case{a, c,  0.0},
      Case{a, r, -1.5},  // slide
      Case{b, a,  0.6},  // slide
      Case{b, x, -0.4},  // slide
      Case{b, c,  0.0},
      Case{b, r, -1.5},  // slide
      Case{c, x, -0.6},  // slide
      Case{c, r, -1.4},  // slide
      Case{c, a,  0.5},  // slide
      Case{c, b,  0.5},  // slide
      Case{x, c, -0.6},  // slide
      Case{x, r, -1.2}   // slide
    }) {
    for (auto seed = 0; seed != 200; ++seed) {
      auto bitgen = std::mt19937_64(seed + 12345);
      SCOPED_TRACE(absl::StrFormat("Seed %d", seed));
      SCOPED_TRACE(absl::StrFormat("%s->%s to %s->%s",
                                   name(tree.at(X).parent), name(X),
                                   name(tree.at(SS).parent), name(SS)));
      run_full_spr_move_test(tree, X, SS, t, bitgen);
    }
  }
}

TEST_F(Spr_move_simple_test, sample_mutational_history) {
  // Sites 2 and 3 are A at all points in `tree`.
  // In the target sequence, site 2 is C and site 3 is A.
  // Thus, in all test cases, site 2 will have a delta and site 3 will not.
  // Those for which a mutation and a 1+-step reversal appear on site 3 are deemed `unusual`.
  // An unusual trajectory that's more than just a mutation followed by a reversal is `super-unusual`.
  // For fixed mu & T, we can calculate the probability of getting such unusual and super-unusual trajectories.
  //
  auto scratch = Scratch_space{};

  auto verbose = false;
  auto num_histories = (estd::is_debug_enabled ? 1'000 : 25'000);
  auto L = tree.num_sites();
  
  auto target_start_seq = Real_sequence{rT, rG, rC, rA};
  auto num_unusual_histories = 0;
  auto num_super_unusual_histories = 0;
  for (auto seed = 0; seed != num_histories; ++seed) {
    auto bitgen = std::mt19937_64(seed + 12345);
    SCOPED_TRACE(absl::StrFormat("Seed %d", seed));

    // Pick a random end_loc
    auto branch = r;
    while (branch == r) { branch = absl::Uniform(absl::IntervalClosedOpen, bitgen, 0, std::ssize(tree)); }
    auto t_end = absl::Uniform(absl::IntervalClosedOpen, bitgen,
                                 tree.branch_begin_t(branch), tree.branch_end_t(branch));
    auto end_loc = Phylo_tree_loc{branch, t_end};

    // Pick T so that history will often have no needless mutations, but not always
    auto T_scale = 1.0 / mu_JC;
    auto T = T_scale;  // mu_JC*T == 1; we assume this when assessing p_unusual and p_super_unusual below

    // Calculate deltas
    auto end_seq = target_start_seq;
    for (auto l = 0; l != L; ++l) {
      end_seq[l] = calc_site_state_at(tree, end_loc, l);
    }
    auto deltas = Site_deltas{scratch};
    for (auto l = 0; l != L; ++l) {
      if (target_start_seq[l] != end_seq[l]) {
        deltas.insert({l, {target_start_seq[l], end_seq[l]}});
      }
    }

    // Go!
    auto mut_history = sample_mutational_history(L, T, mu_JC, deltas, scratch, bitgen);
    adjust_mutational_history(mut_history, deltas, tree, end_loc, scratch);

    if (verbose) {
      std::cerr << "From ";
      for (const auto& s : target_start_seq) { std::cerr << to_char(s); }
      std::cerr << " to ";
      for (const auto& s : end_seq) { std::cerr << to_char(s); }
      std::cerr << absl::StreamFormat(", mu_JC=%.4g, T=%.4g, history: Start -> ", mu_JC, T);
      for (const auto& m : mut_history) {
        std::cerr << absl::StreamFormat("%s -> ", absl::FormatStreamed(m));
      }
      std::cerr << " End\n";
    }
    
    // Check timing requirements
    EXPECT_THAT(are_mutations_sorted(mut_history), testing::IsTrue());
    EXPECT_THAT(are_mutation_times_in_range(mut_history, end_loc.t - T, end_loc.t), testing::IsTrue());
    
    // Apply history to target_start_seq and check we reach end_seq
    auto seq = target_start_seq;
    for (const auto& m : mut_history) {
      EXPECT_THAT(seq[m.site], testing::Eq(m.from));
      seq[m.site] = m.to;
    }
    EXPECT_THAT(seq, testing::Eq(end_seq));

    // Was the history `unusual`, as defined above?
    auto num_muts_site_3 = 0;
    for (const auto& m : mut_history) {
      if (m.site == 3) { ++num_muts_site_3; }
    }
    auto is_unusual = num_muts_site_3 > 0;
    auto is_super_unusual = num_muts_site_3 > 2;
    if (is_unusual) { ++num_unusual_histories; }
    if (is_super_unusual) { ++num_super_unusual_histories; }

    if (verbose) {
      if (is_super_unusual) { std::cerr << "  !!!  ^^^  !!!\n"; }
      else if (is_unusual) { std::cerr << "  ^^^\n"; }
    }
  }
  
  auto p_unusual = num_unusual_histories / double(num_histories);
  auto err_p_unusual = std::sqrt(p_unusual * (1-p_unusual) / num_histories);

  auto p_super_unusual = num_super_unusual_histories / double(num_histories);
  auto err_p_super_unusual = std::sqrt(p_super_unusual * (1-p_super_unusual) / num_histories);
  
  // A bit of theory to know what to expect for p_unusual and p_super_unusual.
  //
  //            /  0  1/3 1/3 1/3 \|  
  //            | 1/3  0  1/3 1/3 |
  // Let P_ab = | 1/3 1/3  0  1/3 | be the transition matrix from state a to state b in a jump
  //            \ 1/3 1/3 1/3  0  /|
  // 
  // of the Jukes-Cantor CTMC.  Naturally, (P^m)_{ab} is the probability that a state a changes to
  // a state b after m mutations.
  //
  // Let p_m be the probability that a mutational history starting in state A and of length T
  // has exactly m mutations *and* ends in state A.  After a bit of thinking, we see that:
  //
  //    p_m = (1 / m!) (mu T)^m e^(-mu T) (P^m)_AA.
  //
  // [Sanity check: p_1 == 0]
  //
  // It's hard to derive analytical expressions for p_m, but it's first M entries are readily
  // calculated by the following Python snippet:
  //
  //    import math
  //    import scipy as sp
  //    import numpy as np
  //    muT = 1.0   # Or anything else...
  //    M = 13      # Or anything else...
  //    P = np.array([[0, 1/3, 1/3, 1/3], [1/3, 0, 1/3, 1/3], [1/3, 1/3, 0, 1/3], [1/3, 1/3, 1/3, 0]])
  //    pp = [
  //      1/math.factorial(m) *
  //      muT**m *
  //      math.exp(-muT) *
  //      sp.linalg.fractional_matrix_power(P, m)[0][0]
  //      for m in range(M)]
  //
  // The probability that such a trajectory results for any number of mutations is naturally:
  //
  //     p* = sum_{m=0}^infty p_m.
  //
  // Hence, the conditional probability that an mutational history is unusual given that it is an A~>A history
  // is given by:
  //
  //     p_unusual = (p* - p0) / p*
  //
  // and similarly, the probability that it is super unusual is given by:
  //
  //     p_super_unusual = (p* - p0 - p1 - p2) / p*
  //
  // For mu * T = 1.0, as above, the first few values of p_m output by the above Python snippet are:
  //
  //     p_0 = 0.36787944117144233
  //     p_1 = 0.0
  //     p_2 = 0.061313240195240384
  //     p_3 = 0.013625164487831196
  //     p_4 = 0.0039740063089507656
  //     p_5 = 0.0007569535826572886
  //     p_6 = 0.00012826157928359615
  //     p_7 = 1.8222956619527322e-05
  //     p_8 = 2.2820414997439205e-06
  //     p_9 = 2.534056509973613e-07
  //     p_10 = 2.5345715621097866e-08
  //     p_11 = 2.3039998891494503e-09
  //     p_12 = 1.920043262181459e-10
  //
  // With those values, we get:
  //
  //     p_* = 0.4476978535708957
  //     p_unusual = 0.17828634147519676
  //     p_super_unusual = 0.04133406505439621

  auto expected_p_unusual = 0.17828634147519676;
  auto expected_p_super_unusual = 0.04133406505439621;

  EXPECT_THAT(p_unusual, testing::DoubleNear(expected_p_unusual, 2*err_p_unusual));
  EXPECT_THAT(p_super_unusual, testing::DoubleNear(expected_p_super_unusual, 2*err_p_super_unusual));
}

}  // namespace delphy
