#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "site_deltas.h"

#include "phylo_tree_calc.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Site_deltas_complex_test : public testing::Test {
 protected:
  //
  // The tree that we build:
  //
  // Time:             -1.0          0.0        1.0        2.0        3.0
  //                   
  //                                  +-- T0C -- a (CAAA)
  //                                  |
  //                     +---- A0T ---+ x (TAAA)
  //                     |            |
  // (AAAA) -----------+ r (AAAA)     +-------- A1G ------- b (TGAA)
  //                     |
  //                     +----------A0T-------- T0G ------------------ c (GAAA)
  // 

  Real_sequence ref_sequence{rA, rA, rA, rA};
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
  
  Site_deltas_complex_test() {
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
    tree.at(r).t_min = -std::numeric_limits<float>::max();
    tree.at(r).t_max = +std::numeric_limits<float>::max();
    tree.at(r).t = -1.0;

    tree.at(x).parent = r;
    tree.at(x).children = {a, b};
    tree.at(x).name = "x";
    tree.at(x).t_min = -std::numeric_limits<float>::max();
    tree.at(x).t_max = +std::numeric_limits<float>::max();
    tree.at(x).t = 0.0;
    tree.at(x).mutations = {Mutation{rA, 0, rT, -0.5}};
    
    tree.at(a).parent = x;
    tree.at(a).children = {};
    tree.at(a).name = "a";
    tree.at(a).t = tree.at(a).t_min = tree.at(a).t_max = 1.0;
    tree.at(a).mutations = {Mutation{rT, 0, rC, 0.5}};

    tree.at(b).parent = x;
    tree.at(b).children = {};
    tree.at(b).name = "b";
    tree.at(b).t = tree.at(a).t_min = tree.at(a).t_max = 2.0;
    tree.at(b).mutations = {Mutation{rA, 1, rG, 1.0}};

    tree.at(c).parent = r;
    tree.at(c).children = {};
    tree.at(c).name = "c";
    tree.at(c).t = tree.at(a).t_min = tree.at(a).t_max = 3.0;
    tree.at(c).mutations = {Mutation{rA, 0, rT, 0.0}, Mutation{rT, 0, rG, 1.0}};

    assert_phylo_tree_integrity(tree);
  }
};

TEST(Site_deltas_test_test, push_pop_front_site_deltas) {
  auto scratch = Scratch_space{};

  // Slowly displace the start of site_deltas from end to start and back:
  //
  //                                              site_deltas (initial)
  //                                                  |
  //              A0T    T1A     A2G     T0A    A1C  /-\               |
  //  start +------x------x-------x-------x------x-----+ end
  //  (ATA)  \____(1)____(2)_____(3)_____(4)____(5)___/  (ACG)
  //                              |
  //                          site_deltas (final)
  
  auto site_deltas = Site_deltas{scratch};

  EXPECT_THAT(site_deltas, testing::IsEmpty());

  push_front_site_deltas(Seq_delta{1, rA, rC}, site_deltas);  // (5)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rA, rC)));

  push_front_site_deltas(Seq_delta{0, rT, rA}, site_deltas);  // (4)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rA, rC)));

  push_front_site_deltas(Seq_delta{2, rA, rG}, site_deltas);  // (3)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rA, rC),
      site_deltas_entry(2, rA, rG)));

  push_front_site_deltas(Seq_delta{1, rT, rA}, site_deltas);  // (2)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rT, rC),
      site_deltas_entry(2, rA, rG)));

  push_front_site_deltas(Seq_delta{0, rA, rT}, site_deltas);  // (1)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rT, rC),
      site_deltas_entry(2, rA, rG)));

  pop_front_site_deltas(Seq_delta{0, rA, rT}, site_deltas);  // rev (1)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rT, rC),
      site_deltas_entry(2, rA, rG)));

  pop_front_site_deltas(Seq_delta{1, rT, rA}, site_deltas);  // rev (2)
  
  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rA, rC),
      site_deltas_entry(2, rA, rG)));

  pop_front_site_deltas(Seq_delta{2, rA, rG}, site_deltas);  // rev (3)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rA, rC)));
  
  pop_front_site_deltas(Seq_delta{0, rT, rA}, site_deltas);  // rev (4)
  
  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rA, rC)));

  pop_front_site_deltas(Seq_delta{1, rA, rC}, site_deltas);  // rev (5)

  EXPECT_THAT(site_deltas, testing::IsEmpty());
}

TEST(Site_deltas_test, push_pop_back_site_deltas) {
  auto scratch = Scratch_space{};

  // Slowly displace the start of site_deltas from start to end and back:
  //
  //      site_deltas (initial)
  //          |
  //         /-\  A0T    T1A     A2G     T0A    A1C
  //  start +------x------x-------x-------x------x-----+ end
  //  (ATA)  \____(1)____(2)_____(3)_____(4)____(5)___/  (ACG)
  //                              |
  //                          site_deltas (final)
  
  auto site_deltas = Site_deltas{scratch};

  EXPECT_THAT(site_deltas, testing::IsEmpty());

  push_back_site_deltas(Seq_delta{0, rA, rT}, site_deltas);  // (1)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT)));

  push_back_site_deltas(Seq_delta{1, rT, rA}, site_deltas);  // (2)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT),
      site_deltas_entry(1, rT, rA)));

  push_back_site_deltas(Seq_delta{2, rA, rG}, site_deltas);  // (3)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT),
      site_deltas_entry(1, rT, rA),
      site_deltas_entry(2, rA, rG)));

  push_back_site_deltas(Seq_delta{0, rT, rA}, site_deltas);  // (4)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rT, rA),
      site_deltas_entry(2, rA, rG)));

  push_back_site_deltas(Seq_delta{1, rA, rC}, site_deltas);  // (5)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rT, rC),
      site_deltas_entry(2, rA, rG)));

  pop_back_site_deltas(Seq_delta{1, rA, rC}, site_deltas);  // rev (5)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(1, rT, rA),
      site_deltas_entry(2, rA, rG)));

  pop_back_site_deltas(Seq_delta{0, rT, rA}, site_deltas);  // rev (4)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT),
      site_deltas_entry(1, rT, rA),
      site_deltas_entry(2, rA, rG)));

  pop_back_site_deltas(Seq_delta{2, rA, rG}, site_deltas);  // rev (3)

  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT),
      site_deltas_entry(1, rT, rA)));
  
  pop_back_site_deltas(Seq_delta{1, rT, rA}, site_deltas);  // rev (2)
  
  EXPECT_THAT(site_deltas, testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT)));

  pop_back_site_deltas(Seq_delta{0, rA, rT}, site_deltas);  // rev (1)

  EXPECT_THAT(site_deltas, testing::IsEmpty());
}

auto apply_site_deltas(const Site_deltas& site_deltas, const Real_sequence& seq) -> Real_sequence {
  auto result = Real_sequence{seq};  // Copy
  for (const auto& [l, delta] : site_deltas) {
    CHECK_NE(delta.from, delta.to) << l;
    CHECK_EQ(delta.from, result[l]);
    result[l] = delta.to;
  }
  return result;
}

TEST_F(Site_deltas_complex_test, displace_site_delta_starts) {
  // A few manual sanity checks, and then a few property tests
  auto scratch = Scratch_space{};
  
  for (const auto& i : {a,b,c,x}) {
    auto site_deltas = Site_deltas{scratch};
    
    displace_site_deltas_start_upwards(tree, site_deltas, tree.node_loc(i), tree.node_loc(r));
    
    EXPECT_THAT(apply_site_deltas(site_deltas, tree.ref_sequence),
                testing::ElementsAreArray(view_of_sequence_at(tree, i)));
    
    displace_site_deltas_start_downwards(tree, site_deltas, tree.node_loc(r), tree.node_loc(i), scratch);
    
    EXPECT_THAT(site_deltas, testing::IsEmpty());
  }
}

TEST(Site_deltas_test, append_site_deltas) {
  auto scratch = Scratch_space{};
  
  auto site_deltas_1 = Site_deltas{scratch};
  site_deltas_1.insert({0, {rA, rT}});
  site_deltas_1.insert({2, {rC, rG}});
  
  auto site_deltas_2 = Site_deltas{scratch};
  site_deltas_2.insert({0, {rT, rA}});
  site_deltas_2.insert({1, {rG, rA}});
  site_deltas_2.insert({2, {rG, rT}});

  append_site_deltas(site_deltas_1, std::as_const(site_deltas_2));

  EXPECT_THAT(site_deltas_1, testing::UnorderedElementsAre(
      site_deltas_entry(1, rG, rA),
      site_deltas_entry(2, rC, rT)));
}

TEST_F(Site_deltas_complex_test, calc_site_deltas_between) {
  auto scratch = Scratch_space{};

  // a -> *
  
  EXPECT_THAT(calc_site_deltas_between(tree, a, a, scratch), testing::IsEmpty());
  EXPECT_THAT(calc_site_deltas_between(tree, a, b, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rC, rT),
      site_deltas_entry(1, rA, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, a, c, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rC, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, a, x, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rC, rT)));
  EXPECT_THAT(calc_site_deltas_between(tree, a, r, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rC, rA)));
  
  // b -> *

  EXPECT_THAT(calc_site_deltas_between(tree, b, a, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rC),
      site_deltas_entry(1, rG, rA)));
  EXPECT_THAT(calc_site_deltas_between(tree, b, b, scratch), testing::IsEmpty());
  EXPECT_THAT(calc_site_deltas_between(tree, b, c, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rG),
      site_deltas_entry(1, rG, rA)));
  EXPECT_THAT(calc_site_deltas_between(tree, b, x, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(1, rG, rA)));
  EXPECT_THAT(calc_site_deltas_between(tree, b, r, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA),
      site_deltas_entry(1, rG, rA)));
  
  // c -> *

  EXPECT_THAT(calc_site_deltas_between(tree, c, a, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rG, rC)));
  EXPECT_THAT(calc_site_deltas_between(tree, c, b, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rG, rT),
      site_deltas_entry(1, rA, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, c, c, scratch), testing::IsEmpty());
  EXPECT_THAT(calc_site_deltas_between(tree, c, x, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rG, rT)));
  EXPECT_THAT(calc_site_deltas_between(tree, c, r, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rG, rA)));
  
  // x -> *

  EXPECT_THAT(calc_site_deltas_between(tree, x, a, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rC)));
  EXPECT_THAT(calc_site_deltas_between(tree, x, b, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(1, rA, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, x, c, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, x, x, scratch), testing::IsEmpty());
  EXPECT_THAT(calc_site_deltas_between(tree, x, r, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rT, rA)));
  
  // r -> *

  EXPECT_THAT(calc_site_deltas_between(tree, r, a, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rC)));
  EXPECT_THAT(calc_site_deltas_between(tree, r, b, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT),
      site_deltas_entry(1, rA, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, r, c, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rG)));
  EXPECT_THAT(calc_site_deltas_between(tree, r, x, scratch), testing::UnorderedElementsAre(
      site_deltas_entry(0, rA, rT)));
  EXPECT_THAT(calc_site_deltas_between(tree, r, r, scratch), testing::IsEmpty());

  // A few tricky tree locs

  EXPECT_THAT(calc_site_deltas_between(tree, Phylo_tree_loc{r, -2.0}, Phylo_tree_loc{c, 0.5}, scratch),
              testing::UnorderedElementsAre(
                  site_deltas_entry(0, rA, rT)));

  EXPECT_THAT(calc_site_deltas_between(tree, Phylo_tree_loc{c, 0.5}, Phylo_tree_loc{x, 0.0}, scratch),
              testing::IsEmpty());
  
  EXPECT_THAT(calc_site_deltas_between(tree, Phylo_tree_loc{b, 1.5}, Phylo_tree_loc{c, 0.5}, scratch),
              testing::UnorderedElementsAre(
                  site_deltas_entry(1, rG, rA)));
}

}  // namespace delphy
