#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "spr_study.h"
#include "phylo_tree_calc.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

class Spr_study_test : public testing::Test {
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

  static constexpr Node_index r = 0;
  static constexpr Node_index x = 1;
  static constexpr Node_index a = 2;
  static constexpr Node_index b = 3;
  static constexpr Node_index c = 4;

  Spr_study_test() {
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

TEST_F(Spr_study_test, trivial_builder) {
  // No work => empty result
  auto scratch = Scratch_space{};
  auto builder = Spr_study_builder{tree, a, 5.0, Scratch_interval_set{scratch}, scratch};

  EXPECT_THAT(builder.result, testing::IsEmpty());
}

TEST_F(Spr_study_test, full_spr_study_a) {
  auto scratch = Scratch_space{};

  // Start with a CANN attached at x
  auto site_deltas = calc_site_deltas_between(tree, x, a, scratch);
  auto missing_at_X = reconstruct_missing_sites_at(tree, a, scratch);

  auto t_X = 1.5;  // A little after 1.0 to include parts of lower branches
  
  auto builder = Spr_study_builder{tree, a, t_X, missing_at_X, scratch};
  builder.seed_fill_from(b, 0, std::move(site_deltas), true);

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = b, .mut_idx = 1, .t_min = -0.5, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 2},   // Note t_max != 2.0
      Candidate_region{.branch = b, .mut_idx = 0, .t_min = -1.0, .t_max = -0.5, .min_muts = 1},
      Candidate_region{.branch = r, .mut_idx = 1, .t_min = -std::numeric_limits<double>::max(), .t_max = -1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 0, .t_min = -1.0, .t_max =  0.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 1, .t_min =  0.0, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 1})); // Note t_max != 3.0
}

TEST_F(Spr_study_test, full_spr_study_a_no_root) {
  auto scratch = Scratch_space{};

  // Start with a CANN attached at x
  auto site_deltas = calc_site_deltas_between(tree, x, a, scratch);
  auto missing_at_X = reconstruct_missing_sites_at(tree, a, scratch);

  auto t_X = 1.5;  // A little after 1.0 to include parts of lower branches
  
  auto builder = Spr_study_builder{tree, a, t_X, missing_at_X, scratch};
  builder.seed_fill_from(b, 0, std::move(site_deltas), false);  // Don't allow regrafting above the root

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = b, .mut_idx = 1, .t_min = -0.5, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 2},   // Note t_max != 2.0
      Candidate_region{.branch = b, .mut_idx = 0, .t_min = -1.0, .t_max = -0.5, .min_muts = 1},
      //Candidate_region{.branch = r, .mut_idx = 1, .t_min = -std::numeric_limits<double>::max(), .t_max = -1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 0, .t_min = -1.0, .t_max =  0.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 1, .t_min =  0.0, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 1})); // Note t_max != 3.0
}

TEST_F(Spr_study_test, study_a_up_to_1_mut_away) {
  auto scratch = Scratch_space{};

  // Start with a CANN attached at x
  auto site_deltas = calc_site_deltas_between(tree, x, a, scratch);
  auto missing_at_X = reconstruct_missing_sites_at(tree, a, scratch);

  auto t_X = 1.5;  // A little after 1.0 to include parts of lower branches
  
  auto builder = Spr_study_builder{tree, a, t_X, missing_at_X, scratch};

  builder.max_muts_from_start = 1;
  
  builder.seed_fill_from(b, 0, std::move(site_deltas), true);

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = b, .mut_idx = 1, .t_min = -0.5, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 2},   // Note t_max != 2.0
      Candidate_region{.branch = b, .mut_idx = 0, .t_min = -1.0, .t_max = -0.5, .min_muts = 1},
      Candidate_region{.branch = r, .mut_idx = 1, .t_min = -std::numeric_limits<double>::max(), .t_max = -1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 0, .t_min = -1.0, .t_max =  0.0, .min_muts = 1}));
      //Candidate_region{.branch = c, .mut_idx = 1, .t_min =  0.0, .t_max =  1.0, .min_muts = 1},
      //Candidate_region{.branch = c, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 1})); // Note t_max != 3.0
}

TEST_F(Spr_study_test, full_spr_study_x) {
  auto scratch = Scratch_space{};

  // Start with a TANN attached upstream of c
  auto site_deltas = calc_site_deltas_between(tree, c, x, scratch);
  auto missing_at_X = reconstruct_missing_sites_at(tree, x, scratch);

  auto builder = Spr_study_builder{tree, x, tree.at(x).t, missing_at_X, scratch};
  builder.seed_fill_from(c, std::ssize(tree.at(c).mutations), std::move(site_deltas), true);

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = c, .mut_idx = 3, .t_min = -std::numeric_limits<double>::max(), .t_max = 0.0, .min_muts = 1}));
}

TEST_F(Spr_study_test, full_spr_study_c) {
  auto scratch = Scratch_space{};

  // Start with a GNAN attached upstream of x
  auto site_deltas = calc_site_deltas_between(tree, x, c, scratch);
  auto missing_at_X = reconstruct_missing_sites_at(tree, c, scratch);

  auto builder = Spr_study_builder{tree, c, tree.at(c).t, missing_at_X, scratch};
  builder.seed_fill_from(x, std::ssize(tree.at(x).mutations), std::move(site_deltas), true);

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = x, .mut_idx = 2, .t_min = -std::numeric_limits<double>::max(), .t_max = 0.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 0, .t_min =  0.0, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 1, .t_min =  1.0, .t_max =  2.0, .min_muts = 1},
      Candidate_region{.branch = a, .mut_idx = 0, .t_min =  0.0, .t_max =  0.5, .min_muts = 1},
      Candidate_region{.branch = a, .mut_idx = 1, .t_min =  0.5, .t_max =  1.0, .min_muts = 1}));
}

TEST_F(Spr_study_test, full_study_new_seq) {
  auto scratch = Scratch_space{};

  // Where should we attach a TAAN?  Start above the root (AAAN)
  auto site_deltas = Site_deltas{{
      site_deltas_entry(0, rA, rT)},
    scratch};
  auto missing_at_X = Scratch_interval_set{{Site_interval{3, 4}}, scratch};

  auto t_X = 1.5;  // A little after 1.0 to include parts of lower branches
  
  auto builder = Spr_study_builder{tree, k_no_node, t_X, missing_at_X, scratch};

  builder.seed_fill_from(tree.root, std::ssize(tree.at_root().mutations), std::move(site_deltas), true);

  EXPECT_THAT(builder.result, testing::UnorderedElementsAre(
      Candidate_region{.branch = a, .mut_idx = 0, .t_min =  0.0, .t_max =  0.5, .min_muts = 0},
      Candidate_region{.branch = a, .mut_idx = 1, .t_min =  0.5, .t_max =  1.0, .min_muts = 1},
      Candidate_region{.branch = b, .mut_idx = 0, .t_min =  0.0, .t_max =  1.0, .min_muts = 0},
      Candidate_region{.branch = b, .mut_idx = 1, .t_min =  1.0, .t_max =  1.5, .min_muts = 1},   // Note t_max != 2.0
      Candidate_region{.branch = x, .mut_idx = 1, .t_min = -0.5, .t_max =  0.0, .min_muts = 0},
      Candidate_region{.branch = x, .mut_idx = 0, .t_min = -1.0, .t_max = -0.5, .min_muts = 1},
      Candidate_region{.branch = r, .mut_idx = 1, .t_min = -std::numeric_limits<double>::max(), .t_max = -1.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 0, .t_min = -1.0, .t_max =  0.0, .min_muts = 1},
      Candidate_region{.branch = c, .mut_idx = 1, .t_min =  0.0, .t_max =  1.0, .min_muts = 0},
      Candidate_region{.branch = c, .mut_idx = 2, .t_min =  1.0, .t_max =  1.5, .min_muts = 1})); // Note t_max != 3.0
}

}  // namespace delphy
