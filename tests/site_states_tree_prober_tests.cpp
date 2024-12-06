#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "site_states_tree_prober.h"

namespace delphy {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Site_states_tree_prober_test, invalid_site) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rA, rC, rG};
  tree.root = 0;
  
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};
  
  for (const auto& site : {-1, 5}) {
    EXPECT_THROW((probe_site_states_on_tree(tree, pop_model, site, 0.0, 1.0, 10)), std::out_of_range);
  }
}

TEST(Site_states_tree_prober_test, invalid_timelines) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rA, rC, rG};
  tree.root = 0;
  tree.at_root().t = 0.8;  // > t_start = 0.5
  
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  EXPECT_THROW((probe_site_states_on_tree(tree, pop_model, 1, 0.5, -0.5, 10)), std::invalid_argument);
}

TEST(Site_states_tree_prober_test, trivial) {
  auto tree = Phylo_tree{1};
  tree.ref_sequence = Real_sequence{rA, rC};
  tree.root = 0;
  
  auto pop = 0.2;   // N_e * rho
  auto pop_model = Const_pop_model{pop};

  auto t_start = 0.0;
  auto t_end = 1e-5;
  auto num_cells = 1;

  auto site = Site_index{0};
  auto results = probe_site_states_on_tree(tree, pop_model, site, t_start, t_end, num_cells);

  // If the time domain is empty, results should just reflect the state of the root sequence at that site
  EXPECT_THAT(results[index_of(rA)].at(t_start), testing::DoubleNear(1.0, 1e-6));
  EXPECT_THAT(results[index_of(rC)].at(t_start), testing::DoubleNear(0.0, 1e-6));
}

static auto site_states_tree_prober_test_body(const Pop_model& pop_model) -> void {
  //
  // The tree that we build:
  //
  // Time:    -1.0        0.0        1.0        2.0        3.0
  //          
  //                       +-- T0C -- a
  //                       |
  //            +-- A0T ---+ x
  //            |          |
  //          r +          +-------- T0G ------- b
  //            |
  //            +------------------- A0G ------------------ c
  // 

  auto ref_sequence = Real_sequence{rA, rC};
  auto tree = Phylo_tree{5};

  const auto r = Node_index{0};
  const auto x = Node_index{1};
  const auto a = Node_index{2};
  const auto b = Node_index{3};
  const auto c = Node_index{4};

  tree.root = r;
  tree.ref_sequence = ref_sequence;
  
  tree.at(r).parent = k_no_node;
  tree.at(r).children = {x, c};
  tree.at(r).t_min = -std::numeric_limits<float>::max();
  tree.at(r).t_max = +std::numeric_limits<float>::max();
  tree.at(r).t = -1.0;

  tree.at(x).parent = r;
  tree.at(x).children = {a, b};
  tree.at(x).t_min = -std::numeric_limits<float>::max();
  tree.at(x).t_max = +std::numeric_limits<float>::max();
  tree.at(x).t = 0.0;
  tree.at(x).mutations = {Mutation{rA, 0, rT, 0.5}};

  tree.at(a).parent = x;
  tree.at(a).children = {};
  tree.at(a).t = tree.at(a).t_min = tree.at(a).t_max = 1.0;
  tree.at(a).mutations = {Mutation{rT, 0, rC, 0.5}};

  tree.at(b).parent = x;
  tree.at(b).children = {};
  tree.at(b).t = tree.at(b).t_min = tree.at(b).t_max = 2.0;
  tree.at(b).mutations = {Mutation{rT, 0, rG, 1.0}};

  tree.at(c).parent = r;
  tree.at(c).children = {};
  tree.at(c).t = tree.at(c).t_min = tree.at(c).t_max = 3.0;
  tree.at(c).mutations = {Mutation{rA, 0, rG, 1.0}};

  // We'd like to sample every 0.1 time units from the root to the end
  // Initially fill results with sentinel values to confirm they are all set
  auto t_start = -1.2;
  auto t_end = 3.0;
  auto t_step = 0.2;
  auto num_t_cells = static_cast<int>(std::round((t_end - t_start) / t_step));

  // Do it!
  auto site = Site_index{0};
  auto results = probe_site_states_on_tree(tree, pop_model, site, t_start, t_end, num_t_cells);

  // Check that everything is sensible
  for (auto cell = 0; cell != num_t_cells; ++cell) {
    auto tot_p = 0.0;
    for (const auto& state : k_all_real_seq_letters) {
      auto p = results[index_of(state)].at_cell(cell);
      EXPECT_THAT(p, testing::Ge(-1e6));
      EXPECT_THAT(p, testing::Le(1+1e6));
      tot_p += p;
    }
    EXPECT_THAT(tot_p, testing::DoubleNear(1.0, 1e-6));
  }
  
  // Visual check
  //for (auto state : k_all_real_seq_letters) {
  //  std::cout << to_char(state) << ": ";
  //  for (auto cell = 0; cell != num_t_cells; ++cell) {
  //    std::cout << absl::StreamFormat("%.1f, ", results[index_of(state)].at_cell(cell));
  //  }
  //  std::cout << "\n";
  //}
}

TEST(Site_states_tree_prober_test, typical_const_pop_model) {
  site_states_tree_prober_test_body(Const_pop_model{1.0});
}

TEST(Site_states_tree_prober_test, typical_exp_pop_model) {
  site_states_tree_prober_test_body(Exp_pop_model{0.0, 1.0, 0.1});
}

}  // namespace delphy
