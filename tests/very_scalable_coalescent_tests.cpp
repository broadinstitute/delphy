#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <cmath>
#include <random>

#include "very_scalable_coalescent.h"

namespace delphy {

TEST(Very_scalable_coalescent_test, cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;
  
  EXPECT_THAT(cell_for(7.5, t_ref, t_step), testing::Eq(0));
  EXPECT_THAT(cell_for(6.1, t_ref, t_step), testing::Eq(0));
  EXPECT_THAT(cell_for(5.9, t_ref, t_step), testing::Eq(1));
  EXPECT_THAT(cell_for(2.5, t_ref, t_step), testing::Eq(2));
  EXPECT_THAT(cell_for(3.5, t_ref, t_step), testing::Eq(2));

  EXPECT_THAT(cell_lbound(0, t_ref, t_step), testing::DoubleNear(6.0, 1e-6));
  EXPECT_THAT(cell_ubound(0, t_ref, t_step), testing::DoubleNear(8.0, 1e-6));
  EXPECT_THAT(cell_lbound(1, t_ref, t_step), testing::DoubleNear(4.0, 1e-6));
  EXPECT_THAT(cell_ubound(1, t_ref, t_step), testing::DoubleNear(6.0, 1e-6));
}

TEST(Very_scalable_coalescent_test, add_interval_single_whole_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(2.0, 4.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  2.0,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_whole_cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(2.0, 6.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  2.0,  2.0,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_single_partial_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.0, 3.5, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.5, 0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_partial_cells) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.5, 5.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  1.0,  0.5,  0.0));
}

TEST(Very_scalable_coalescent_test, add_interval_two_partial_cells_and_one_whole_cell) {
  const auto t_ref = 8.0;
  const auto t_step = 2.0;

  using namespace Very_scalable_coalescent_prior;

  auto cells = std::vector<double>(4, 0.0);

  //                                     (6,8] (4,6] (2,4] (0,2]
  EXPECT_THAT(cells, testing::ElementsAre(0.0,  0.0,  0.0,  0.0));
  add_interval(3.5, 7.0, +2.0, cells, t_ref, t_step);
  EXPECT_THAT(cells, testing::ElementsAre(1.0,  2.0,  0.5,  0.0));
}

TEST(Very_scalable_coalescent_test, delta_log_prior) {
  auto n0 = 20.0;   // N_e * rho, in days
  auto g = 0.1;     // growing by a factor of `e` every 10 days

  auto pop_model = std::make_shared<Exp_pop_model>(0.0, n0, g);

  auto num_tips = 10;
  auto num_nodes = 2*num_tips - 1;

  // Kludge: Very_scalable_coalescent_prior_part doesn't actually use the subtree topology,
  // just the node times and whether or not they are inner nodes
  auto tree = Phylo_tree{num_nodes};
  auto t_ref = 0.0;
  for (auto i = 0; i != num_tips - 1; ++i) {
    tree.at(i).children = {1,2};  // "Inner node"
    tree.at(i).t = t_ref;
  }
  for (auto i = num_tips - 1; i != num_nodes; ++i) {
    tree.at(i).children = {};     // "Tip"
    tree.at(i).t = t_ref;
  }
  
  auto bitgen = std::mt19937(12345);
  auto t_step = 0.17841;
  auto prior = Very_scalable_coalescent_prior_part{
    pop_model, tree, bitgen, true, t_ref, t_step,
    {8.5, 4.7},  // k_bar_p
    {-0.5672, +6.141897},  // k_twiddle_bar_p
    {+1.8963, -3.527041},  // k_twiddle_bar
    {pop_model->pop_integral(t_ref - 1*t_step, t_ref - 0*t_step),
     pop_model->pop_integral(t_ref - 2*t_step, t_ref - 1*t_step)
    },  // popsize_bar
    {4, 2}   // num_active_parts
  };

  auto ts = std::vector{-50.0, -40.0, -40.0, -30.0, -30.0, -30.0, -30.0, -20.0, -10.0};
  auto prng = std::mt19937{std::random_device{}()};
  std::ranges::shuffle(ts, prng);
  auto log_prior_now = prior.calc_partial_log_prior();
  for (auto i = 0; i != num_tips - 1; ++i) {
    auto delta = prior.calc_delta_partial_log_prior_after_displace_coalescence(t_ref, ts[i]);
    tree.at(i).t = ts[i];
    prior.coalescence_displaced(t_ref, ts[i]);
    auto new_log_prior = prior.calc_partial_log_prior();
    auto expected_delta = new_log_prior - log_prior_now;

    EXPECT_NEAR(delta, expected_delta, 1e-8);
    log_prior_now = new_log_prior;
  }
  for (auto i = 0; i != num_tips - 1; ++i) {
    prior.coalescence_displaced(ts[i], t_ref);
    tree.at(i).t = t_ref;
  }

  for (auto i = 0; i != num_tips - 1; ++i) {
    prior.coalescence_displaced(t_ref, -60.0);
    tree.at(i).t = t_ref;
  }
  auto i = num_tips - 1;
  for (const auto& t : {-45.0, -35.0, -30.0, -30.0, -25.0}) {
    tree.at(i).t = t;
    prior.tip_displaced(t_ref, t);
    ++i;
  }

  log_prior_now = prior.calc_partial_log_prior();
  for (auto i = 0; i != num_tips - 1; ++i) {
    auto delta = prior.calc_delta_partial_log_prior_after_displace_coalescence(t_ref, ts[i]);
    tree.at(i).t = ts[i];
    prior.coalescence_displaced(t_ref, ts[i]);
    auto new_log_prior = prior.calc_partial_log_prior();

    if (abs(new_log_prior - (log_prior_now + delta)) > 1e-8) {
      prior.coalescence_displaced(ts[i], -60.0);
      tree.at(i).t = -60.0;
      delta = prior.calc_delta_partial_log_prior_after_displace_coalescence(-60.0, ts[i]);
      prior.coalescence_displaced(-60.0, ts[i]);
      tree.at(i).t = ts[i];
      new_log_prior = prior.calc_partial_log_prior();
      EXPECT_NEAR(new_log_prior, log_prior_now + delta, 1e-8);
    }
    log_prior_now = new_log_prior;
  }
}

}  // namespace delphy
