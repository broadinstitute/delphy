#include "gtest/gtest.h"

#include <cmath>
#include <random>

#include "scalable_coalescent.h"

namespace delphy {

TEST(Scalable_coalescent_test, log_prior) {
  auto pop = 20.0;   // N_e * rho, in days
  auto pop_model = Const_pop_model{pop};  // 20.0 = N_e * rho, in days

  auto num_tips = 10;
  auto num_nodes = 2*num_tips - 1;
  auto prior = Scalable_coalescent_prior{pop_model, num_nodes, 0.0, 1.0};
  for (auto i = 0; i != num_tips - 1; ++i) { prior.mark_as_coalescence(i); }
  for (auto i = num_tips - 1; i != num_nodes; ++i) { prior.mark_as_tip(i); }

  EXPECT_NEAR(prior.calc_log_prior(), -(num_tips - 1) * std::log(pop), 1e-8);

  prior.displace_coalescence(0, -50);
  prior.displace_coalescence(1, -40);
  prior.displace_coalescence(2, -40);
  prior.displace_coalescence(3, -30);
  prior.displace_coalescence(4, -30);
  prior.displace_coalescence(5, -30);
  prior.displace_coalescence(6, -30);
  prior.displace_coalescence(7, -20);
  prior.displace_coalescence(8, -10);

  auto expected_log_prior = 0.0
      - (-40.0 - (-50.0)) * (2 * 1) / 2 / pop
      - (-30.0 - (-40.0)) * (4 * 3) / 2 / pop
      - (-20.0 - (-30.0)) * (8 * 7) / 2 / pop
      - (-10.0 - (-20.0)) * (9 * 8) / 2 / pop
      - (-00.0 - (-10.0)) * (10 * 9) / 2 / pop
      - (num_tips - 1) * std::log(pop);
  EXPECT_NEAR(prior.calc_log_prior(), expected_log_prior, 1e-8);

  auto i = num_tips - 1;
  for (const auto& t : {-45.0, -35.0, -30.0, -30.0, -25.0}) {
    prior.displace_tip(i, t);
    ++i;
  }

  expected_log_prior = 0.0                     // ( . ,-50) = 1 lineage
      - (-45.0 - (-50.0)) * (2 * 1) / 2 / pop  // [-50,-45) = 2 lineages
      - (-40.0 - (-45.0)) * (1 * 0) / 2 / pop  // [-45,-40) = 1 lineage
      - (-35.0 - (-40.0)) * (3 * 2) / 2 / pop  // [-40,-35) = 3 lineages
      - (-30.0 - (-35.0)) * (2 * 1) / 2 / pop  // [-35,-30) = 2 lineages
      - (-25.0 - (-30.0)) * (4 * 3) / 2 / pop  // [-30,-25) = 4 lineages
      - (-20.0 - (-25.0)) * (3 * 2) / 2 / pop  // [-25,-20) = 3 lineages
      - (-10.0 - (-20.0)) * (4 * 3) / 2 / pop  // [-20,-10) = 4 lineages
      - (-00.0 - (-10.0)) * (5 * 4) / 2 / pop  // [-10,  0) = 5 lineages
      - (num_tips - 1) * std::log(pop);
  EXPECT_NEAR(prior.calc_log_prior(), expected_log_prior, 1e-8);
}

TEST(Scalable_coalescent_test, delta_log_prior) {
  auto n0 = 20.0;   // N_e * rho, in days
  auto g = 0.1;     // growing by a factor of `e` every 10 days

  auto pop_model = Exp_pop_model{0.0, n0, g};

  auto num_tips = 10;
  auto num_nodes = 2*num_tips - 1;
  auto prior = Scalable_coalescent_prior{pop_model, num_nodes, 0.0, 0.17841};
  for (auto i = 0; i != num_tips - 1; ++i) { prior.mark_as_coalescence(i); }
  for (auto i = num_tips - 1; i != num_nodes; ++i) { prior.mark_as_tip(i); }

  auto ts = std::vector{-50.0, -40.0, -40.0, -30.0, -30.0, -30.0, -30.0, -20.0, -10.0};
  auto prng = std::mt19937{std::random_device{}()};
  std::ranges::shuffle(ts, prng);
  auto log_prior_now = prior.calc_log_prior();
  for (auto i = 0; i != num_tips - 1; ++i) {
    auto delta = prior.calc_delta_log_prior_after_displace_coalescence(i, ts[i]);
    prior.displace_coalescence(i, ts[i]);
    auto new_log_prior = prior.calc_log_prior();
    auto expected_delta = new_log_prior - log_prior_now;

    EXPECT_NEAR(delta, expected_delta, 1e-8);
    log_prior_now = new_log_prior;
  }
  for (auto i = 0; i != num_tips - 1; ++i) {
    prior.displace_coalescence(i, 0.0);
  }

  for (auto i = 0; i != num_tips - 1; ++i) {
    prior.displace_coalescence(i, -60.0);
  }
  auto i = num_tips - 1;
  for (const auto& t : {-45.0, -35.0, -30.0, -30.0, -25.0}) {
    prior.displace_tip(i, t);
    ++i;
  }

  log_prior_now = prior.calc_log_prior();
  for (auto i = 0; i != num_tips - 1; ++i) {
    auto delta = prior.calc_delta_log_prior_after_displace_coalescence(i, ts[i]);
    prior.displace_coalescence(i, ts[i]);
    auto new_log_prior = prior.calc_log_prior();

    if (abs(new_log_prior - (log_prior_now + delta)) > 1e-8) {
      prior.displace_coalescence(i, -60.0);
      delta = prior.calc_delta_log_prior_after_displace_coalescence(i, ts[i]);
      prior.displace_coalescence(i, ts[i]);
      new_log_prior = prior.calc_log_prior();
      EXPECT_NEAR(new_log_prior, log_prior_now + delta, 1e-8);
    }
    log_prior_now = new_log_prior;
  }
}

}  // namespace delphy
