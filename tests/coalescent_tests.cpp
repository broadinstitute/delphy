#include "gtest/gtest.h"

#include <cmath>
#include <random>

#include "coalescent.h"

namespace delphy {

TEST(Coalescent_test, counting) {
  auto pop = 20.0;   // N_e * rho, in days
  auto pop_model = Const_pop_model{pop};  // 20.0 = N_e * rho, in days

  auto num_tips = 10;
  auto prior = Coalescent_prior{&pop_model, num_tips};

  EXPECT_EQ(prior.count_tips_between(0.0, 0.0), 0);
  EXPECT_EQ(prior.count_tips_between(-1.0, 1.0), num_tips);
  EXPECT_EQ(prior.count_tips_between(-1.0, 0.0), 0);
  EXPECT_EQ(prior.count_tips_between(0.0, 1.0), num_tips);

  EXPECT_EQ(prior.count_coalescences_between(0.0, 0.0), 0);
  EXPECT_EQ(prior.count_coalescences_between(-1.0, 1.0), num_tips - 1);
  EXPECT_EQ(prior.count_coalescences_between(-1.0, 0.0), 0);
  EXPECT_EQ(prior.count_coalescences_between(0.0, 1.0), num_tips - 1);

  prior.displace_tip(0.0, 5.0);
  prior.displace_tip(0.0, 3.0);
  prior.displace_coalescence(0.0, 2.0);

  EXPECT_EQ(prior.count_tips_between(2.0, 6.0), 2);
  EXPECT_EQ(prior.count_tips_between(2.0, 4.0), 1);
  EXPECT_EQ(prior.count_coalescences_between(1.0, 3.0), 1);
  EXPECT_EQ(prior.count_coalescences_between(2.5, 3.5), 0);

  EXPECT_ANY_THROW(prior.displace_tip(55.0, 6.0));
  EXPECT_ANY_THROW(prior.displace_coalescence(75.0, 42.0));
}

TEST(Coalescent_test, log_prior) {
  auto pop = 20.0;   // N_e * rho, in days
  auto pop_model = Const_pop_model{pop};  // 20.0 = N_e * rho, in days

  auto num_tips = 10;
  auto prior = Coalescent_prior{&pop_model, num_tips};

  EXPECT_NEAR(prior.calc_log_prior(), -(static_cast<double>(num_tips) - 1) * std::log(pop), 1e-8);

  for (auto t : {-50, -40, -40, -30, -30, -30, -30, -20, -10}) {
    prior.displace_coalescence(0.0, t);
  }

  auto expected_log_prior = 0
      - (-40.0 - (-50.0)) * (2 * 1) / 2 / pop
      - (-30.0 - (-40.0)) * (4 * 3) / 2 / pop
      - (-20.0 - (-30.0)) * (8 * 7) / 2 / pop
      - (-10.0 - (-20.0)) * (9 * 8) / 2 / pop
      - (-00.0 - (-10.0)) * (10 * 9) / 2 / pop
      - (num_tips - 1) * std::log(pop);
  EXPECT_NEAR(prior.calc_log_prior(), expected_log_prior, 1e-8);

  for (auto t : {-45, -35, -30, -30, -25}) {
    prior.displace_tip(0.0, t);
  }

  expected_log_prior = 0                           // ( . ,-50) = 1 lineage
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

TEST(Coalescent_test, delta_log_prior) {
  auto n0 = 20.0;   // N_e * rho, in days
  auto g = 0.1;     // growing by a factor of `e` every 10 days

  auto pop_model = Exp_pop_model{0.0, n0, g};

  auto num_tips = 10;
  auto prior = Coalescent_prior{&pop_model, num_tips};

  auto ts = std::vector{-50.0, -40.0, -40.0, -30.0, -30.0, -30.0, -30.0, -20.0, -10.0};
  auto prng = std::mt19937{std::random_device{}()};
  std::ranges::shuffle(ts, prng);
  auto log_prior_now = prior.calc_log_prior();
  for (auto t : ts) {
    auto delta = prior.calc_delta_log_prior_after_displace_coalescence(0.0, t);
    prior.displace_coalescence(0.0, t);
    auto new_log_prior = prior.calc_log_prior();

    EXPECT_NEAR(new_log_prior, log_prior_now + delta, 1e-8);
    log_prior_now = new_log_prior;
  }
  for (auto t : ts) {
    prior.displace_coalescence(t, 0.0);
  }

  for (auto i = 0; i != num_tips - 1; ++i) {
    prior.displace_coalescence(0.0, -60.0);
  }
  for (auto t : {-45.0, -35.0, -30.0, -30.0, -25.0}) {
    prior.displace_tip(0.0, t);
  }

  log_prior_now = prior.calc_log_prior();
  for (auto t : ts) {
    auto delta = prior.calc_delta_log_prior_after_displace_coalescence(-60.0, t);
    prior.displace_coalescence(-60.0, t);
    auto new_log_prior = prior.calc_log_prior();

    if (abs(new_log_prior - (log_prior_now + delta)) > 1e-8) {
      prior.displace_coalescence(t, -60.0);
      delta = prior.calc_delta_log_prior_after_displace_coalescence(-60.0, t);
      prior.displace_coalescence(-60.0, t);
      new_log_prior = prior.calc_log_prior();
      EXPECT_NEAR(new_log_prior, log_prior_now + delta, 1e-8);
    }
    log_prior_now = new_log_prior;
  }
}

}  // namespace delphy
