#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <cmath>
#include <random>
#include <numbers>

#include "pop_model.h"

namespace delphy {

TEST(Pop_model_test, const_pop_model_invalid) {
  // Don't support nonpositive populations
  EXPECT_THROW((Const_pop_model{0.0}), std::invalid_argument);
  EXPECT_THROW((Const_pop_model{-1.0}), std::invalid_argument);
}

TEST(Pop_model_test, const_pop_model_normal) {
  auto pop = 20.0;   // N_e * rho, in days

  auto pop_model = Const_pop_model{pop};  // 20.0 = N_e * rho, in days

  EXPECT_EQ(pop_model.pop(), pop);
  EXPECT_EQ(pop_model.pop_at_time(0.0), pop);
  EXPECT_EQ(pop_model.pop_at_time(5.0), pop);
  EXPECT_EQ(pop_model.intensity_integral(2.0, 5.0), (5.0 - 2.0) / pop);
}

TEST(Pop_model_test, const_pop_model_print) {
  auto pop = 20.0;
  auto pop_model = Const_pop_model{pop};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat("Const_pop_model{pop=%g}", pop)));
}

TEST(Pop_model_test, exp_pop_model_invalid) {
  // Don't support nonpositive populations
  EXPECT_THROW((Exp_pop_model{0.0, 0.0, 1.0}), std::invalid_argument);
  EXPECT_THROW((Exp_pop_model{0.0, -1.0, 1.0}), std::invalid_argument);
}

TEST(Pop_model_test, exp_pop_model_normal) {
  auto n0 = 10.0;                // N_e(t=0) * rho, in days
  auto g = std::numbers::ln2;    // e^(ln(2)*t) = 2^t, i.e., doubles every day

  auto pop_model = Exp_pop_model{0.0, n0, g};

  EXPECT_EQ(pop_model.pop_at_t0(), n0);
  EXPECT_EQ(pop_model.growth_rate(), g);

  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(20.0, 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(-2.0), testing::DoubleNear(2.5, 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(-3.0), testing::DoubleNear(1.25, 1e-6));

  EXPECT_THAT(pop_model.pop_integral( 0.0,  1.0), testing::DoubleNear((1/g) * (20.0 - 10.0 ), 1e-6));
  EXPECT_THAT(pop_model.pop_integral(-2.0,  0.0), testing::DoubleNear((1/g) * (10.0 -  2.5 ), 1e-6));
  EXPECT_THAT(pop_model.pop_integral(-3.0, -2.0), testing::DoubleNear((1/g) * ( 2.5 -  1.25), 1e-6));

  EXPECT_THAT(pop_model.intensity_integral( 0.0,  1.0), testing::DoubleNear(
      (1/n0) * (1/g) * (1.0 -  0.5), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral( 1.0,  2.0), testing::DoubleNear(
      (1/n0) * (1/g) * (0.5 - 0.25), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(-3.0, -2.0), testing::DoubleNear(
      (1/n0) * (1/g) * (8.0 -  4.0), 1e-6));
}

TEST(Pop_model_test, exp_pop_model_print) {
  auto t0 = 3.25;
  auto n0 = 10.0;
  auto g = 2.0;

  auto pop_model = Exp_pop_model{t0, n0, g};

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat(
      "Exp_pop_model{t0=%g, n0=%g, g=%g}", t0, n0, g)));
}

TEST(Pop_model_test, skygrid_pop_model_invalid) {
  // No empty ranges
  EXPECT_THAT(([&] { Skygrid_pop_model{{}, {}, Skygrid_pop_model::Type::k_staircase}; }),
              testing::Throws<std::invalid_argument>());

  // x (knot times) and gamma (knot values) must have the same size
  EXPECT_THAT(([&] { Skygrid_pop_model{{1.0, 2.0}, {3.0}, Skygrid_pop_model::Type::k_staircase}; }),
              testing::Throws<std::invalid_argument>());

  // No duplicated time points
  EXPECT_THAT(([&] { Skygrid_pop_model{{1.0, 1.0}, {2.0, 3.0}, Skygrid_pop_model::Type::k_staircase}; }),
              testing::Throws<std::invalid_argument>());
}

TEST(Pop_model_test, skygrid_pop_model_accessors) {
  auto pop_model = Skygrid_pop_model{
    { 1.0, 2.0, 4.0},
    {-4.0, 7.0, 3.0},
    Skygrid_pop_model::Type::k_staircase
  };

  EXPECT_THAT(pop_model.x(), testing::ElementsAre(1.0, 2.0, 4.0));
  EXPECT_THAT(pop_model.gamma(), testing::ElementsAre(-4.0, 7.0, 3.0));
  EXPECT_THAT(pop_model.type(), testing::Eq(Skygrid_pop_model::Type::k_staircase));

  EXPECT_THAT(pop_model.M(), testing::Eq(2));
  EXPECT_THAT(pop_model.T(), testing::Eq(pop_model.x()[2]));
  EXPECT_THAT(pop_model.K(), testing::Eq(pop_model.x()[2] - pop_model.x()[0]));
  EXPECT_THAT(pop_model.x(1), testing::Eq(pop_model.x()[1]));
  EXPECT_THAT(pop_model.gamma(1), testing::Eq(pop_model.gamma()[1]));

  EXPECT_THAT(([&] { pop_model.x(-1); }), testing::Throws<std::out_of_range>());
  EXPECT_THAT(([&] { pop_model.x(5); }), testing::Throws<std::out_of_range>());
  EXPECT_THAT(([&] { pop_model.gamma(-1); }), testing::Throws<std::out_of_range>());
  EXPECT_THAT(([&] { pop_model.gamma(5); }), testing::Throws<std::out_of_range>());

  EXPECT_THAT(pop_model.interval_containing_t(0.0), testing::Eq(0));
  EXPECT_THAT(pop_model.interval_containing_t(1.0), testing::Eq(0));
  EXPECT_THAT(pop_model.interval_containing_t(1.5), testing::Eq(1));
  EXPECT_THAT(pop_model.interval_containing_t(2.0), testing::Eq(1));
  EXPECT_THAT(pop_model.interval_containing_t(2.5), testing::Eq(2));
  EXPECT_THAT(pop_model.interval_containing_t(4.0), testing::Eq(2));
  EXPECT_THAT(pop_model.interval_containing_t(4.5), testing::Eq(3));
}

TEST(Pop_model_test, skygrid_pop_model_normal_staircase) {
  auto pop_model = Skygrid_pop_model{
    { 1.0,  2.0,  4.0,  8.0},
    {-4.0, +7.0, +3.0, +1.0},
    Skygrid_pop_model::Type::k_staircase
  };

  EXPECT_THAT(pop_model.pop_at_time(0.0), testing::DoubleNear(std::exp(-4.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(std::exp(-4.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(1.5), testing::DoubleNear(std::exp(+7.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(2.0), testing::DoubleNear(std::exp(+7.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(3.0), testing::DoubleNear(std::exp(+3.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(4.0), testing::DoubleNear(std::exp(+3.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(5.0), testing::DoubleNear(std::exp(+1.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(6.0), testing::DoubleNear(std::exp(+1.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(7.0), testing::DoubleNear(std::exp(+1.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(8.0), testing::DoubleNear(std::exp(+1.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(9.0), testing::DoubleNear(std::exp(+1.0), 1e-6));

  EXPECT_THAT(pop_model.log_N(0.0), testing::DoubleNear(-4.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(1.0), testing::DoubleNear(-4.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(1.5), testing::DoubleNear(+7.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(2.0), testing::DoubleNear(+7.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(3.0), testing::DoubleNear(+3.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(4.0), testing::DoubleNear(+3.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(5.0), testing::DoubleNear(+1.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(6.0), testing::DoubleNear(+1.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(7.0), testing::DoubleNear(+1.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(8.0), testing::DoubleNear(+1.0, 1e-6));
  EXPECT_THAT(pop_model.log_N(9.0), testing::DoubleNear(+1.0, 1e-6));

  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 0), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 0), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 1), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 1), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 2), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 2), testing::DoubleNear(1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 3), testing::DoubleNear(0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 3), testing::DoubleNear(1, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 3), testing::DoubleNear(1, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 3), testing::DoubleNear(1, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 3), testing::DoubleNear(1, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 0), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 1), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 2), testing::DoubleNear(0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 3), testing::DoubleNear(1, 1e-6));
  
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(0), testing::FieldsAre(
      testing::DoubleNear(-std::numeric_limits<double>::infinity(), 1e-6),
      testing::DoubleNear(1.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(1), testing::FieldsAre(
      testing::DoubleNear(1.0, 1e-6),
      testing::DoubleNear(2.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(2), testing::FieldsAre(
      testing::DoubleNear(2.0, 1e-6),
      testing::DoubleNear(4.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(3), testing::FieldsAre(
      testing::DoubleNear(4.0, 1e-6),
      testing::DoubleNear(+std::numeric_limits<double>::infinity(), 1e-6)));

  // A. Interval wholly before SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(0.0, 0.5), testing::DoubleNear(
      std::exp(-4.0)*(0.5-0.0), 1e-6));

  // B. Interval covering before SkyGrid time range and part of first cell
  EXPECT_THAT(pop_model.pop_integral(0.5, 1.5), testing::DoubleNear(
      std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(1.5-1.0), 1e-6));
  
  // C. Interval covering before SkyGrid time range, first cell and part of second cell
  EXPECT_THAT(pop_model.pop_integral(0.5, 2.5), testing::DoubleNear(
      std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(2.0-1.0) + std::exp(+3.0)*(2.5-2.0), 1e-6));
  
  // D. Interval completely inside a single cell
  EXPECT_THAT(pop_model.pop_integral(2.5, 3.5), testing::DoubleNear(
      std::exp(+3.0)*(3.5-2.5), 1e-6));
  
  // E. Interval straddling two consecutive cells
  EXPECT_THAT(pop_model.pop_integral(1.5, 2.5), testing::DoubleNear(
      std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(2.5-2.0), 1e-6));
  
  // F. Interval straddling three consecutive cells
  EXPECT_THAT(pop_model.pop_integral(1.5, 4.5), testing::DoubleNear(
      std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(4.0-2.0) + std::exp(+1.0)*(4.5-4.0), 1e-6));

  // G. Interval covering part of second-to-last cell, last cell, and after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(2.5, 9.0), testing::DoubleNear(
      std::exp(+3.0)*(4.0-2.5) + std::exp(+1.0)*(8.0-4.0) + std::exp(+1.0)*(9.0-8.0), 1e-6));
  
  // H. Interval covering part of final cell and after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(6.0, 9.0), testing::DoubleNear(
      std::exp(+1.0)*(8.0-6.0) + std::exp(+1.0)*(9.0-8.0), 1e-6));
  
  // I. Interval wholly after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(9.5, 10.5), testing::DoubleNear(
      std::exp(+1.0)*(10.5-9.5), 1e-6));

  
  // Intensity integrals are all the same, but with all the gamma signs flipped
  EXPECT_THAT(pop_model.intensity_integral(0.0, 0.5), testing::DoubleNear(
      std::exp(+4.0)*(0.5-0.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(0.5, 1.5), testing::DoubleNear(
      std::exp(+4.0)*(1.0-0.5) + std::exp(-7.0)*(1.5-1.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(0.5, 2.5), testing::DoubleNear(
      std::exp(+4.0)*(1.0-0.5) + std::exp(-7.0)*(2.0-1.0) + std::exp(-3.0)*(2.5-2.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(2.5, 3.5), testing::DoubleNear(
      std::exp(-3.0)*(3.5-2.5), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(1.5, 2.5), testing::DoubleNear(
      std::exp(-7.0)*(2.0-1.5) + std::exp(-3.0)*(2.5-2.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(1.5, 4.5), testing::DoubleNear(
      std::exp(-7.0)*(2.0-1.5) + std::exp(-3.0)*(4.0-2.0) + std::exp(-1.0)*(4.5-4.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(2.5, 9.0), testing::DoubleNear(
      std::exp(-3.0)*(4.0-2.5) + std::exp(-1.0)*(8.0-4.0) + std::exp(-1.0)*(9.0-8.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(6.0, 9.0), testing::DoubleNear(
      std::exp(-1.0)*(8.0-6.0) + std::exp(-1.0)*(9.0-8.0), 1e-6));
  EXPECT_THAT(pop_model.intensity_integral(9.5, 10.5), testing::DoubleNear(
      std::exp(-1.0)*(10.5-9.5), 1e-6));

  // Derivatives of log-integrals above

  // A. Interval wholly before SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 0), testing::DoubleNear(1.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 1), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 2), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 3), testing::DoubleNear(0.0, 1e-6));
  
  // B. Interval covering before SkyGrid time range and part of first cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 0), testing::DoubleNear(
      (std::exp(-4.0)*(1.0-0.5))
      / (std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(1.5-1.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 1), testing::DoubleNear(
      (std::exp(+7.0)*(1.5-1.0))
      / (std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(1.5-1.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 2), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 3), testing::DoubleNear(0.0, 1e-6));
  
  // C. Interval covering before SkyGrid time range, first cell and part of second cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 0), testing::DoubleNear(
      (std::exp(-4.0)*(1.0-0.5))
      / (std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(2.0-1.0) + std::exp(+3.0)*(2.5-2.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 1), testing::DoubleNear(
      (std::exp(+7.0)*(2.0-1.0))
      / (std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(2.0-1.0) + std::exp(+3.0)*(2.5-2.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 2), testing::DoubleNear(
      (std::exp(+3.0)*(2.5-2.0))
      / (std::exp(-4.0)*(1.0-0.5) + std::exp(+7.0)*(2.0-1.0) + std::exp(+3.0)*(2.5-2.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 3), testing::DoubleNear(0.0, 1e-6));
  
  // D. Interval completely inside a single cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 1), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 2), testing::DoubleNear(1.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 3), testing::DoubleNear(0.0, 1e-6));
  
  // E. Interval straddling two consecutive cells
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 1), testing::DoubleNear(
      (std::exp(+7.0)*(2.0-1.5))
      / (std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(2.5-2.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 2), testing::DoubleNear(
      (std::exp(+3.0)*(2.5-2.0))
      / (std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(2.5-2.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 3), testing::DoubleNear(0.0, 1e-6));
  
  // F. Interval straddling three consecutive cells
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 1), testing::DoubleNear(
      (std::exp(+7.0)*(2.0-1.5))
      / (std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(4.0-2.0) + std::exp(+1.0)*(4.5-4.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 2), testing::DoubleNear(
      (std::exp(+3.0)*(4.0-2.0))
      / (std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(4.0-2.0) + std::exp(+1.0)*(4.5-4.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 3), testing::DoubleNear(
      (std::exp(+1.0)*(4.5-4.0))
      / (std::exp(+7.0)*(2.0-1.5) + std::exp(+3.0)*(4.0-2.0) + std::exp(+1.0)*(4.5-4.0)),
      1e-6));

  // G. Interval covering part of second-to-last cell, last cell, and after SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 1), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 2), testing::DoubleNear(
      (std::exp(+3.0)*(4.0-2.5))
      / (std::exp(+3.0)*(4.0-2.5) + std::exp(+1.0)*(8.0-4.0) + std::exp(+1.0)*(9.0-8.0)),
      1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 3), testing::DoubleNear(
      (std::exp(+1.0)*(8.0-4.0) + std::exp(+1.0)*(9.0-8.0))
      / (std::exp(+3.0)*(4.0-2.5) + std::exp(+1.0)*(8.0-4.0) + std::exp(+1.0)*(9.0-8.0)),
      1e-6));
  
  // H. Interval covering part of final cell and after SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 1), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 2), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 3), testing::DoubleNear(1.0, 1e-6));
  
  // I. Interval wholly after SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(9.5, 10.5, 0), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(9.5, 10.5, 1), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(9.5, 10.5, 2), testing::DoubleNear(0.0, 1e-6));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(9.5, 10.5, 3), testing::DoubleNear(1.0, 1e-6));
}

TEST(Pop_model_test, skygrid_pop_model_normal_log_linear) {
  auto pop_model = Skygrid_pop_model{
    { 1.0,  2.0,  4.0,  8.0},
    {-4.0, +7.0, +3.0, +1.0},
    Skygrid_pop_model::Type::k_log_linear
  };

  EXPECT_THAT(pop_model.pop_at_time(0.0), testing::DoubleNear(std::exp(-4.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(1.0), testing::DoubleNear(std::exp(-4.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(1.5), testing::DoubleNear(std::exp(+1.5), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(2.0), testing::DoubleNear(std::exp(+7.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(3.0), testing::DoubleNear(std::exp(+5.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(4.0), testing::DoubleNear(std::exp(+3.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(5.0), testing::DoubleNear(std::exp(+2.5), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(6.0), testing::DoubleNear(std::exp(+2.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(7.0), testing::DoubleNear(std::exp(+1.5), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(8.0), testing::DoubleNear(std::exp(+1.0), 1e-6));
  EXPECT_THAT(pop_model.pop_at_time(9.0), testing::DoubleNear(std::exp(+1.0), 1e-6));

  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 0), testing::DoubleNear(   1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(0.0, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 0), testing::DoubleNear(   1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.0, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 0), testing::DoubleNear( 0.5, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 1), testing::DoubleNear( 0.5, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(1.5, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 1), testing::DoubleNear(   1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(2.0, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 1), testing::DoubleNear( 0.5, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 2), testing::DoubleNear( 0.5, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(3.0, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 2), testing::DoubleNear(   1, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(4.0, 3), testing::DoubleNear(   0, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 2), testing::DoubleNear(0.75, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(5.0, 3), testing::DoubleNear(0.25, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 2), testing::DoubleNear( 0.5, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(6.0, 3), testing::DoubleNear( 0.5, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 2), testing::DoubleNear(0.25, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(7.0, 3), testing::DoubleNear(0.75, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(8.0, 3), testing::DoubleNear(   1, 1e-6));
  
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 0), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 1), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 2), testing::DoubleNear(   0, 1e-6));
  EXPECT_THAT(pop_model.d_log_N_d_gamma(9.0, 3), testing::DoubleNear(   1, 1e-6));
  
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(0), testing::FieldsAre(
      testing::DoubleNear(-std::numeric_limits<double>::infinity(), 1e-6),
      testing::DoubleNear(2.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(1), testing::FieldsAre(
      testing::DoubleNear(1.0, 1e-6),
      testing::DoubleNear(4.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(2), testing::FieldsAre(
      testing::DoubleNear(2.0, 1e-6),
      testing::DoubleNear(8.0, 1e-6)));
  EXPECT_THAT(pop_model.support_of_d_log_N_d_gamma(3), testing::FieldsAre(
      testing::DoubleNear(4.0, 1e-6),
      testing::DoubleNear(+std::numeric_limits<double>::infinity(), 1e-6)));

  // A. Interval wholly before SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(0.0, 0.5), testing::DoubleNear(
      std::exp(-4.0)*(0.5-0.0),
      1e-6));
  
  // B. Interval covering before SkyGrid time range and part of first cell
  EXPECT_THAT(pop_model.pop_integral(0.5, 1.5), testing::DoubleNear(
      std::exp(-4.0)*(1.0-0.5) +
      1.0 * std::exp(-4.0) * std::exp((+7.0 - -4.0) * 0.0) * std::expm1((+7.0 - -4.0) * (0.5-0.0)) / (+7.0 - -4.0),
      1e-6));
  
  // C. Interval covering before SkyGrid time range, first cell and part of second cell
  EXPECT_THAT(pop_model.pop_integral(0.5, 2.5), testing::DoubleNear(
      std::exp(-4.0)*(1.0-0.5) +
      1.0 * std::exp(-4.0) * std::exp((+7.0 - -4.0) * 0.0) * std::expm1((+7.0 - -4.0) * (1.0 -0.0)) / (+7.0 - -4.0) +
      2.0 * std::exp(+7.0) * std::exp((+3.0 - +7.0) * 0.0) * std::expm1((+3.0 - +7.0) * (0.25-0.0)) / (+3.0 - +7.0),
      1e-6));
  
  // D. Interval completely inside a single cell
  EXPECT_THAT(pop_model.pop_integral(2.5, 3.5), testing::DoubleNear(
      2.0 * std::exp(+7.0) * std::exp((+3.0 - +7.0) * 0.25) * std::expm1((+3.0 - +7.0) * (0.75-0.25)) / (+3.0 - +7.0),
      1e-6));

  // E. Interval straddling two consecutive cells
  EXPECT_THAT(pop_model.pop_integral(1.5, 2.5), testing::DoubleNear(
      1.0 * std::exp(-4.0) * std::exp((+7.0 - -4.0) * 0.5) * std::expm1((+7.0 - -4.0) * (1.0 -0.5)) / (+7.0 - -4.0) +
      2.0 * std::exp(+7.0) * std::exp((+3.0 - +7.0) * 0.0) * std::expm1((+3.0 - +7.0) * (0.25-0.0)) / (+3.0 - +7.0),
      1e-6));
  
  // F. Interval straddling three consecutive cells
  EXPECT_THAT(pop_model.pop_integral(1.5, 4.5), testing::DoubleNear(
      1.0 * std::exp(-4.0) * std::exp((+7.0 - -4.0) * 0.5) * std::expm1((+7.0 - -4.0) * (1.0  -0.5)) / (+7.0 - -4.0) +
      2.0 * std::exp(+7.0) * std::exp((+3.0 - +7.0) * 0.0) * std::expm1((+3.0 - +7.0) * (1.0  -0.0)) / (+3.0 - +7.0) +
      4.0 * std::exp(+3.0) * std::exp((+1.0 - +3.0) * 0.0) * std::expm1((+1.0 - +3.0) * (0.125-0.0)) / (+1.0 - +3.0),
      1e-6));

  // G. Interval covering part of second-to-last cell, last cell, and after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(2.5, 9.0), testing::DoubleNear(
      2.0 * std::exp(+7.0) * std::exp((+3.0 - +7.0) * 0.25) * std::expm1((+3.0 - +7.0) * (1.0  -0.25)) / (+3.0 - +7.0) +
      4.0 * std::exp(+3.0) * std::exp((+1.0 - +3.0) * 0.0 ) * std::expm1((+1.0 - +3.0) * (1.0  -0.0 )) / (+1.0 - +3.0) +
      std::exp(+1.0)*(9.0-8.0),
      1e-6));

  // H. Interval covering part of final cell and after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(6.0, 9.0), testing::DoubleNear(
      4.0 * std::exp(+3.0) * std::exp((+1.0 - +3.0) * 0.5) * std::expm1((+1.0 - +3.0) * (1.0-0.5)) / (+1.0 - +3.0) +
      std::exp(+1.0)*(9.0-8.0),
      1e-6));
  
  // I. Interval wholly after SkyGrid time range
  EXPECT_THAT(pop_model.pop_integral(9.5, 10.5), testing::DoubleNear(
      std::exp(+1.0)*(10.5-9.5),
      1e-6));
  
  
  // *** Intensity integrals are all the same, but with all the gamma signs flipped
  EXPECT_THAT(pop_model.intensity_integral(0.0, 0.5), testing::DoubleNear(
      std::exp(+4.0)*(0.5-0.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(0.5, 1.5), testing::DoubleNear(
      std::exp(+4.0)*(1.0-0.5) +
      1.0 * std::exp(+4.0) * std::exp((-7.0 - +4.0) * 0.0) * std::expm1((-7.0 - +4.0) * (0.5-0.0)) / (-7.0 - +4.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(0.5, 2.5), testing::DoubleNear(
      std::exp(+4.0)*(1.0-0.5) +
      1.0 * std::exp(+4.0) * std::exp((-7.0 - +4.0) * 0.0) * std::expm1((-7.0 - +4.0) * (1.0 -0.0)) / (-7.0 - +4.0) +
      2.0 * std::exp(-7.0) * std::exp((-3.0 - -7.0) * 0.0) * std::expm1((-3.0 - -7.0) * (0.25-0.0)) / (-3.0 - -7.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(2.5, 3.5), testing::DoubleNear(
      2.0 * std::exp(-7.0) * std::exp((-3.0 - -7.0) * 0.25) * std::expm1((-3.0 - -7.0) * (0.75-0.25)) / (-3.0 - -7.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(1.5, 2.5), testing::DoubleNear(
      1.0 * std::exp(+4.0) * std::exp((-7.0 - +4.0) * 0.5) * std::expm1((-7.0 - +4.0) * (1.0 -0.5)) / (-7.0 - +4.0) +
      2.0 * std::exp(-7.0) * std::exp((-3.0 - -7.0) * 0.0) * std::expm1((-3.0 - -7.0) * (0.25-0.0)) / (-3.0 - -7.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(1.5, 4.5), testing::DoubleNear(
      1.0 * std::exp(+4.0) * std::exp((-7.0 - +4.0) * 0.5) * std::expm1((-7.0 - +4.0) * (1.0  -0.5)) / (-7.0 - +4.0) +
      2.0 * std::exp(-7.0) * std::exp((-3.0 - -7.0) * 0.0) * std::expm1((-3.0 - -7.0) * (1.0  -0.0)) / (-3.0 - -7.0) +
      4.0 * std::exp(-3.0) * std::exp((-1.0 - -3.0) * 0.0) * std::expm1((-1.0 - -3.0) * (0.125-0.0)) / (-1.0 - -3.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(2.5, 9.0), testing::DoubleNear(
      2.0 * std::exp(-7.0) * std::exp((-3.0 - -7.0) * 0.25) * std::expm1((-3.0 - -7.0) * (1.0  -0.25)) / (-3.0 - -7.0) +
      4.0 * std::exp(-3.0) * std::exp((-1.0 - -3.0) * 0.0 ) * std::expm1((-1.0 - -3.0) * (1.0  -0.0 )) / (-1.0 - -3.0) +
      std::exp(-1.0)*(9.0-8.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(6.0, 9.0), testing::DoubleNear(
      4.0 * std::exp(-3.0) * std::exp((-1.0 - -3.0) * 0.5) * std::expm1((-1.0 - -3.0) * (1.0-0.5)) / (-1.0 - -3.0) +
      std::exp(-1.0)*(9.0-8.0),
      1e-6));
  EXPECT_THAT(pop_model.intensity_integral(9.5, 10.5), testing::DoubleNear(
      std::exp(-1.0)*(10.5-9.5),
      1e-6));

  // Derivatives of log-integrals above
  //
  // The analytical expressions here are so messy that getting the test right is almost as
  // hard as writing the original calculation.  So instead we independently test the
  // derivative by comparing it to a numerical estimate
  auto numerical_deriv = [&](double a, double b, int k) -> double {
    const auto& old_gamma = pop_model.gamma();
    auto new_gamma = std::vector{old_gamma};

    const auto delta = 1e-4;
    new_gamma[k] += delta;

    auto new_pop_model = Skygrid_pop_model{
      pop_model.x(),
      new_gamma,
      Skygrid_pop_model::Type::k_log_linear
    };

    auto old_log_int = std::log(pop_model.pop_integral(a, b));
    auto new_log_int = std::log(new_pop_model.pop_integral(a, b));
    return (new_log_int - old_log_int) / delta;
  };

  // A. Interval wholly before SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 0), testing::DoubleNear(numerical_deriv(0.0, 0.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 1), testing::DoubleNear(numerical_deriv(0.0, 0.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 2), testing::DoubleNear(numerical_deriv(0.0, 0.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.0, 0.5, 3), testing::DoubleNear(numerical_deriv(0.0, 0.5, 3), 1e-4));
  
  // B. Interval covering before SkyGrid time range and part of first cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 0), testing::DoubleNear(numerical_deriv(0.5, 1.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 1), testing::DoubleNear(numerical_deriv(0.5, 1.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 2), testing::DoubleNear(numerical_deriv(0.5, 1.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 1.5, 3), testing::DoubleNear(numerical_deriv(0.5, 1.5, 3), 1e-4));
  
  // C. Interval covering before SkyGrid time range, first cell and part of second cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 0), testing::DoubleNear(numerical_deriv(0.5, 2.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 1), testing::DoubleNear(numerical_deriv(0.5, 2.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 2), testing::DoubleNear(numerical_deriv(0.5, 2.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(0.5, 2.5, 3), testing::DoubleNear(numerical_deriv(0.5, 2.5, 3), 1e-4));
  
  // D. Interval completely inside a single cell
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 0), testing::DoubleNear(numerical_deriv(2.5, 3.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 1), testing::DoubleNear(numerical_deriv(2.5, 3.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 2), testing::DoubleNear(numerical_deriv(2.5, 3.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 3.5, 3), testing::DoubleNear(numerical_deriv(2.5, 3.5, 3), 1e-4));

  // E. Interval straddling two consecutive cells
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 0), testing::DoubleNear(numerical_deriv(1.5, 2.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 1), testing::DoubleNear(numerical_deriv(1.5, 2.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 2), testing::DoubleNear(numerical_deriv(1.5, 2.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 2.5, 3), testing::DoubleNear(numerical_deriv(1.5, 2.5, 3), 1e-4));

  // F. Interval straddling three consecutive cells
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 0), testing::DoubleNear(numerical_deriv(1.5, 4.5, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 1), testing::DoubleNear(numerical_deriv(1.5, 4.5, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 2), testing::DoubleNear(numerical_deriv(1.5, 4.5, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(1.5, 4.5, 3), testing::DoubleNear(numerical_deriv(1.5, 4.5, 3), 1e-4));

  // G. Interval covering part of second-to-last cell, last cell, and after SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 0), testing::DoubleNear(numerical_deriv(2.5, 9.0, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 1), testing::DoubleNear(numerical_deriv(2.5, 9.0, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 2), testing::DoubleNear(numerical_deriv(2.5, 9.0, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(2.5, 9.0, 3), testing::DoubleNear(numerical_deriv(2.5, 9.0, 3), 1e-4));

  // H. Interval covering part of final cell and after SkyGrid time range
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 0), testing::DoubleNear(numerical_deriv(6.0, 9.0, 0), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 1), testing::DoubleNear(numerical_deriv(6.0, 9.0, 1), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 2), testing::DoubleNear(numerical_deriv(6.0, 9.0, 2), 1e-4));
  EXPECT_THAT(pop_model.d_log_int_N_d_gamma(6.0, 9.0, 3), testing::DoubleNear(numerical_deriv(6.0, 9.0, 3), 1e-4));
}

TEST(Pop_model_test, skygrid_pop_model_staircase_print) {
  auto pop_model = Skygrid_pop_model{
    { 1.0,  2.0,  4.0,  8.0},
    {-4.0, +7.0, +3.0, +1.0},
    Skygrid_pop_model::Type::k_staircase
  };

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat(
      "Skygrid_pop_model{type=k_staircase, ln N(t)=[at t=%g, %g; at t=%g, %g; at t=%g, %g; at t=%g, %g]}",
      pop_model.x(0), pop_model.gamma(0),
      pop_model.x(1), pop_model.gamma(1),
      pop_model.x(2), pop_model.gamma(2),
      pop_model.x(3), pop_model.gamma(3))));
}

TEST(Pop_model_test, skygrid_pop_model_log_linear_print) {
  auto pop_model = Skygrid_pop_model{
    { 1.0,  2.0,  4.0,  8.0},
    {-4.0, +7.0, +3.0, +1.0},
    Skygrid_pop_model::Type::k_log_linear
  };

  auto ss = std::stringstream{};
  ss << pop_model;

  EXPECT_THAT(ss.str(), testing::StrEq(absl::StrFormat(
      "Skygrid_pop_model{type=k_log_linear, ln N(t)=[at t=%g, %g; at t=%g, %g; at t=%g, %g; at t=%g, %g]}",
      pop_model.x(0), pop_model.gamma(0),
      pop_model.x(1), pop_model.gamma(1),
      pop_model.x(2), pop_model.gamma(2),
      pop_model.x(3), pop_model.gamma(3))));
}

TEST(Pop_model_test, extraneous_NaN_in_Skygrid) {
  // Taken from a real run that triggered this bug
  auto pop_model = Skygrid_pop_model{
    {2.000, 8.350, 14.700, 21.050, 27.400, 33.750, 40.100, 46.450, 52.800, 59.150, 65.500, 71.850, 78.200, 84.550, 90.900, 97.250, 103.600, 109.950, 116.300, 122.650, 129.000},
    {5.410, 5.859, 5.952, 5.975, 6.161, 6.789, 7.032, 7.134, 7.261, 8.745, 9.984, 10.012, 9.037, 8.122, 7.214, 8.392, 9.377, 7.492, 7.385, 8.065, 8.067},
    Skygrid_pop_model::Type::k_log_linear
  };

  auto t_min = 129.000;
  auto t_max = 129.868;

  auto N_c = pop_model.pop_integral(t_min, t_max) / (t_max - t_min);
  EXPECT_THAT(N_c, testing::Not(testing::IsNan()));
}

TEST(Pop_model_test, render_population_curve_const_pop_model) {
  auto pop_model = Const_pop_model{20.0};

  auto result = render_population_curve(pop_model, 0.0, 10.0, 10);

  EXPECT_THAT(result, testing::Each(testing::Eq(20.0)));
}

TEST(Pop_model_test, render_population_curve_exp_pop_model) {
  auto t0 = 0.0;
  auto n0 = 10.0;
  auto g = 2.0;

  auto pop_model = Exp_pop_model{t0, n0, g};

  auto result = render_population_curve(pop_model, 0.0, 10.0, 10);

  // int_t^(t+Delta) dt' N(t') = n0/g (exp(g(t+Delta)) - exp(g t))
  //                           = N(t) * (1/g) * (exp(g Delta) - 1)
  for (auto cell = 0; cell != result.num_cells(); ++cell) {
    auto t = cell_lbound(result, cell);
    auto expected = pop_model.pop_at_time(t) * (1/g) * (std::exp(g * result.cell_size()) - 1.0);
    EXPECT_THAT(result.at_cell(cell), testing::DoubleNear(expected, 1e-6));
  }
}

}  // namespace delphy
